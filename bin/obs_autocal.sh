#! /bin/bash

# This script was orginally written to use slurm job arrays 
# Author: Tim Galvin
# Suggested change is to use mulit-prog, contributions by Pascal Elahi
# Extra comments are to guide why certain changes have been 
# made 

# set -x

usage()
{
# echo "obs_autocal.sh [-d dep] [-a account] [-t] obsnum
#   -p project : project, no default
#   -d dep     : job number for dependency (afterok)
#   -i         : disable the ionospheric metric tests (default = False)
#   -t         : test. Don't submit job, just make the batch file
#                and then return the submission command
#   obsnum     : the obsid to process, or a text file of obsids (newline separated). 
#                A job-array task will be submitted to process the collection of obsids. " 1>&2;
# exit 1;
echo "obs_autocal.sh [-d dep] [-a account] [-t] obsnum
  -p project : project, no default
  -d dep     : job number for dependency (afterok)
  -i         : disable the ionospheric metric tests (default = False)
  -t         : test. Don't submit job, just make the batch file
               and then return the submission command
  file_of_obs: text file of obsids (newline separated). 
               A multi-prog slurm job script will be produced. 
               " 1>&2;
exit 1;
}

pipeuser=${GXUSER}

dep=
tst=
ion=1

# parse args and set options
while getopts ':tia:d:p:' OPTION
do
    case "$OPTION" in
	d)
	    dep=${OPTARG}
	    ;;
    a)
        account=${OPTARG}
        ;;
	p)
	    project=${OPTARG}
	    ;;
	i)
	    ion=
	    ;;
	t)
	    tst=1
	    ;;
	? | : | h)
	    usage
	    ;;
  esac
done
# set the obsid to be the first non option
shift  "$(($OPTIND -1))"
# PJE: changing name of obsnum variable to obslist
# obsnum=$1
obslist=$1

# if obsid or project are empty then just print help
if [[ -z ${obslist} || -z ${project} ]]
then
    usage
fi

if [[ ! -z ${GXACCOUNT} ]]
then
    account="--account=${GXACCOUNT}"
fi

# # Establish job array options
# if [[ -f ${obsnum} ]]
# then
#     numfiles=$(wc -l "${obsnum}" | awk '{print $1}')
#     jobarray="--array=1-${numfiles}"
# else
#     numfiles=1
#     jobarray=''
# fi

# if obs is not text file, exit
if [ ! -f ${obslist} ]
then 
    echo "Was not passed a text file of observation numbers. "
    usage
fi

numfiles=$(wc -l "${obslist}" | awk '{print $1}')

queue="-p ${GXSTANDARDQ}"
datadir="${GXSCRATCH}/$project"

# set dependency
# if [[ ! -z ${dep} ]]
# then
#     if [[ -f ${obsnum} ]]
#     then
#         depend="--dependency=aftercorr:${dep}"
#     else
#         depend="--dependency=afterok:${dep}"
#     fi
# fi
# PJE: Not clear how to best set dependences
# PJE: but if mult-prog grouping is acceptable
# can simply add this to the batch script
if [[ ! -z ${dep} ]]
then
    depend="--dependency=aftercorr:${dep}"
fi

# PJE: Now here are major changes as idea is 
# to generate a script wrapping the commands for 
# each ob and then construct a multi-prog config
# script listing each of these scripts 
# note that for testing I have used $MYGROUP 
# env variable to write the sbatch and conf files 
# This could be changed to GXSCRIPT directory 
script_dir=$MYGROUP
conffile=${script_dir}/multi_prog.${project}.conf
batchscript=${script_dir}/run_multi_prog.${project}.sbatch
echo "# Config file for ${project}" > ${conffile}
echo "#!/bin/bash" > ${batchfile}
for taskid in $(seq ${numfiles})
do
    curobsnum=$(head -n "${taskid}" "${obslist}" | tail -1)
    script="${GXSCRIPT}/autocal_${curobsnum}.sh"
    stdoutfile="${GXLOG}/autocal_${curobsnum}.stdout"
    stderfile="${GXLOG}/autocal_${curobsnum}.stderr"
    cat "${GXBASE}/templates/autocal.tmpl" | sed -e "s:OBSNUM:${curobsnum}:g" \
                                     -e "s:DATADIR:${datadir}:g" \
                                     -e "s:IONOTEST:${ion}:g" \
                                     -e "s:PIPEUSER:${pipeuser}:g > ${stdoutfile} " > "${script}"
    chmod 755 "${script}"
    singularity_script="${GXSCRIPT}/singularity.autocal_${curobsnum}.sh" 
    # sbatch submissions need to start with a shebang
    echo '#!/bin/bash' > ${singularity_script}.sh
    echo "singularity run ${GXCONTAINER} ${script}" >> ${singularity_script}.sh
    chmod 755 ${singularity_script}

    master_thread_id=$((${GXNCPULINE}*${taskid}))
    echo "${master_thread_id} ${singularity_script} ${runargs}" >> ${conffile}
done 
# submission 
totalmem=$((${GXABSMEMORY}*${numfiles}))
totalcpus=$((${GXNCPULINE}*${numfiles}))
# here adding requirements to single batch file 
echo "#SBATCH --begin=now+5minutes " >> ${batchscript}
echo "#SBATCH --export=ALL " >> ${batchscript}
echo "#SBATCH --time=02:00:00 " >> ${batchscript}
echo "#SBATCH --mem=${totalmem}G " >> ${batchscript}
echo "#SBATCH -M ${GXCOMPUTER} " >> ${batchscript}
echo "#SBATCH --output=${output} " >> ${batchscript}
echo "#SBATCH --error=${error} " >> ${batchscript}
echo "#SBATCH --ntasks=${totalcpus} " >> ${batchscript}
echo "#SBATCH ${account} " >> ${batchscript}
echo "#SBATCH ${GXTASKLINE} " >> ${batchscript}
echo "#SBATCH ${depend} " >> ${batchscript}
echo "#SBATCH ${queue} " >> ${batchscript}

# here following original code and constructing a submission string
sub="sbatch "
sub+="--begin=now+5minutes "
sub+="--export=ALL "
sub+="--time=02:00:00 "
sub+="--mem=${totalmem}G "
sub+="-M ${GXCOMPUTER} "
sub+="--output=${output} "
sub+="--error=${error} "
sub+="--ntasks=${totalcpus} "
sub+="${account} "
sub+="${GXTASKLINE} "
sub+="${depend} "
sub+="${queue} "
sub+="${batchfile} "
    
echo ${sub}

# PJE: OLD CODE below

# script="${GXSCRIPT}/autocal_${obsnum}.sh"

# cat "${GXBASE}/templates/autocal.tmpl" | sed -e "s:OBSNUM:${obsnum}:g" \
#                                      -e "s:DATADIR:${datadir}:g" \
#                                      -e "s:IONOTEST:${ion}:g" \
#                                      -e "s:PIPEUSER:${pipeuser}:g" > "${script}"


# output="${GXLOG}/autocal_${obsnum}.o%A"
# error="${GXLOG}/autocal_${obsnum}.e%A"

# if [[ -f ${obsnum} ]]
# then
#    output="${output}_%a"
#    error="${error}_%a"
# fi

# chmod 755 "${script}"

# # sbatch submissions need to start with a shebang
# echo '#!/bin/bash' > ${script}.sbatch
# echo "singularity run ${GXCONTAINER} ${script}" >> ${script}.sbatch

# sub="sbatch --begin=now+5minutes --export=ALL  --time=02:00:00 --mem=${GXABSMEMORY}G -M ${GXCOMPUTER} --output=${output} --error=${error}"
# sub="${sub} ${GXNCPULINE} ${account} ${GXTASKLINE} ${jobarray} ${depend} ${queue} ${script}.sbatch"
# if [[ ! -z ${tst} ]]
# then
#     echo "script is ${script}"
#     echo "submit via:"
#     echo "${sub}"
#     exit 0
# fi

# # submit job
# jobid=($(${sub}))
# jobid=${jobid[3]}

# echo "Submitted ${script} as ${jobid} . Follow progress here:"

# for taskid in $(seq ${numfiles})
#     do
#     # rename the err/output files as we now know the jobid
#     obserror=$(echo "${error}" | sed -e "s/%A/${jobid}/" -e "s/%a/${taskid}/")
#     obsoutput=$(echo "${output}" | sed -e "s/%A/${jobid}/" -e "s/%a/${taskid}/")

#     if [[ -f ${obsnum} ]]
#     then
#         obs=$(sed -n -e "${taskid}"p "${obsnum}")
#     else
#         obs=$obsnum
#     fi

#     if [ "${GXTRACK}" = "track" ]
#     then
#         # record submission
#         ${GXCONTAINER} track_task.py queue --jobid="${jobid}" --taskid="${taskid}" --task='calibrate' --submission_time="$(date +%s)" --batch_file="${script}" \
#                             --obs_id="${obs}" --stderr="${obserror}" --stdout="${obsoutput}"
#     fi

#     echo "$obsoutput"
#     echo "$obserror"
# done
