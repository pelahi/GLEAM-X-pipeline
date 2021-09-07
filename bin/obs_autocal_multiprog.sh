#! /bin/bash

# This script was orginally written to use slurm job arrays 
# Author: Tim Galvin
# Suggested change is to use mulit-prog, contributions by Pascal Elahi
# Extra comments are to guide why certain changes have been 
# made 

# set -x

usage()
{
echo "obs_autocal.sh [-d dep] [-a account] [-t] obsnum
  -p project  : project, no default
  -s profile  : full path to gleam profile script file to source settings,
                no default
  -d dep      : job number for dependency (afterok)
  -i          : disable the ionospheric metric tests (default = False)
  -t          : test. Don't submit job, just make the batch file
                and then return the submission command
  -g          : for debugging scripts, none submitted 
  file_of_obs : text file of obsids (newline separated). 
                A multi-prog slurm job script will be produced. 
               " 1>&2;
exit 1;
}

pipeuser=${GXUSER}
dep=
tst=
ion=1
debugging_scripts=0

# parse args and set options
while getopts ':tiag:d:p:s:' OPTION
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
	s)
	    profile_script=${OPTARG}
	    ;;
	i)
	    ion=
	    ;;
	t)
	    tst=1
	    ;;
	g)
        debugging_scripts=1
        ;;
	? | : | h)
	    usage
	    ;;
  esac
done
# set the obsid to be the first non option
shift  "$(($OPTIND -1))"
obslist=$1

# if obsid or project are empty then just print help
if [[ -z ${obslist} || -z ${project} || -z ${profile_script} ]]
then
    usage
fi

if [[ ! -z ${GXACCOUNT} ]]
then
    account="--account=${GXACCOUNT}"
fi

# if obs is not text file, exit
if [ ! -f ${obslist} ]
then 
    echo "Was not passed a text file of observation numbers. "
    usage
fi

numfiles=$(wc -l "${obslist}" | awk '{print $1}')

queue="-p ${GXSTANDARDQ}"
datadir="${GXSCRATCH}/$project"

# PJE: Not clear how to best set dependences
# PJE: but if mult-prog grouping is acceptable
# can simply add this to the batch script
if [[ ! -z ${dep} ]]
then
    depend="--dependency=aftercorr:${dep}"
fi

# PJE: Now here are major changes as idea is 
# Idea to generate a script wrapping the commands for 
# each ob and then construct a multi-prog config
# script listing each of these scripts 
script_dir=$GXSCRIPT/
conffile=${script_dir}/multi_prog.${project}.conf
batchscript=${script_dir}/run_multi_prog.${project}.sbatch
echo "# Config file for ${project}" > ${conffile}
echo "#!/bin/bash" > ${batchscript}
for taskid in $(seq ${numfiles})
do
    curobsnum=$(head -n "${taskid}" "${obslist}" | tail -1)
    script="${GXSCRIPT}/autocal_multiprog_${curobsnum}.sh"
    stdoutfile="${GXLOG}/autocal_multiprog_${curobsnum}.stdout"
    stderrfile="${GXLOG}/autocal_multiprog_${curobsnum}.stderr"
    # set the OMP_NUM_THREADS  
    cat "${GXBASE}/templates/autocal_multiprog.tmpl" | sed -e "s:OBSNUM:${curobsnum}:g" \
                                     -e "s:TASKID:${taskid}:g" \
                                     -e "s:DATADIR:${datadir}:g" \
                                     -e "s:IONOTEST:${ion}:g" \
                                     -e "s:PIPEUSER:${pipeuser}:g " >> "${script}"
    chmod 755 "${script}"
    singularity_script="${GXSCRIPT}/singularity.autocal_multiprog_${curobsnum}.sh" 
    # sbatch submissions need to start with a shebang
    echo '#!/bin/bash' > ${singularity_script}
    echo "singularity run ${GXCONTAINER} ${script} 1>${stdoutfile} 2>${stderrfile} " >> ${singularity_script}
    chmod 755 ${singularity_script}

    master_thread_id=$((${GXNPCPUS}*(${taskid}-1)))
    echo "${master_thread_id} ${singularity_script} ${runargs}" >> ${conffile}
    idle_thread_start=$((${master_thread_id}+1))
    idle_thread_end=$((${master_thread_id}+${GXNPCPUS}-1))
    echo "${idle_thread_start}-${idle_thread_end} echo \"task %d running as a OMP/pthread\"" >> ${conffile}
done 
# submission 
totalmem=$((${GXABSMEMORY}*${numfiles}))
totalcpus=$((${GXNPCPUS}*${numfiles}))
output=${GXLOG}/slurm-${project}.o%A
error=${GXLOG}/slurm-${project}.e%A
# here adding requirements to single batch file 
echo "#SBATCH --begin=now+5minutes " >> ${batchscript}
echo "#SBATCH --export=ALL " >> ${batchscript}
echo "#SBATCH --time=02:00:00 " >> ${batchscript}
#echo "#SBATCH --mem=${totalmem}GB " >> ${batchscript}
echo "#SBATCH --mem=${GXABSMEMORY}GB " >> ${batchscript}
echo "#SBATCH -M ${GXCOMPUTER} " >> ${batchscript}
echo "#SBATCH --output=${output} " >> ${batchscript}
echo "#SBATCH --error=${error} " >> ${batchscript}
echo "#SBATCH --ntasks=${totalcpus} " >> ${batchscript}
echo "#SBATCH ${account} " >> ${batchscript}
echo "#SBATCH ${queue} " >> ${batchscript}
if [ ! -z ${depend} ]
then 
    echo "#SBATCH ${depend} " >> ${batchscript}
fi 
echo "# load the singularity module " >> ${batchscript}
echo "module singularity " >> ${batchscript}
echo "source ${profile_script}" >> ${batchscript}
if [ ${debugging_scripts} -eq 1 ]
then
    echo "echo \"Environment \"" >> ${batchscript}
    echo "env" >> ${batchscript}
    
fi
echo "srun --multi-prog ${conffile}" >> ${batchscript}

# for debubing purposes
if [ ${debugging_scripts} -eq 1 ]
then
    echo "Just testing script production" 
    exit 1;
fi

#submit script 
sub="sbatch ${batchscript}"
echo ${sub}
jobid=($(${sub}))
jobid=${jobid[3]}

echo "Submitted ${script} as ${jobid} . Follow progress here:"
for taskid in $(seq ${numfiles})
do
    curobsnum=$(head -n "${taskid}" "${obslist}" | tail -1)
    stdoutfile="${GXLOG}/autocal_${curobsnum}.stdout"
    stderrfile="${GXLOG}/autocal_${curobsnum}.stderr"
    if [ "${GXTRACK}" = "track" ]
    then
        # record submission
        ${GXCONTAINER} track_task.py queue --jobid="${jobid}" \
            --taskid="${taskid}" \
            --task='calibrate' \
            --submission_time="$(date +%s)" \
            --batch_file="${script}" \
            --obs_id="${obs}" \
            --stderr="${stderrfile}" \
            --stdout="${stdoutfile}"
    fi
done
