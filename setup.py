from distutils.core import setup

reqs = [
    "astropy",
    "numpy",
    "scipy",
    "matplotlib",
    "mwa_pb_lookup",
    "calplots",
    "casacore",
    "healpy",
    "requests",
    "mysqlconnector",
    "psutil",
]

scripts = [
    "gleam_x/bin/pyhead.py",
    "gleam_x/bin/alt_az_corrector.py",
    "gleam_x/bin/aocal_diff.py",
    "gleam_x/bin/aocal_phaseref.py",
    "gleam_x/bin/beam_value_at_radec.py",
    "gleam_x/bin/calc_pointing.py",
    "gleam_x/bin/clip_clean_components.py",
    "gleam_x/bin/crop_catalogue.py",
    "gleam_x/bin/dd_flux_mod.py",
    "gleam_x/bin/generate_beam_list.py",
    "gleam_x/bin/generate_weight_map.py",
    "gleam_x/bin/iono_update.py",
    "gleam_x/bin/ms_flag_by_uvdist.py",
    "gleam_x/bin/multiply.py",
    "gleam_x/bin/new_fk5_template.py",
    "gleam_x/bin/polyfit_snapshots.py",
    "gleam_x/bin/psf_combine_axes.py",
    "gleam_x/bin/psf_create.py",
    "gleam_x/bin/psf_projected.py",
    "gleam_x/bin/psf_select.py",
    "gleam_x/bin/pyhead.py",
    "gleam_x/bin/track_task.py",
    "gleam_x/bin/vo2model.py",
    "gleam_x/bin/fits_trim.py",
    "gleam_x/bin/check_assign_solutions.py",
    "gleam_x/db/check_sources_vs_obsids.py",
    "gleam_x/db/check_src_fov.py",
    "gleam_x/db/import_observations_from_db.py",
    "gleam_x/utils/download_obsid_list.py",
    "gleam_x/utils/flag_tiles_bad_dipoles.py",
]

setup(
    name="gleam_x",
    version="0.1",
    author="Natasha Hurley-Walker, Paul Hancock, Tim Galvin",
    description="Python scripts to support the processing of GLEAM-X data.",
    url="https://github.com/nhurleywalker/GLEAM-X-pipeline",
    long_description=open("README.md").read(),
    packages=["gleam_x", "gleam_x.bin", "gleam_x.db", "gleam_x.utils"],
    requires=reqs,
    scripts=scripts,
)
