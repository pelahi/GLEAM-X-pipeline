#! /usr/bin/env python

import logging
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.stats.circstats import circmean
from astropy.coordinates import SkyCoord
from argparse import ArgumentParser
from typing import Iterable, Union
from gleam_x.utils.obsid_ops import read_obsids_file

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


def filter_files_from_obsid(
    files: Iterable[str], obsids: Iterable[int]
) -> Iterable[str]:
    """Given a set of file paths (assuming GLEAM-X directories) and a set of obsids, filter
    out all paths that do not have an obsid beloning to the set of obsids. 

    Args:
        files (Iterable[str]): List of files to filter out
        obsids (Iterable[int]): List of acceptable obsids

    Returns:
        Iterable[str]: Path to files whose obsid is included in the obsid set
    """
    filter_files = []
    for f in files:
        for obsid in obsids:
            if str(obsid) in f:
                filter_files.append(f)
                break

    return filter_files


def calculate_mean(
    files: Iterable[str],
    print_pos: bool = False,
    filter_obsids: str = None,
    ra_field: str = "CRVAL1",
    dec_field: str = "CRVAL2",
) -> Union[float, float]:
    """Calculate the mean position based on RA/Dec values extracted from a set of FITS files

    Args:
        files (Iterable[str]): FITS files to process
        print_pos (bool, optional): Output results to stdout. Defaults to False.
        filter_obsids (str, optional): Path to standard GLEAm-X obsids text file (new line delimited). Defaults to None.
        ra_field (str, optional): RA field in FITS header. Defaults to "CRVAL1".
        dec_field (str, optional): Dec field in FITS header. Defaults to "CRVAL2".

    Returns:
        Union[float, float]: Mean RA and Dec position
    """

    if filter_obsids is not None:
        logger.debug(f"Number of input files: {len(files)}")

        obsids = read_obsids_file(filter_obsids)
        logger.debug(f"Number of obsids to include: {len(obsids)}")

        files = filter_files_from_obsid(files, obsids)
        logger.debug(f"Number of valid files: {len(files)}")

    positions = []
    for file in files:
        with fits.open(file, mode="readonly", memmap=True) as in_fits:
            if ra_field not in in_fits[0].header:
                logger.debug(f"{ra_field} not in {file}")
                continue
            if dec_field not in in_fits[0].header:
                logger.debug(f"{dec_field} not in {file}")
                continue

            try:
                ra = in_fits[0].header[ra_field]
                dec = in_fits[0].header[dec_field]

                positions.append((file, ra, dec))
            except:
                logger.warn(f"Position could not be read for {file}")

    if len(positions) == 0:
        logger.debug("No positions were captured. ")
        return np.nan, np.nan

    ras = [i[1] for i in positions]
    decs = [i[2] for i in positions]

    sky_pos = SkyCoord(ras, decs, unit=(u.deg, u.deg))

    # The use of circmean is accurate only when delta(dec) is small. For the
    # GLEAM-X declination strips this is the case.
    mean_ra = circmean(sky_pos.ra)
    mean_dec = np.mean(sky_pos.dec)
    mean_pos = SkyCoord(mean_ra, mean_dec, unit=(u.deg, u.deg))

    if print_pos:
        print(f"{mean_pos.ra.deg} {mean_pos.dec.deg}")

    return mean_ra, mean_dec


if __name__ == "__main__":

    parser = ArgumentParser(
        description="Given a set of FITS files, extract their RA/Dec center and then calculate the mean position from this set"
    )

    parser.add_argument(
        "files",
        nargs="+",
        type=str,
        help="Path to the (many) FITS files to extract the position from ",
    )

    parser.add_argument(
        "--filter-obsids",
        type=str,
        default=None,
        help="Provided a GLEAM-X obsid text file (new line delimited), filter out files that do not contain a desired obsid. ",
    )

    parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Enable extra output when running",
    )

    parser.add_argument(
        "--ra-field",
        default="CRVAL1",
        type=str,
        help="Field name that contains the RA value",
    )
    parser.add_argument(
        "--dec-field",
        default="CRVAL2",
        type=str,
        help="Field name that contains the Dec value",
    )

    args = parser.parse_args()

    if args.verbose is True:
        logger.setLevel(logging.DEBUG)

    calculate_mean(
        args.files,
        print_pos=True,
        filter_obsids=args.filter_obsids,
        ra_field=args.ra_field,
        dec_field=args.dec_field,
    )
