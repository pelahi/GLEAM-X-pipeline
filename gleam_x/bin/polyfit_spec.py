#! /usr/bin/env python

import os
import logging
from typing import Iterable, Union, Tuple
from glob import glob

import numpy as np
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, vstack, hstack
from argparse import ArgumentParser

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())


EXTS = ["", "_bkg", "_rms", "_weight"]


class FitsObseration:
    """Helper class to manage the fits images and their associated
    products for this rescaling
    """

    def __init__(self, fitsimage: str, ref_cata: Union[str, Table]):
        logger.info(f"Loading information for {fitsimage}")
        self.ref_cata = ref_cata

        self.fitsimage = fitsimage
        self.obsid = self.fitsimage[:10]
        self.components = self.fitsimage.replace(".fits", "_comp.fits")
        self.path, self.fitsname = os.path.split(self.fitsimage)

        self.metafits = glob(f"{self.path}/{self.obsid}*metafits*")[0]
        self.metaheader = fits.getheader(self.metafits)

        self.header = fits.getheader(self.fitsimage)
        self.time = Time(int(self.obsid), format="gps")
        self.centfreq = self.header["CRVAL3"] * 1e6  # Hertz to MegaHerta

        self.table = Table.read(self.components)
        self.table["RACEN"] = self.metaheader["RA"]
        self.table["DECCEN"] = self.metaheader["DEC"]
        self.table["SNR"] = self.table["int_flux"] / self.table["local_rms"]

        self.table["hourangle"] = self.table["RAJ2000"] - self.table["RACEN"]

        self._cross_match_to_ref()
        self._evaluate_sed_model()

    def _cross_match_to_ref(self):
        logger.debug("Performing cross-match to reference catalogue. ")
        ref_tab = (
            Table.read(self.ref_cata)
            if isinstance(self.ref_cata, str)
            else self.ref_cata
        )

        # Load in coordinates
        self.ref_sky = SkyCoord(ref_tab.RAJ2000 * u.deg, ref_tab.DEJ2000 * u.deg)
        self.comp_sky = SkyCoord(self.table.RAJ2000 * u.deg, self.table.DEJ2000 * u.deg)

        m_idx, sep, _ = match_coordinates_sky(
            self.comp_sky, self.ref_sky, nthneighbor=1
        )
        sep_mask = sep < (45 * u.arcsecond)
        logger.debug(f"Found {np.sum(sep_mask)} matches... ")

        self.match_table = hstack(self.table[sep_mask], ref_tab[m_idx][sep_mask])

    def _evaluate_sed_model(self):
        logger.debug("Evaluating SED and calculating ratio")
        self.table[f"model_flux_200"] = (
            self.table["S_200"] * (self.centfreq / 200) ** self.table["alpha"]
        )
        self.table["logratio"] = np.log10(
            self.table[f"mode_flux_200"] / self.table["int_flux"]
        )

    def correct_table(self, dec_model, ra_model):
        logger.debug(f"Correcting table with derived models")
        outpath = self.components.replace(".fits", "_rescaled.fits")

        cols = [
            "background",
            "local_rms",
            "peak_flux",
            "err_peak_flux",
            "int_flux",
            "err_int_flux",
            "residual_mean",
            "residual_std",
        ]
        dec_corr = 10 ** dec_model(self.table["DEJ2000"])
        ra_corr = (
            10 ** ra_model(self.table["RAJ2000"]) if ra_model is not None else None
        )
        for col in cols:
            self.table[col] *= dec_corr
            if ra_corr is not None:
                self.table[col]

        logger.info(f"Saving adjusted table to {outpath}")
        self.table.write(outpath)


def fit_polynomial_model(
    table: Table, x_label: str, y_label: str, w_label: str, order: int = 5
) -> np.polynomial.Polynomial:
    logger.debug(f"Fitting {order} degree polynomial to ({x_label}, {y_label})")
    poly_result = np.polynomial.Polynomial.fit(
        table[x_label], table[y_label], w=table[w_label], deg=order
    )
    logger.debug(f"Polynomial for {y_label} fit finished")

    return poly_result


def sigma_clip_table(table: Table, label: str = "logratio", sigma: float = 1) -> Table:
    logger.debug(f"Running the sigma clip on {label} with {sigma}")
    label_median = np.nanmedian(table[label])
    label_std = np.nanstd(table[label])

    mask = (table[label] > label_median - sigma * label_std) & (
        table[label] < label_median + sigma * label_std
    )

    logger.debug(f"{np.sum(mask)} of {len(table)} were retained after sigma clipping")
    return table[mask]


def filter_sources(table: Table, snr_thres: float = None, no_srcs: int = None) -> Table:

    if snr_thres is not None:
        logger.debug(f"Apply a minimum SNR threshold cut of {snr_thres}")
        snr_mask = table["SNR"] > snr_thres
        table = table[snr_mask]
        logger.debug(f"{np.sum(snr_mask)} sources were above {snr_thres}")

    elif no_srcs is not None:
        logger.debug(f"Saving only the brightest {min(no_srcs, len(table))}")
        order = np.argsort(table["SNR"])[-no_srcs:]
        table = table[order]

    return table


def apply_dec_model(
    table: Table, model: np.polynomial.Polynomial, freq: float
) -> Table:
    table["int_flux_dec_norm"] = table["int_flux"] * 10 ** model(table["DEJ2000"])
    table["logratio_dec_norm"] = np.log10(
        table["S_200"] * (freq / 200) ** table["alpha"] / table["logratio_dec_norm"]
    )


def apply_ra_model(table: Table, model: np.polynomial.Polynomial, freq: float) -> Table:
    table["int_flux_ra_norm"] = table["int_flux"] * 10 ** model(table["DEJ2000"])
    table["logratio_ra_norm"] = np.log10(
        table["S_200"] * (freq / 200) ** table["alpha"] / table["logratio_ra_norm"]
    )


def save_coefficents(models: Iterable) -> None:
    dec_model, ra_model = models
    logger.info("save_coefficents stub file called. Need to expand. ")

    pass


def save_table(table: Table, outpath: str) -> None:
    logger.info(f"Saving night table to {outpath}")
    table.write(outpath)


def derive_correction_models(
    table: Table, obs_freq: float, apply_ra: bool = False, save_coeffs: bool = True
):
    dec_model = fit_polynomial_model(table, "DEJ2000", "logratio", "SNR")
    table = apply_dec_model(table, dec_model, obs_freq)

    ra_model = None
    if apply_ra:
        ra_model = fit_polynomial_model(table, "hourangle", "logratio_dec_norm", "SNR")

    if save_coeffs:
        save_coefficents([dec_model, ra_model])

    return dec_model, ra_model


def correct_night_comps(
    obs: Iterable[Table],
    dec_model: np.polynomial.Polynomial,
    ra_model: np.polynomial.Polynomial = None,
):
    logger.infor("Saving corrected component tables")
    for o in obs:
        o.correct_table(dec_model, ra_model)


def create_correction_screens(
    obs: Iterable[FitsObseration],
    dec_model: np.polynomial.Polynomial,
    ra_model: np.polynomial.Polynomial,
):
    screens = {}
    for o in obs:
        for e in EXTS:
            f = o.fitsname.replace(".fits", f"{e}.fits")
            with fits.open(f, memmap=True) as curr_fits:
                data_shape = curr_fits[0].data.shape
                w = WCS(curr_fits[0].header).celestial

                if data_shape in screens.keys():
                    continue

                idxs = np.indices(data_shape)
                ra, dec = w.pixel_to_world(idxs)

                corr = 10 ** dec_model(dec)

                if ra_model is not None:
                    ra = ra - o["RACEN"]
                    corr *= 10 ** ra_model(ra)

                screens[data_shape] = corr

    return screens


def apply_correction_screens(obs: Iterable[FitsObseration]):

    pass


def derive_apply_spatial_corrections(
    filelist: str,
    ref_cata: str,
    snr_thres: float = None,
    no_srcs: int = None,
    save_night: str = None,
    apply_ra: bool = False,
    correct_all: bool = True,
):

    images = np.loadtxt(filelist)
    logger.info(f"{len(images)} found in {filelist}")

    obs = [FitsObseration(f, ref_cata) for f in images]
    night_table = vstack([o.table for o in obs])

    night_table = sigma_clip_table(night_table, 3)
    night_table = filter_sources(night_table, snr_thres=snr_thres, no_srcs=no_srcs)

    if isinstance(save_night, str):
        save_table(night_table, save_night)

    dec_model, ra_model = derive_correction_models(
        night_table, obs[0].centfreq, apply_ra=apply_ra, save_coeffs=True
    )
    save_table(night_table, "processed_matched_catalogue.fits")

    if correct_all:
        correct_night_comps(obs, dec_model, ra_model)

        screens = create_correction_screens(obs, dec_model=dec_model, ra_model=ra_model)

        # TODO: Apply corrections to all images


if __name__ == "__main__":
    try:
        refcat = f"{os.environ['GXBASE']}//models/GGSM_sparse_unresolved.fits"
    except:
        refcat = None

    parser = ArgumentParser(
        description="Accumulate source statistics across a night to derive a correction screen as a function of RA and Dec"
    )

    parser.add_argument(
        "filelist",
        type=str,
        help="Path to the GLEAM-X images (nominally those produced from fits_warp: *pb_warp.fits) from which to base processing on. Typical naming conventions are assumed throughout when searching from components and imaged.",
    )
    parser.add_argument(
        "--ra-field",
        type=str,
        default="RAJ2000",
        help="The field containing RA position information",
    )
    parser.add_argument(
        "--dec-field",
        type=str,
        default="DEJ2000",
        help="The field containing Dec position information",
    )
    parser.add_argument(
        "-s",
        "--snr",
        default=10,
        type=int,
        help="The minimum SNR for sources for them to be included in the ensemble to fit",
    )
    parser.add_argument(
        "-n",
        "--nsrcs",
        default=10000,
        type=int,
        help="The maximum number of sources to include in the fitting",
    )

    parser.add_argument(
        "--reference-cata",
        type=str,
        default=refcat,
        help="Path to the reference catalogue to use",
    )
    parser.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Enable debug logging when running",
    )

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

