from __future__ import annotations

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS


try:
    from importlib.resources import files
except ImportError:
    from importlib_resources import files

from .config import source_filename


def load_sources_table() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """
    Load the source table containing the sources within the field of view.
    It uses the source_filename variable from the config.py file.
    """
    # load table containing sources within FoV
    # sources within FOV
    data_file = files("prodige_core").joinpath(source_filename)
    sources_tab = np.loadtxt(data_file, dtype="U",
                             delimiter=",", comments="#", skiprows=1)
    # source name
    name_list = sources_tab[:, 0]
    # source RA
    RA_list = sources_tab[:, 1]
    # source DEC
    Dec_list = sources_tab[:, 2]
    # source color
    color_list = sources_tab[:, 3]
    # source vlsr
    vlsr_list = sources_tab[:, 4].astype(float)
    # source outflow PA
    outflowPA = sources_tab[:, 5].astype(float)
    # label offset PA
    label_offsetPA = sources_tab[:, 6].astype(float)

    return (
        name_list,
        RA_list,
        Dec_list,
        color_list,
        vlsr_list,
        outflowPA,
        label_offsetPA,
    )


# regions dictionary storing the region name, the central RA and DEC (FK5), and
# the size of the region (height and width wiht units)
# these are supposed to be used for plotting purposes
fig_width_def = 6.0
fig_height_def = 6.0
width_def = 40.0 * u.arcsec
height_def = 40.0 * u.arcsec

region_dic = {
    # these are the name of the mosaicked regions
    "L1448N": {
        "RA0": "3:25:36.44",
        "Dec0": "30:45:18.3",
        "height": 33 * u.arcsec,
        "width": 30 * u.arcsec,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    # The following sources are the original pointings
    "B5-IRS1": {
        "RA0": "03:47:41.591",
        "Dec0": "32:51:43.672",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "L1448IRS3A": {
        "RA0": "03:25:36.499",
        "Dec0": "30:45:21.880",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "L1448-IRS3B": {
        "RA0": "03:25:36.379",
        "Dec0": "30:45:14.728",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "L1448C": {
        "RA0": "03:25:38.875",
        "Dec0": "30:44:05.283",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-2": {
        "RA0": "03:32:17.928",
        "Dec0": "30:49:47.825",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-5": {
        "RA0": "03:31:20.939",
        "Dec0": "30:45:30.273",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-8": {
        "RA0": "03:44:43.982",
        "Dec0": "32:01:35.210",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-18": {
        "RA0": "03:29:11.258",
        "Dec0": "31:18:31.073",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-22": {
        "RA0": "03:25:22.409",
        "Dec0": "30:45:13.258",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-29": {
        "RA0": "03:33:17.877",
        "Dec0": "31:09:31.817",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-30": {
        "RA0": "03:33:27.303",
        "Dec0": "31:07:10.160",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-50": {
        "RA0": "03:29:07.768",
        "Dec0": "31:21:57.128",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "Per-emb-62": {
        "RA0": "03:44:12.977",
        "Dec0": "32:01:35.419",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "B1-bN": {
        "RA0": "03:33:21.209",
        "Dec0": "31:07:43.665",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "B1-bS": {
        "RA0": "03:33:21.355",
        "Dec0": "31:07:26.372",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "HH211": {
        "RA0": "03:43:56.806",
        "Dec0": "32:00:50.202",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "IC348MMS": {
        "RA0": "03:43:57.065",
        "Dec0": "32:03:04.788",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "IRAS2A": {
        "RA0": "03:28:55.569",
        "Dec0": "31:14:37.025",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "IRAS2B": {
        "RA0": "03:28:57.374",
        "Dec0": "31:14:15.772",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "IRAS4B": {
        "RA0": "03:29:12.016",
        "Dec0": "31:13:08.031",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "IRAS4C": {
        "RA0": "03:29:13.548",
        "Dec0": "31:13:58.150",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "L1448NW": {
        "RA0": "03:25:35.671",
        "Dec0": "30:45:34.193",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "SVS13A": {
        "RA0": "03:29:03.764",
        "Dec0": "31:16:03.808",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
    "SVS13B": {
        "RA0": "03:29:03.078",
        "Dec0": "31:15:51.740",
        "height": height_def,
        "width": width_def,
        "fig_width": fig_width_def,
        "fig_height": fig_height_def,
    },
}


def load_cutout(
    file_in: str, source: str = "L1448N", is_hdu: bool = False
) -> fits.hdu.hdulist.HDUList:
    """
    Convenience function to load a FITS file, or an existing HDU,
    and to generate a cutout following the requested center and sizes.

    params:

    file_in : FITS file name to be loaded (if is_hdu=False). If is_hdu=True,
              the file_in represents a valid HDU to be used.
    source : String parameter. It represents the source name to be used to
             define the center of the cutout. This will be loaded from the
             region dictionary.
    is_hdu : Boolean parameter. It controls what is passed to the function.

    Return:
    It returns the HDU of the cutout defined by the global variables: position and cutout_size.
    """
    validate_source_id(source)
    position = SkyCoord(
        region_dic[source]["RA0"] + " " + region_dic[source]["Dec0"],
        unit=(u.hourangle, u.deg),
    )
    cutout_size = u.Quantity(
        (region_dic[source]["height"], region_dic[source]["width"])
    )

    if is_hdu == False:
        hdu = fits.open(file_in)[0]
    else:
        hdu = file_in.copy()
    # Make the cutout, including the WCS
    if hdu.header["NAXIS"] != 2:
        wcs_original = WCS(hdu.header).dropaxis(2)
    else:
        wcs_original = WCS(hdu.header)
    cutout = Cutout2D(
        np.squeeze(hdu.data), position=position, size=cutout_size, wcs=wcs_original
    )
    # update HDU with new data and updated header information
    hdu.data = cutout.data
    hdu.header.update(cutout.wcs.to_header())
    # list of keys to be removed (if they exist)
    key_list = [
        "CTYPE3",
        "CRVAL3",
        "CRPIX3",
        "CDELT3",
        "CUNIT3",
        "NAXIS3",
        "CROTA3",
        "SPECSYS",
        "VELREF",
    ]
    for key in key_list:
        if key in hdu.header:
            del hdu.header[key]
    return hdu


def validate_source_id(source: str) -> bool:
    """
    Convenience function to validate the source name.
    """
    if source not in region_dic:
        raise ValueError("Source not in source_id list")
    return True


def get_figsize(source: str) -> tuple[float, float]:
    """
    Convenience function to get the figure size for a given source.
    """
    validate_source_id(source)
    fig_width = region_dic[source]["fig_width"]
    fig_height = region_dic[source]["fig_height"]
    return fig_width, fig_height


def get_region_center(source: str) -> tuple[u.Quantity, u.Quantity]:
    """
    Convenience function to get the region center for a given source.
    """
    validate_source_id(source)
    position = SkyCoord(
        region_dic[source]["RA0"] + " " + region_dic[source]["Dec0"],
        unit=(u.hourangle, u.deg),
    )
    ra0 = position.ra.deg
    dec0 = position.dec.deg
    return ra0, dec0


def get_region_names() -> list[str]:
    """
    Convenience function to get the available region names.
    """
    # return region_dic.keys()
    return list(region_dic)


def get_region_vlsr(source: str) -> float:
    """
    Convenience function to get the vlsr for a given source.
    """
    source_name, _, _, _, vlsr_source, _, _ = load_sources_table()
    validate_source_id(source)
    idx = np.where(source_name == source)
    return vlsr_source[idx][0]


def get_outflow_information() -> (
    tuple[list[str], list[str], list[str], np.ndarray, np.ndarray]
):
    """
    Convenience function to get the list of region center for a given source.
    """
    sources_outflowPA, RA_source, Dec_source, _, _, outflowPA, label_offsetPA = (
        load_sources_table()
    )
    return sources_outflowPA, RA_source, Dec_source, outflowPA, label_offsetPA
