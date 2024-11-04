from __future__ import annotations
import pytest
import numpy as np
from astropy.io import fits
import prodige_core.source_catalogue
from astropy import units as u


@pytest.fixture
def sample_image() -> fits.PrimaryHDU:
    def make_sample_image(is_2d: bool = True) -> fits.PrimaryHDU:
        if is_2d:
            data = np.ones((101, 101))
        else:
            data = np.ones((1, 101, 101))
        ra0, dec0 = prodige_core.source_catalogue.get_region_center("B1-bS")
        hdu = fits.PrimaryHDU(data=data)
        hdu.header["CRVAL1"] = ra0
        hdu.header["CRVAL2"] = dec0
        hdu.header["CRPIX1"] = 51
        hdu.header["CRPIX2"] = 51
        hdu.header["CDELT1"] = 40.0 * u.arcsec.to(u.deg) / 20
        hdu.header["CDELT2"] = 40.0 * u.arcsec.to(u.deg) / 20
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["RESTFREQ"] = (72.78382e9, "Hz")
        hdu.header["BUNIT"] = ("Jy/Beam", "Brightness unit")
        return hdu
    return make_sample_image
