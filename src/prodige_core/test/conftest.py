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
            data = np.ones((501, 501))
        else:
            data = np.ones((1, 501, 501))
        ra0, dec0 = prodige_core.source_catalogue.get_region_center("B1-bS")
        hdu = fits.PrimaryHDU(data=data)
        hdu.header["CRVAL1"] = ra0
        hdu.header["CRVAL2"] = dec0
        hdu.header["CRPIX1"] = 251
        hdu.header["CRPIX2"] = 251
        hdu.header["CDELT1"] = 40.0 * u.arcsec.to(u.deg) / 200
        hdu.header["CDELT2"] = 40.0 * u.arcsec.to(u.deg) / 200
        hdu.header["CUNIT1"] = "deg"
        hdu.header["CUNIT2"] = "deg"
        hdu.header["CTYPE1"] = "RA---TAN"
        hdu.header["CTYPE2"] = "DEC--TAN"
        hdu.header["EQUINOX"] = 2000.0
        hdu.header["RADESYS"] = ("FK5", 'Coordinate system')
        hdu.header["RESTFREQ"] = (72.78382e9, "Hz")
        hdu.header["BUNIT"] = ("mJy/Beam", "Brightness unit")
        hdu.header["BMAJ"] = 0.26E-3
        hdu.header["BMIN"] = 0.16E-3
        hdu.header["BPA"] = 22.0
        return hdu
    return make_sample_image
