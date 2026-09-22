from __future__ import annotations

from collections.abc import Callable
from typing import Protocol

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits

import prodige_core.source_catalogue


class SampleImageFactory(Protocol):
    def __call__(self, is_2d: bool = True) -> fits.PrimaryHDU: ...


@pytest.fixture
def sample_image() -> SampleImageFactory:
    def make_sample_image(is_2d: bool = True) -> fits.PrimaryHDU:
        if is_2d:
            data = np.ones((501, 501))
        else:
            data = np.ones((1, 501, 501))
        ra0, dec0 = prodige_core.source_catalogue.get_region_center("B1-bS")
        header = fits.Header()
        header["CRVAL1"] = ra0.value
        header["CRVAL2"] = dec0.value
        header["CRPIX1"] = 251
        header["CRPIX2"] = 251
        header["CDELT1"] = 40.0 * u.arcsec.to(u.deg) / 200
        header["CDELT2"] = 40.0 * u.arcsec.to(u.deg) / 200
        header["CUNIT1"] = "deg"
        header["CUNIT2"] = "deg"
        header["CTYPE1"] = "RA---TAN"
        header["CTYPE2"] = "DEC--TAN"
        header["EQUINOX"] = 2000.0
        header["RADESYS"] = ("FK5", "Coordinate system")
        header["RESTFREQ"] = (72.78382e9, "Hz")
        header["BUNIT"] = ("mJy/Beam", "Brightness unit")
        header["BMAJ"] = 0.26e-3
        header["BMIN"] = 0.16e-3
        header["BPA"] = 22.0
        hdu = fits.PrimaryHDU(data=data, header=header)
        return hdu

    return make_sample_image
