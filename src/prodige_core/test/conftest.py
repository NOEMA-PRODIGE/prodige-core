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


class SampleImageFactoryVel(Protocol):
    def __call__(self, is_vlsr: bool = True) -> fits.PrimaryHDU: ...


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
        header["CDELT1"] = -40.0 * u.arcsec.to(u.deg) / 200  # type: ignore
        header["CDELT2"] = 40.0 * u.arcsec.to(u.deg) / 200  # type: ignore
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


@pytest.fixture
def sample_image_line() -> SampleImageFactoryVel:
    def make_sample_image_line(is_vlsr: bool = True) -> fits.PrimaryHDU:
        data = np.ones((501, 501))
        if is_vlsr:
            data = np.ones((501, 501))
            bunit_str = ("km/s", "Centroid Velocity")
        else:
            data = np.ones((1, 501, 501))
            bunit_str = ("mJy/Beam km/s", "Integrated intensity unit")
        ra0, dec0 = prodige_core.source_catalogue.get_region_center("B1-bS")
        header = fits.Header()
        header["CRVAL1"] = ra0.value  # type: ignore
        header["CRVAL2"] = dec0.value  # type: ignore
        header["CRPIX1"] = 251
        header["CRPIX2"] = 251
        header["CDELT1"] = -40.0 * u.arcsec.to(u.deg) / 200  # type: ignore
        header["CDELT2"] = 40.0 * u.arcsec.to(u.deg) / 200  # type: ignore
        header["CUNIT1"] = "deg"
        header["CUNIT2"] = "deg"
        header["CTYPE1"] = "RA---TAN"
        header["CTYPE2"] = "DEC--TAN"
        header["EQUINOX"] = 2000.0
        header["RADESYS"] = ("FK5", "Coordinate system")
        header["RESTFREQ"] = (72.78382e9, "Hz")
        header["BUNIT"] = bunit_str
        header["BMAJ"] = 0.26e-3
        header["BMIN"] = 0.16e-3
        header["BPA"] = 22.0
        hdu = fits.PrimaryHDU(data=data, header=header)
        return hdu

    return make_sample_image_line
