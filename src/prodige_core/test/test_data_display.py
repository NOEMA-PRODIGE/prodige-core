from __future__ import annotations
import numpy as np

import prodige_core.data_display
from astropy.io import fits
from astropy import units as u
from astropy.io import fits

import pytest


def test_pb_telecope_good_frequency() -> None:
    assert (
        prodige_core.data_display.pb_telecope(72.78382 * u.GHz, telescope="NOEMA")
        == 64.1 * u.arcsec
    )
    assert (
        prodige_core.data_display.pb_telecope(345 * u.GHz, telescope="SMA")
        == 36.0 * u.arcsec
    )
    assert (
        prodige_core.data_display.pb_telecope(1 * u.GHz, telescope="VLA")
        == 45.0 * u.arcmin
    )
    assert (
        prodige_core.data_display.pb_telecope(300 * u.GHz, telescope="ALMA")
        == 19.0 * u.arcsec
    )
    with pytest.raises(ValueError):
        prodige_core.data_display.pb_telecope(72.78382 * u.GHz, telescope="TEST")


def test_validate_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)


def test_validate_determine_noise_map_bad_input() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.determine_noise_map("test")


def test_get_contour_params() -> None:
    steps_arr, line_style = prodige_core.get_contour_params(50.0, 1.0)
    assert (steps_arr == [-5.0, 5.0, 10.0, 20.0, 40.0]).all()
    assert line_style == ["dotted"] + ["solid"] * 4


def test_get_frequency() -> None:
    hdr = fits.Header()
    hdr["RESTFREQ"] = (72.78382e9, "Hz")
    assert prodige_core.data_display.get_frequency(hdr) == 72.78382
    with pytest.raises(ValueError):
        prodige_core.data_display.get_frequency(72.78382 * u.GHz)


def test_get_wavelength() -> None:
    hdr = fits.Header()
    hdr["RESTFREQ"] = (250.0e9, "Hz")
    assert prodige_core.data_display.get_wavelength(hdr) == 1.2 * u.mm
    with pytest.raises(ValueError):
        prodige_core.data_display.get_wavelength(72.78382 * u.m)


def test_noise_map() -> None:
    """noise level is well calculated, while also resistent to adding 4 NaNs"""
    rms = 0.1
    data_2d = np.random.normal(0, rms, (500, 500))
    assert (
        pytest.approx(prodige_core.data_display.determine_noise_map(data_2d), rel=0.05)
        == rms
    )
    data_2d[0, -1] = np.nan
    data_2d[-1, -1] = np.nan
    data_2d[0, 0] = np.nan
    data_2d[-1, 0] = np.nan
    assert (
        pytest.approx(prodige_core.data_display.determine_noise_map(data_2d), rel=0.05)
        == rms
    )


def test_filename_continuum() -> None:
    assert (
        prodige_core.data_display.filename_continuum("test", "li", mosaic=True)
        == "test_CD_li_cont_rob1-selfcal.fits"
    )
    assert (
        prodige_core.data_display.filename_continuum("test", "li", mosaic=False)
        == "test_CD_li_cont_rob1-selfcal-pbcor.fits"
    )

def test_load_continuum_data(tmp_path) -> None:
    dir = tmp_path / "sub"
    dir.mkdir()
    file_link = dir / "test_image.fits"
    # data = np.ones((101, 101))
    rms = 0.1
    data = np.random.normal(0, rms, (501, 501))
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
    hdu.header["BUNIT"] = "mJy/beam"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_continuum_data(
        file_link, "B1-bS"
    )
    assert (hd["NAXIS1"] == 200) and (hd["NAXIS2"] == 200)
    assert (rms == pytest.approx(rms_out, rel=0.05))
    
    hdu.header["BUNIT"] = "Jy/beam"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out2, _ = prodige_core.data_display.load_continuum_data(
        file_link, "B1-bS"
    )
    assert (rms*1e3 == pytest.approx(rms_out2, rel=0.05))