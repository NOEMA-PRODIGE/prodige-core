from __future__ import annotations
import numpy as np

import prodige_core.data_display
from astropy import units as u
from astropy.io import fits

import pytest


def test_pb_telecope_good_frequency() -> None:
    assert prodige_core.data_display.pb_telecope(
        72.78382*u.GHz, telescope='NOEMA') == 64.1 * u.arcsec
    assert prodige_core.data_display.pb_telecope(
        345*u.GHz, telescope='SMA') == 36.0 * u.arcsec
    assert prodige_core.data_display.pb_telecope(
        1*u.GHz, telescope='VLA') == 45.0 * u.arcmin
    assert prodige_core.data_display.pb_telecope(
        300*u.GHz, telescope='ALMA') == 19.0*u.arcsec


def test_pb_telecope_bad_telescope() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.pb_telecope(
            72.78382 * u.GHz, telescope='TEST')


def test_validate_frequency_bad_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)


def test_validate_frequency_good_frequency() -> None:
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)


def test_validate_determine_noise_map() -> None:
    rms = 0.5
    # make sure the noise map is valid
    # create a fake noise map with the same rms value and size 500x500
    noise_map = np.random.normal(0, rms, (500, 500))
    assert pytest.approx(
        prodige_core.data_display.determine_noise_map(noise_map), rel=1e-1) == rms


def test_validate_determine_noise_map_bad_input() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.determine_noise_map('test')


def test_get_contour_params() -> None:
    steps_arr, line_style = prodige_core.get_contour_params(50.0, 1.0)
    assert (steps_arr == [-5., 5., 10., 20., 40.]
            ).all() and (line_style == ['dotted'] + ['solid']*4)


def test_filename_continuum() -> None:
    assert prodige_core.data_display.filename_continuum(
        'test', 'li', mosaic=True) == "test_CD_li_cont_rob1-selfcal.fits"
    assert prodige_core.data_display.filename_continuum(
        'test', 'li', mosaic=False) == "test_CD_li_cont_rob1-selfcal-pbcor.fits"


def test_get_frequency() -> None:
    hdu = fits.PrimaryHDU()
    hdu.header["RESTFREQ"] = (72.78382e9, 'Hz')
    assert prodige_core.data_display.get_frequency(hdu.header) == 72.78382


def test_get_wavelength() -> None:
    hdu = fits.PrimaryHDU()
    hdu.header["RESTFREQ"] = (250.0e9, 'Hz')
    assert prodige_core.data_display.get_wavelength(hdu.header) == 1.2*u.mm
