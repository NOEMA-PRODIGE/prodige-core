from __future__ import annotations
import numpy as np

import prodige_core.data_display
from astropy.io import fits
from astropy import units as u

import pytest


def test_pb_telecope_good_frequency() -> None:
    assert prodige_core.data_display.pb_telecope(72.78382*u.GHz, telescope='NOEMA') == 64.1 * u.arcsec
    assert prodige_core.data_display.pb_telecope(345*u.GHz, telescope='SMA') == 36.0 * u.arcsec
    assert prodige_core.data_display.pb_telecope(1*u.GHz, telescope='VLA') == 45.0 * u.arcmin
    assert prodige_core.data_display.pb_telecope(300*u.GHz, telescope='ALMA') == 19.0*u.arcsec
    
def test_pb_telecope_bad_telescope() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.pb_telecope(72.78382 * u.GHz, telescope='TEST')

def test_validate_frequency_bad_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)

def test_validate_frequency_good_frequency() -> None:
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)

def test_get_wavelength_bad_wavelength() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.get_wavelength(72.78382 * u.m)

def test_get_wavelength() -> None:
    hdr = fits.Header()
    wave0 = 1e-3*u.m
    print(wave0)
    hdr['RESTFREQ'] = wave0.to(u.Hz, equivalencies=u.spectral()).value
    print(hdr)
    assert np.around(wave0.to(u.mm), decimals=1) == prodige_core.data_display.get_wavelength(hdr)

def test_get_frequency_bad_frequency() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.get_frequency(72.78382 * u.GHz)

def test_get_frequency() -> None:
    hdr = fits.Header()
    freq0 = 72.78382 * u.GHz
    hdr['RESTFREQ'] = freq0.to(u.Hz).value
    assert freq0.to(u.GHz).value == prodige_core.data_display.get_frequency(hdr)

def test_noise_map() -> None:
    """ noise level is well calculated, while also resistent to adding 4 NaNs"""
    rms = 0.1
    data_2d = np.random.randn(200, 200) * rms
    assert pytest.approx(prodige_core.data_display.determine_noise_map(data_2d), rel=0.05) == rms
    # data_2d[0, -1] = np.nan
    # data_2d[-1, -1] = np.nan
    data_2d[0, 0] = np.nan
    # data_2d[-1, 0] = np.nan
    assert pytest.approx(prodige_core.data_display.determine_noise_map(data_2d), rel=0.05) == rms

def test_get_contour_params() -> None:
    lev_param = [-5.0, 5.0, 10.0]
    line_par = ['dotted'] + ['solid']*2
    ret_lev_param, ret_line_par = prodige_core.data_display.get_contour_params(5.1, 1.0)
    assert  (ret_lev_param == lev_param).all()
    assert ret_line_par == line_par

def test_get_filename() -> None:
    datafile_single = "B1_CD_li_cont_rob1-selfcal-pbcor.fits"
    datafile_mosaic = "B1_CD_li_cont_rob1-selfcal.fits"
    assert prodige_core.data_display.filename_continuum('B1', 'li', mosaic=False) == datafile_single
    assert prodige_core.data_display.filename_continuum('B1', 'li', mosaic=True) == datafile_mosaic