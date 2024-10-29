from __future__ import annotations
import numpy as np

import prodige_core.data_display
from astropy.io import fits
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
# def test_pb_telecope_bad_telescope() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.pb_telecope(
            72.78382 * u.GHz, telescope='TEST')

def test_validate_frequency() -> None:
# def test_validate_frequency_bad_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)
# def test_validate_frequency_good_frequency() -> None:
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)





def test_validate_determine_noise_map_bad_input() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.determine_noise_map('test')


def test_get_contour_params() -> None:
    steps_arr, line_style = prodige_core.get_contour_params(50.0, 1.0)
    assert (steps_arr == [-5., 5., 10., 20., 40.]
            ).all() and (line_style == ['dotted'] + ['solid']*4)


def test_get_frequency() -> None:
    hdr = fits.Header()
    hdr['RESTFREQ'] = (72.78382e9, 'Hz')
    assert prodige_core.data_display.get_frequency(hdr) == 72.78382
# def test_get_frequency_bad_frequency() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.get_frequency(72.78382 * u.GHz)

    
def test_get_wavelength() -> None:
    hdr = fits.Header()
    hdr['RESTFREQ'] =  = (250.0e9, 'Hz')
    assert prodige_core.data_display.get_wavelength(hdr) == 1.2*u.mm
# 
# def test_get_wavelength_bad_wavelength() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.get_wavelength(72.78382 * u.m)


def test_noise_map() -> None:
    """ noise level is well calculated, while also resistent to adding 4 NaNs"""
    rms = 0.1
    # data_2d = np.random.randn(200, 200) * rms
    data_2d = np.random.normal(0, rms, (500, 500))
    assert pytest.approx(prodige_core.data_display.determine_noise_map(data_2d), rel=0.05) == rms
    data_2d[0, -1] = np.nan
    data_2d[-1, -1] = np.nan
    data_2d[0, 0] = np.nan
    data_2d[-1, 0] = np.nan
    assert pytest.approx(prodige_core.data_display.determine_noise_map(data_2d), rel=0.05) == rms
# def test_validate_determine_noise_map() -> None:
#     rms = 0.5
    # make sure the noise map is valid
    # create a fake noise map with the same rms value and size 500x500
    # noise_map = np.random.normal(0, rms, (500, 500))
    # assert pytest.approx(
    #     prodige_core.data_display.determine_noise_map(noise_map), rel=1e-1) == rms

def test_get_contour_params() -> None:
    lev_param = [-5.0, 5.0, 10.0]
    line_par = ['dotted'] + ['solid']*2
    ret_lev_param, ret_line_par = prodige_core.data_display.get_contour_params(5.1, 1.0)
    assert  (ret_lev_param == lev_param).all()
    assert ret_line_par == line_par



def test_filename_continuum() -> None:
    assert prodige_core.data_display.filename_continuum(
        'test', 'li', mosaic=True) == "test_CD_li_cont_rob1-selfcal.fits"
    assert prodige_core.data_display.filename_continuum(
        'test', 'li', mosaic=False) == "test_CD_li_cont_rob1-selfcal-pbcor.fits"

