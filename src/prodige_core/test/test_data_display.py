from __future__ import annotations
import numpy as np
import os

import prodige_core.data_display
from astropy.io import fits
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.testing.decorators import image_comparison

import pytest


def test_pb_telecope_good_frequency() -> None:
    assert (
        prodige_core.data_display.pb_telecope(
            72.78382 * u.GHz, telescope="NOEMA")
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
        prodige_core.data_display.pb_telecope(
            72.78382 * u.GHz, telescope="TEST")


def test_validate_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)


def test_validate_determine_noise_map_bad_input() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.determine_noise_map("test")


def test_get_contour_params() -> None:
    steps_arr, line_style, do_contours = prodige_core.get_contour_params(
        50.0, 1.0)
    assert (steps_arr == [-5.0, 5.0, 10.0, 20.0, 40.0]).all()
    assert line_style == ["dotted"] + ["solid"] * 4
    assert do_contours == True


def test_get_frequency(sample_image) -> None:
    hdu = sample_image(is_2d=True)
    hdr = hdu.header
    assert prodige_core.data_display.get_frequency(hdr) == 72.78382
    with pytest.raises(ValueError):
        prodige_core.data_display.get_frequency(72.78382 * u.GHz)


def test_get_wavelength(sample_image) -> None:
    hdu = sample_image(is_2d=True)
    assert prodige_core.data_display.get_wavelength(hdu.header) == 4.1 * u.mm
    with pytest.raises(ValueError):
        prodige_core.data_display.get_wavelength(72.78382 * u.m)


def test_noise_map() -> None:
    """noise level is well calculated, while also resistent to adding 4 NaNs"""
    rms = 0.1
    data_2d = np.random.normal(0, rms, (500, 500))
    assert (
        pytest.approx(
            prodige_core.data_display.determine_noise_map(data_2d), rel=0.05)
        == rms
    )
    data_2d[0, -1] = np.nan
    data_2d[-1, -1] = np.nan
    data_2d[0, 0] = np.nan
    data_2d[-1, 0] = np.nan
    assert (
        pytest.approx(
            prodige_core.data_display.determine_noise_map(data_2d), rel=0.05)
        == rms
    )


def test_filename_continuum() -> None:
    assert (
        prodige_core.data_display.filename_continuum("test", "li", mosaic=True)
        == "test_CD_li_cont_rob1-selfcal.fits"
    )
    assert (
        prodige_core.data_display.filename_continuum(
            "test", "li", mosaic=False)
        == "test_CD_li_cont_rob1-selfcal-pbcor.fits"
    )


def test_load_continuum_data(tmp_path, sample_image) -> None:
    dir = tmp_path / "sub"
    dir.mkdir()
    file_link = os.path.join(os.fspath(dir), "test_image.fits")
    file_link2 = os.path.join(os.fspath(dir), "test_image2.fits")
    hdu = sample_image(is_2d=True)
    rms = 0.1
    data = np.random.normal(0, rms, hdu.data.shape)  # (501, 501))
    hdu.data = data
    hdu.header["BUNIT"] = "mJy/beam"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_continuum_data(
        file_link, "B1-bS")
    assert (hd["NAXIS1"] == 200) and (hd["NAXIS2"] == 200)
    assert (hd['BUNIT'].casefold() == 'mJy/beam'.casefold())
    assert rms == pytest.approx(rms_out, rel=0.05)

    hdu.header["BUNIT"] = "Jy/beam"
    hdu.writeto(file_link2, overwrite=True)
    _, rms_out2, _ = prodige_core.data_display.load_continuum_data(
        file_link2, "B1-bS")
    assert (hd['BUNIT'].casefold() == 'mJy/beam'.casefold())
    assert rms * 1e3 == pytest.approx(rms_out2, rel=0.05)


@image_comparison(baseline_images=['example_map'], remove_text=True,
                  extensions=['png'], style='mpl20', tol=10)
def test_plot_continuum(tmp_path, sample_image) -> None:
    dir = tmp_path
    dir.mkdir(exist_ok=True)
    file_name = prodige_core.data_display.filename_continuum(
        'B1-bS', 'li', False)
    file_link = os.path.join(os.fspath(dir), file_name)  # "test_image.fits")
    hdu = sample_image(is_2d=True)
    hdu.header['RESTFREQ'] = 216.7230e9
    rms = 0.1
    seed = 122807528840384100672342137672332424406
    rng = np.random.default_rng(seed)
    data = rng.standard_normal(hdu.data.shape) * rms
    hdu.data = data
    hdu.writeto(file_link, overwrite=True)

    prodige_core.data_display.plot_continuum(
        'B1-bS', 'li', os.fspath(dir)+'/', mosaic=False, vmin=-0.5, vmax=2.0, save_fig=False, do_marker=True, do_annotation=True, do_outflow=True)
    
