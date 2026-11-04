from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import matplotlib
import numpy as np
import pytest
from astropy import units as u
from astropy.io.fits import Header
from astropy.utils.exceptions import AstropyUserWarning
from matplotlib import ft2font
from matplotlib.figure import Figure

import prodige_core.data_display

from .conftest import (  # Import the Protocol type
    SampleImageFactory,
    SampleImageFactoryVel,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.ion()

IMAGE_COMPARISON_DPI = 150
IMAGE_COMPARISON_MATPLOTLIB = "3.11.2"
IMAGE_COMPARISON_FREETYPE = "2.14.3"


def test_image_comparison_rendering_versions() -> None:
    """Keep the PNG baselines tied to one deterministic rendering stack."""
    assert matplotlib.__version__ == IMAGE_COMPARISON_MATPLOTLIB
    assert ft2font.__freetype_version__ == IMAGE_COMPARISON_FREETYPE


def test_pb_telecope_good_frequency() -> None:
    assert (
        prodige_core.data_display.pb_telecope(72.78382 * u.GHz, telescope="NOEMA")  # type: ignore
        == 64.1 * u.arcsec
    )
    assert (
        prodige_core.data_display.pb_telecope(345 * u.GHz, telescope="SMA")  # type: ignore
        == 36.0 * u.arcsec  # type: ignore
    )
    assert (
        prodige_core.data_display.pb_telecope(1 * u.GHz, telescope="VLA")  # type: ignore
        == 45.0 * u.arcmin  # type: ignore
    )
    assert (
        prodige_core.data_display.pb_telecope(300 * u.GHz, telescope="ALMA")  # type: ignore
        == 19.0 * u.arcsec  # type: ignore
    )
    assert (
        prodige_core.data_display.pb_telecope(1 * u.GHz, telescope="30m")  # type: ignore
        == 2460.0 * u.arcsec  # type: ignore
    )
    with pytest.raises(ValueError):
        prodige_core.data_display.pb_telecope(
            72.78382 * u.GHz,
            telescope="TEST",  # type: ignore
        )
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.pb_telecope(1 * u.m, telescope="30m")  # type: ignore


def test_validate_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)  # type: ignore
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)  # type: ignore


def test_validate_determine_noise_map_bad_input() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.determine_noise_map("test")  # type: ignore


def test_get_contour_params() -> None:
    steps_arr, line_style, do_contours = prodige_core.get_contour_params(50.0, 1.0)
    assert (steps_arr == [-5.0, 5.0, 10.0, 20.0, 40.0]).all()
    assert line_style == ["dotted"] + ["solid"] * 4
    assert do_contours


@pytest.mark.parametrize(
    ("maximum", "noise"),
    [(np.nan, 1.0), (1.0, 0.0), (1.0, np.nan), (1.0, np.inf)],
)
def test_get_contour_params_invalid_inputs(maximum: float, noise: float) -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.get_contour_params(maximum, noise)


def test_get_contour_params_zero_maximum() -> None:
    levels, styles, do_contours = prodige_core.data_display.get_contour_params(0.0, 1.0)
    assert (levels == [0.0]).all()
    assert styles == ["solid"]
    assert not do_contours


def test_get_frequency(sample_image: SampleImageFactory) -> None:
    hdu = sample_image(is_2d=True)
    hdr: Header = hdu.header
    assert prodige_core.data_display.get_frequency(hdr) == 72.78382
    with pytest.raises((TypeError, ValueError)):
        prodige_core.data_display.get_frequency(72.78382 * u.GHz)  # type: ignore


def test_get_wavelength(sample_image: SampleImageFactory) -> None:
    hdu = sample_image(is_2d=True)
    assert prodige_core.data_display.get_wavelength(hdu.header) == 4.1
    with pytest.raises((TypeError, ValueError)):
        prodige_core.data_display.get_wavelength(72.78382 * u.m)  # type: ignore


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
    with pytest.warns(
        AstropyUserWarning,
        match="Input data contains invalid values (NaNs or infs)*",
    ):
        assert (
            pytest.approx(
                prodige_core.data_display.determine_noise_map(data_2d), rel=0.05
            )
            == rms
        )


@pytest.mark.parametrize("data", [np.zeros((4, 4)), np.full((4, 4), np.nan)])
def test_determine_noise_map_degenerate_data(data: np.ndarray) -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.determine_noise_map(data)


@pytest.mark.parametrize(
    ("panels", "ncols", "labels", "message"),
    [
        ([], 1, None, "panels must contain"),
        ([("B1-bS", "lo")], 0, None, "ncols must be"),
        ([("B1-bS", "lo")], 1, [], "labels must contain"),
    ],
)
def test_plot_continuum_grid_validates_inputs(
    panels: list[tuple[str, str]],
    ncols: int,
    labels: list[str] | None,
    message: str,
) -> None:
    with pytest.raises(ValueError, match=message):
        prodige_core.data_display.plot_continuum_grid(
            panels, data_directory="", ncols=ncols, labels=labels, save_fig=False
        )


def test_line_grid_wrappers_build_expected_panels(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[
        tuple[list[prodige_core.data_display.LinePanelSpec], dict[str, object]]
    ] = []

    def fake_plot_pairs(
        panels: list[prodige_core.data_display.LinePanelSpec],
        data_directory: str,
        **kwargs: object,
    ) -> tuple[None, list[object]]:
        assert data_directory == "data/"
        calls.append((panels, kwargs))
        return None, []

    monkeypatch.setattr(
        prodige_core.data_display, "plot_source_line_pairs", fake_plot_pairs
    )
    prodige_core.data_display.plot_lines_for_source(
        "B1-bS", ["DCO+_3-2", "H2CO_l21"], "data/", map_kind="velocity"
    )
    prodige_core.data_display.plot_sources_for_line(
        ["B1-bS", "L1448N"], "DCO+_3-2", "data/"
    )
    prodige_core.data_display.plot_source_line_grid(
        ["B1-bS", "L1448N"], ["DCO+_3-2", "H2CO_l21"], "data/"
    )

    assert [spec.map_kind for spec in calls[0][0]] == ["velocity", "velocity"]
    assert [spec.source for spec in calls[1][0]] == ["B1-bS", "L1448N"]
    assert [(spec.source, spec.line) for spec in calls[2][0]] == [
        ("B1-bS", "DCO+_3-2"),
        ("B1-bS", "H2CO_l21"),
        ("L1448N", "DCO+_3-2"),
        ("L1448N", "H2CO_l21"),
    ]
    assert calls[2][1]["ncols"] == 2


def test_plot_source_line_pairs_moment0_and_velocity(
    tmp_path: Path, sample_image_line: SampleImageFactoryVel
) -> None:
    moment_file = tmp_path / prodige_core.data_display.filename_line_TdV(
        "B1-bS", "DCO+_3-2"
    )
    velocity_file = tmp_path / prodige_core.data_display.filename_line_vlsr(
        "B1-bS", "DCO+_3-2"
    )
    moment_hdu = sample_image_line(is_vlsr=False)
    moment_hdu.writeto(moment_file)
    velocity_hdu = sample_image_line(is_vlsr=True)
    velocity_hdu.writeto(velocity_file)

    fig, axs = prodige_core.data_display.plot_source_line_pairs(
        [
            prodige_core.data_display.LinePanelSpec("B1-bS", "DCO+_3-2"),
            prodige_core.data_display.LinePanelSpec(
                "B1-bS", "DCO+_3-2", map_kind="velocity"
            ),
        ],
        os.fspath(tmp_path) + "/",
        ncols=2,
        do_annotation=False,
        show_colorbar=False,
        save_fig=False,
    )
    assert len(axs) == 2
    assert fig.axes == axs
    plt.close(fig)


def test_filename_continuum() -> None:
    assert (
        prodige_core.data_display.filename_continuum("test", "li", mosaic=True)
        == "test_CD_li_cont_rob1-selfcal.fits"
    )
    assert (
        prodige_core.data_display.filename_continuum("test", "li", mosaic=False)
        == "test_CD_li_cont_rob1-selfcal-pbcor.fits"
    )


def test_filename_line_TdV() -> None:
    assert (
        prodige_core.data_display.filename_line_TdV("test", "line", mosaic=True)
        == "test_CD_line_TdV.fits"
    )
    assert (
        prodige_core.data_display.filename_line_TdV("test", "line", mosaic=False)
        == "test_CD_line_TdV.fits"
    )


def test_filename_line_vlsr() -> None:
    assert (
        prodige_core.data_display.filename_line_vlsr("test", "line", mosaic=True)
        == "test_CD_line_Vlsr.fits"
    )
    assert (
        prodige_core.data_display.filename_line_vlsr("test", "line", mosaic=False)
        == "test_CD_line_Vlsr.fits"
    )


def test_load_continuum_data(tmp_path: Path, sample_image: SampleImageFactory) -> None:
    dir = tmp_path / "sub"
    dir.mkdir()
    file_link = os.path.join(os.fspath(dir), "test_image.fits")

    hdu = sample_image(is_2d=True)
    rms = 0.1
    data = np.random.normal(0, rms, hdu.data.shape)  # (501, 501))
    hdu.data = data
    hdu.header["BUNIT"] = "mJy/beam"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_continuum_data(file_link, "B1-bS")
    assert (hd["NAXIS1"] == 200) and (hd["NAXIS2"] == 200)
    assert str(hd["BUNIT"]).casefold() == "mJy/beam".casefold()
    assert rms == pytest.approx(rms_out, rel=0.05)

    hdu.header["BUNIT"] = "Jy/beam"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_continuum_data(file_link, "B1-bS")
    assert str(hd["BUNIT"]).casefold() == "mJy/beam".casefold()
    assert rms * 1e3 == pytest.approx(rms_out, rel=0.05)


def test_load_line_TdV(tmp_path: Path, sample_image_line: SampleImageFactory) -> None:
    dir = tmp_path / "sub"
    dir.mkdir()
    file_link = os.path.join(os.fspath(dir), "test_noise.fits")

    hdu = sample_image_line(is_vlsr=False)
    rms = 0.1
    data = np.random.normal(0, rms, hdu.data.shape)  # (501, 501))
    hdu.data = data
    hdu.header["BUNIT"] = "mJy/beam km/s"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_line_TdV(file_link, "B1-bS")
    assert (hd["NAXIS1"] == 200) and (hd["NAXIS2"] == 200)
    assert str(hd["BUNIT"]).casefold() == "mJy/beam km/s".casefold()
    assert rms == pytest.approx(rms_out, rel=0.05)

    hdu.header["BUNIT"] = "Jy/beam km/s"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_line_TdV(file_link, "B1-bS")
    assert str(hd["BUNIT"]).casefold() == "mJy/beam km/s".casefold()
    assert rms * 1e3 == pytest.approx(rms_out, rel=0.05)

    hdu.header["BUNIT"] = "K km/s"
    hdu.writeto(file_link, overwrite=True)
    _, rms_out, hd = prodige_core.data_display.load_line_TdV(file_link, "B1-bS")
    assert str(hd["BUNIT"]).casefold() == "K km/s".casefold()
    assert rms == pytest.approx(rms_out, rel=0.05)


def _plot_continuum_helper(
    tmp_path: Path,
    sample_image: SampleImageFactory,
    vmin: float | None,
    vmax: float | None,
    do_annotation: bool,
    do_offsets: bool,
) -> Figure:
    """Helper function for continuum plotting tests."""
    dir = tmp_path
    dir.mkdir(exist_ok=True)
    file_name = prodige_core.data_display.filename_continuum(
        "B1-bS", "li", mosaic=False
    )
    file_link = os.path.join(os.fspath(dir), file_name)
    hdu = sample_image(is_2d=True)
    hdu.header["RESTFREQ"] = 216.7230e9
    rms = 0.1
    seed = 122807528840384100672342137672332424406
    rng = np.random.default_rng(seed)
    data = rng.standard_normal(hdu.data.shape) * rms
    hdu.data = data
    hdu.writeto(file_link, overwrite=True)

    with patch.dict(prodige_core.data_display.pyplot_params, {"text.usetex": False}):
        prodige_core.data_display.plot_continuum(
            "B1-bS",
            "li",
            os.fspath(dir) + "/",
            mosaic=False,
            vmin=vmin,
            vmax=vmax,
            save_fig=False,
            do_marker=True,
            do_annotation=do_annotation,
            do_offsets=do_offsets,
            do_outflow=True,
        )
    return plt.gcf()


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_data_display",
    style={"text.usetex": False},
    filename="example_map.png",
    tolerance=10,
    remove_text=True,
    savefig_kwargs={"dpi": IMAGE_COMPARISON_DPI},
)
def test_plot_continuum(tmp_path: Path, sample_image: SampleImageFactory) -> Figure:
    return _plot_continuum_helper(tmp_path, sample_image, -0.5, 2.0, True, False)


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_data_display",
    style={"text.usetex": False},
    filename="example_map_offset.png",
    tolerance=10,
    remove_text=True,
    savefig_kwargs={"dpi": IMAGE_COMPARISON_DPI},
)
def test_plot_continuum_offset(
    tmp_path: Path, sample_image: SampleImageFactory
) -> Figure:
    return _plot_continuum_helper(tmp_path, sample_image, None, None, False, True)


def _plot_line_mom0_helper(
    tmp_path: Path,
    sample_image_line: SampleImageFactoryVel,
    vmin: float | None,
    vmax: float | None,
    do_annotation: bool,
    do_offsets: bool,
) -> Figure:
    """Helper function for line moment 0 plotting tests."""
    dir = tmp_path
    dir.mkdir(exist_ok=True)
    file_name = prodige_core.data_display.filename_line_TdV(
        "B1-bS", "DCO+_3-2", mosaic=False
    )
    file_link = os.path.join(os.fspath(dir), file_name)
    hdu = sample_image_line(is_vlsr=False)
    hdu.header["RESTFREQ"] = 230.538e9
    rms = 0.1
    seed = 122807528840384100672342137672332424406
    rng = np.random.default_rng(seed)
    data = hdu.data + rng.standard_normal(hdu.data.shape) * rms
    hdu.data = data
    hdu.header["BUNIT"] = "mJy/beam km/s"
    hdu.writeto(file_link, overwrite=True)

    with patch.dict(prodige_core.data_display.pyplot_params, {"text.usetex": False}):
        prodige_core.data_display.plot_line_mom0(
            "B1-bS",
            "DCO+_3-2",
            "li",
            os.fspath(dir) + "/",
            mosaic=False,
            vmin=vmin,
            vmax=vmax,
            save_fig=False,
            do_marker=True,
            do_annotation=do_annotation,
            do_offsets=do_offsets,
            do_outflow=True,
        )
    return plt.gcf()


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_data_display",
    filename="example_line_mom0.png",
    style={"text.usetex": False},
    tolerance=10,
    remove_text=True,
    savefig_kwargs={"dpi": IMAGE_COMPARISON_DPI},
)
def test_plot_line_mom0(
    tmp_path: Path, sample_image_line: SampleImageFactoryVel
) -> Figure:
    return _plot_line_mom0_helper(tmp_path, sample_image_line, -0.3, 1.3, True, False)


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_data_display",
    filename="example_line_mom0_offset.png",
    style={"text.usetex": False},
    tolerance=10,
    remove_text=True,
    savefig_kwargs={"dpi": IMAGE_COMPARISON_DPI},
)
def test_plot_line_mom0_offset(
    tmp_path: Path, sample_image_line: SampleImageFactoryVel
) -> Figure:
    return _plot_line_mom0_helper(tmp_path, sample_image_line, None, None, False, True)


def _plot_line_vlsr_helper(
    tmp_path: Path,
    sample_image_line: SampleImageFactoryVel,
    vmin: float | None,
    vmax: float | None,
    do_annotation: bool,
    do_offsets: bool,
) -> Figure:
    """Helper function for line moment 0 plotting tests."""
    dir = tmp_path
    dir.mkdir(exist_ok=True)
    file_name_vlsr = prodige_core.data_display.filename_line_vlsr(
        "B1-bS", "DCO+_3-2", mosaic=False
    )
    file_link_vlsr = os.path.join(os.fspath(dir), file_name_vlsr)
    hdu = sample_image_line(is_vlsr=True)
    print(r"max and min of vlsr data", np.nanmax(hdu.data), np.nanmin(hdu.data))
    hdu.writeto(file_link_vlsr, overwrite=True)

    file_name_TdV = prodige_core.data_display.filename_line_TdV(
        "B1-bS", "DCO+_3-2", mosaic=False
    )
    file_link_TdV = os.path.join(os.fspath(dir), file_name_TdV)
    hdu = sample_image_line(is_vlsr=False)
    rms = 0.1
    seed = 122807528840384100672342137672332424406
    rng = np.random.default_rng(seed)
    data = hdu.data + rng.standard_normal(hdu.data.shape) * rms
    hdu.data = data
    hdu.writeto(file_link_TdV, overwrite=True)

    with patch.dict(prodige_core.data_display.pyplot_params, {"text.usetex": False}):
        prodige_core.data_display.plot_line_vlsr(
            "B1-bS",
            "DCO+_3-2",
            os.fspath(dir) + "/",
            color_nan="0.9",
            mosaic=False,
            vmin=vmin,
            vmax=vmax,
            save_fig=False,
            do_marker=True,
            do_annotation=do_annotation,
            do_offsets=do_offsets,
            do_outflow=True,
        )
    return plt.gcf()


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_data_display",
    filename="example_line_velocity.png",
    tolerance=10,
    remove_text=True,
    savefig_kwargs={"dpi": IMAGE_COMPARISON_DPI},
)
def test_plot_line_vlsr(
    tmp_path: Path, sample_image_line: SampleImageFactoryVel
) -> Figure:
    return _plot_line_vlsr_helper(tmp_path, sample_image_line, 5.75, 7.75, True, False)


@pytest.mark.mpl_image_compare(
    baseline_dir="baseline_images/test_data_display",
    filename="example_line_velocity_offset.png",
    tolerance=10,
    remove_text=True,
    savefig_kwargs={"dpi": IMAGE_COMPARISON_DPI},
)
def test_plot_line_vlsr_offset(
    tmp_path: Path, sample_image_line: SampleImageFactoryVel
) -> Figure:
    return _plot_line_vlsr_helper(tmp_path, sample_image_line, None, None, False, True)
