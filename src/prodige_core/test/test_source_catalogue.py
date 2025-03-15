from __future__ import annotations
import numpy as np
from contextlib import nullcontext as does_not_raise

from astropy import units as u
from astropy.io import fits
import prodige_core.source_catalogue
from prodige_core.source_catalogue import region_dic
import pytest


@pytest.mark.parametrize("source_id, expected_raise",
                         [
                             ("test", pytest.raises(ValueError)),
                             ("B1-bS", does_not_raise()),],)
def test_validate_source_id(source_id, expected_raise) -> None:
    with expected_raise:
        prodige_core.source_catalogue.validate_source_id(
            source_id) is not None


def test_get_outflow_information() -> None:
    sources_outflowPA, _, _, _, _ = (
        prodige_core.source_catalogue.get_outflow_information()
    )
    assert len(sources_outflowPA) == 76
    assert (sources_outflowPA[0] == "IRS3A") and (
        sources_outflowPA[-1] == "SVS13C")


def test_get_region_names() -> None:
    source_name = prodige_core.source_catalogue.get_region_names()
    # test that the first source is 'IRS3A'
    assert (source_name[0] == "L1448N") and (source_name[-1] == "SVS13B")
    assert len(source_name) == len(list(region_dic))


def test_load_cutout(sample_image) -> None:
    with pytest.raises(ValueError):
        prodige_core.source_catalogue.load_cutout("test.fits", source="test")
    # dir = tmp_path / "sub"
    # dir.mkdir()
    # file_link = dir / "test_image.fits"
    hdu_2d = sample_image(is_2d=True)
    ra0, dec0 = prodige_core.source_catalogue.get_region_center("B1-bS")
    hdu_3d = sample_image(is_2d=False)
    hdu_2d_new = prodige_core.source_catalogue.load_cutout(
        hdu_2d, source="B1-bS", is_hdu=True
    )
    hdu_new_3d = prodige_core.source_catalogue.load_cutout(
        hdu_3d, source="B1-bS", is_hdu=True
    )
    assert (hdu_2d_new.header["NAXIS1"] == 200) and (
        hdu_2d_new.header["NAXIS2"] == 200)
    assert (hdu_2d_new.header["CRVAL1"] == pytest.approx(ra0)) and (
        hdu_2d_new.header["CRVAL2"] == pytest.approx(dec0)
    )
    assert (hdu_2d_new.header) == (hdu_new_3d.header)


def test_get_region_center() -> None:
    ra0, dec0 = prodige_core.source_catalogue.get_region_center("L1448N")
    assert (ra0 == pytest.approx((3 + (25 + 36.44 / 60.0) / 60.0) * 15.0)) and (
        dec0 == pytest.approx(30 + (45 + 18.3 / 60.0) / 60.0)
    )


def test_get_figsize() -> None:
    # Default case should be 6.0 x 6.0
    fig_size = prodige_core.source_catalogue.get_figsize("Per-emb-2")
    assert fig_size == (6.0, 6.0)

@pytest.mark.parametrize("source_id, expected_raise",
                         [
                             ("test", pytest.raises(ValueError)),
                             ("B1-bS", does_not_raise()),],)
def test_get_region_vlsr(source_id, expected_raise) -> None:
    v_lsr = prodige_core.source_catalogue.get_region_vlsr("B1-bS")
    assert v_lsr == 6.75
    with expected_raise:
        prodige_core.source_catalogue.get_region_vlsr(
            source_id) is not None
