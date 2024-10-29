from __future__ import annotations
import numpy as np

from astropy import units as u
from astropy.io import fits
import prodige_core.source_catalogue
from prodige_core.source_catalogue import region_dic
import pytest


def test_validate_source_id() -> None:
    with pytest.raises(ValueError):
        prodige_core.source_catalogue.validate_source_id("test")
    # def test_validate_source_id_source() -> None:
    # with pytest.raises(ValueError):
    assert prodige_core.source_catalogue.validate_source_id("B1-bS")


def test_get_outflow_information() -> None:
    sources_outflowPA, _, _, _, _ = (
        prodige_core.source_catalogue.get_outflow_information()
    )
    assert len(sources_outflowPA) == 76
    assert (sources_outflowPA[0] == "IRS3A") and (sources_outflowPA[-1] == "SVS13C")


def test_get_region_names_first() -> None:
    source_name = prodige_core.source_catalogue.get_region_names()
    # test that the first source is 'IRS3A'
    assert (source_name[0] == "L1448N") and (source_name[-1] == "SVS13B")
    assert len(source_name) == len(list(region_dic))


def test_load_cutout() -> None:
    with pytest.raises(ValueError):
        prodige_core.source_catalogue.load_cutout("test.fits", source="test")
    # dir = tmp_path / "sub"
    # dir.mkdir()
    # file_link = dir / "test_image.fits"
    data = np.ones((101, 101))
    data2 = data[np.newaxis, :, :]
    ra0, dec0 = prodige_core.source_catalogue.get_region_center("B1-bS")
    hdu = fits.PrimaryHDU(data=data)
    hdu2 = fits.PrimaryHDU(data=data2)
    hdu.header["CRVAL1"] = ra0
    hdu.header["CRVAL2"] = dec0
    hdu.header["CRPIX1"] = 51
    hdu.header["CRPIX2"] = 51
    hdu.header["CDELT1"] = 40.0 * u.arcsec.to(u.deg) / 20
    hdu.header["CDELT2"] = 40.0 * u.arcsec.to(u.deg) / 20
    hdu.header["CUNIT1"] = "deg"
    hdu.header["CUNIT2"] = "deg"
    hdu.header["CTYPE1"] = "RA---TAN"
    hdu.header["CTYPE2"] = "DEC--TAN"
    hdu2.header = hdu.header.copy()
    # hdu.writeto(file_link, overwrite=True)
    hdu_new = prodige_core.source_catalogue.load_cutout(
        hdu, source="B1-bS", is_hdu=True
    )
    hdu_new2 = prodige_core.source_catalogue.load_cutout(
        hdu2, source="B1-bS", is_hdu=True
    )
    assert (hdu_new.header["NAXIS1"] == 20) and (hdu_new.header["NAXIS2"] == 20)
    assert (hdu_new.header["CRVAL1"] == pytest.approx(ra0)) and (
        hdu_new.header["CRVAL2"] == pytest.approx(dec0)
    )
    assert (hdu_new.header) == (hdu_new2.header)


def test_get_region_center() -> None:
    ra0, dec0 = prodige_core.source_catalogue.get_region_center("L1448N")
    assert (ra0 == pytest.approx((3 + (25 + 36.44 / 60.0) / 60.0) * 15.0)) and (
        dec0 == pytest.approx(30 + (45 + 18.3 / 60.0) / 60.0)
    )


def test_get_figsize() -> None:
    # Default case should be 6.0 x 6.0
    fig_size = prodige_core.source_catalogue.get_figsize("Per-emb-2")
    assert fig_size == (6.0, 6.0)
