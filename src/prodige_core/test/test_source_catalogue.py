import prodige_core.source_catalogue

from prodige_core.source_catalogue import region_dic
import pytest


def test_get_outflow_information_number_source() -> None:
    list_source_info = prodige_core.source_catalogue.get_outflow_information()
    # test that there are 16 sources in the region dictionary
    assert len(list_source_info[0]) == 16


def test_get_outflow_information_first_last_sources() -> None:
    sources_outflowPA, _, _, _, _ = prodige_core.source_catalogue.get_outflow_information()
    # print(source_name[0])
    # test that the first source is 'IRS3A'
    assert (sources_outflowPA[0] == 'IRS3A') and (
        sources_outflowPA[-1] == 'B1-bS')


def test_get_region_names_first_last_sources() -> None:
    source_name = prodige_core.source_catalogue.get_region_names()
    # test that the first source is 'IRS3A'
    assert (source_name[0] == 'L1448N') and (source_name[-1] == 'B1-bS')


def test_get_region_names_number_source() -> None:
    source_name = prodige_core.source_catalogue.get_region_names()
    # test that there are 16 sources in the region dictionary
    assert len(source_name) == len(list(region_dic))


def test_load_cutout_nofile() -> None:
    with pytest.raises(ValueError):
        prodige_core.source_catalogue.load_cutout('test.fits', source='test')
