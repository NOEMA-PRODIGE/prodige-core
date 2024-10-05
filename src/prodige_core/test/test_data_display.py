from __future__ import annotations

import prodige_core.data_display
from astropy import units as u

import pytest


def test_pb_telecope_good_frequency() -> None:
    assert prodige_core.data_display.pb_telecope(72.78382*u.GHz, telescope='NOEMA') == 64.1 * u.arcsec \
        and prodige_core.data_display.pb_telecope(345*u.GHz, telescope='SMA') == 36.0 * u.arcsec \
        and prodige_core.data_display.pb_telecope(1*u.GHz, telescope='VLA') == 45.0 * u.arcmin \
        and prodige_core.data_display.pb_telecope(300*u.GHz, telescope='ALMA') == 19.0*u.arcsec
    
def test_pb_telecope_bad_telescope() -> None:
    with pytest.raises(ValueError):
        prodige_core.data_display.pb_telecope(72.78382 * u.GHz, telescope='TEST')

def test_validate_frequency_bad_frequency() -> None:
    with pytest.raises(u.UnitsError):
        prodige_core.data_display.validate_frequency(72.78382 * u.m)

def test_validate_frequency_good_frequency() -> None:
    assert prodige_core.data_display.validate_frequency(72.78382 * u.GHz)
