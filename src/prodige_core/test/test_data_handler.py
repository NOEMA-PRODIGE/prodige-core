import os
import numpy as np

from astropy.io import fits
from spectral_cube import SpectralCube
import prodige_core.data_handler

import pytest

@pytest.fixture
def fits_files(tmp_path):
    # Create temporary FITS files for testing
    filenames = []
    for i in range(3):
        data = np.random.random((10, 40, 40))
        hdu = fits.PrimaryHDU(data)
        hdu.header['CTYPE1'] = 'RA---TAN'
        hdu.header['CTYPE2'] = 'DEC--TAN'
        hdu.header['CTYPE3'] = 'VRAD'
        hdu.header['CUNIT1'] = 'deg'
        hdu.header['CUNIT2'] = 'deg'
        hdu.header['CUNIT3'] = 'km/s'
        hdu.header['CRPIX1'] = 20
        hdu.header['CRPIX2'] = 20
        hdu.header['CRPIX3'] = 5
        hdu.header['CRVAL1'] = 0.0
        hdu.header['CRVAL2'] = 0.0
        hdu.header['CRVAL3'] = 0.0
        hdu.header['CDELT1'] = -0.01 * (1 + i * 0.1)
        hdu.header['CDELT2'] = 0.01 * (1 + i * 0.1)
        hdu.header['CDELT3'] = 0.1
        hdu.header['BMAJ'] = 0.1 + i * 0.01
        hdu.header['BMIN'] = 0.1 + i * 0.01
        hdu.header['BPA'] = 45.0
        filename = tmp_path / f'test_{i}.fits'
        hdu.writeto(filename)
        filenames.append(str(filename))
    return filenames


def test_common_beam_files(fits_files):
    prodige_core.data_handler.common_beam_files(fits_files)
    # largest beam is the last cube
    last_smoothed_file = fits_files[-1].replace('.fits', '_smooth.fits')
    last_smoothed_cube = SpectralCube.read(last_smoothed_file)
    for f in fits_files:
        smoothed_file = f.replace('.fits', '_smooth.fits')
        assert os.path.exists(smoothed_file)

        original_cube = SpectralCube.read(f)
        smoothed_cube = SpectralCube.read(smoothed_file)

        assert smoothed_cube.shape == original_cube.shape
        assert smoothed_cube.beam == last_smoothed_cube.beam


def test_common_beam_files_invalid_extension(fits_files):
    invalid_files = [f.replace('.fits', '.txt') for f in fits_files]
    with pytest.raises(ValueError, match='All files must be FITS files.'):
        prodige_core.data_handler.common_beam_files(invalid_files)


def test_common_beam_files_nonexistent_file(fits_files):
    nonexistent_files = fits_files + ['nonexistent_file.fits']
    with pytest.raises(FileNotFoundError, match='File not found: nonexistent_file.fits'):
        prodige_core.data_handler.common_beam_files(nonexistent_files)


@pytest.fixture
def template_file(tmp_path):
    # Create a temporary template FITS file for testing
    data = np.random.random((10, 40, 40))
    hdu = fits.PrimaryHDU(data)
    hdu.header['CTYPE1'] = 'RA---TAN'
    hdu.header['CTYPE2'] = 'DEC--TAN'
    hdu.header['CTYPE3'] = 'VRAD'
    hdu.header['CUNIT1'] = 'deg'
    hdu.header['CUNIT2'] = 'deg'
    hdu.header['CUNIT3'] = 'km/s'
    hdu.header['CRPIX1'] = 20
    hdu.header['CRPIX2'] = 20
    hdu.header['CRPIX3'] = 5
    hdu.header['CRVAL1'] = 0.0
    hdu.header['CRVAL2'] = 0.0
    hdu.header['CRVAL3'] = 0.0
    hdu.header['CDELT1'] = -0.01
    hdu.header['CDELT2'] = 0.01
    hdu.header['CDELT3'] = 0.1
    hdu.header['NAXIS1'] = 40
    hdu.header['NAXIS2'] = 40
    filename = tmp_path / 'template.fits'
    hdu.writeto(filename)
    return str(filename)


def test_regrid_cubes_from_files(fits_files, template_file):
    prodige_core.data_handler.regrid_cubes_from_files(
        fits_files, fits_files[-1])
    for f in fits_files:
        regridded_file = f.replace('.fits', '_regrid.fits')
        assert os.path.exists(regridded_file)

        original_cube = SpectralCube.read(f)
        regridded_cube = SpectralCube.read(regridded_file)

        assert regridded_cube.shape == original_cube.shape
        assert regridded_cube.header['CRVAL1'] == fits.getheader(template_file)[
            'CRVAL1']
        assert regridded_cube.header['CRVAL2'] == fits.getheader(template_file)[
            'CRVAL2']


def test_regrid_cubes_from_files_invalid_extension(fits_files, template_file):
    invalid_files = [f.replace('.fits', '.txt') for f in fits_files]
    with pytest.raises(ValueError, match='All files must be FITS files.'):
        prodige_core.data_handler.regrid_cubes_from_files(
            invalid_files, template_file)


def test_regrid_cubes_from_files_nonexistent_file(fits_files, template_file):
    nonexistent_files = fits_files + ['nonexistent_file.fits']
    with pytest.raises(FileNotFoundError, match='File not found: nonexistent_file.fits'):
        prodige_core.data_handler.regrid_cubes_from_files(
            nonexistent_files, template_file)
