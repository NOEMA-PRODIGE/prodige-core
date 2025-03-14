import os
import numpy as np
import pytest
from astropy.io import fits
from spectral_cube import SpectralCube
from prodige_core.data_handler import common_beam_files


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
        hdu.header['CRPIX1'] = 20
        hdu.header['CRPIX2'] = 20
        hdu.header['CRPIX3'] = 5
        hdu.header['CRVAL1'] = 0.0
        hdu.header['CRVAL2'] = 0.0
        hdu.header['CRVAL3'] = 0.0
        hdu.header['CDELT1'] = -0.01
        hdu.header['CDELT2'] = 0.01
        hdu.header['CDELT3'] = 0.1
        hdu.header['BMAJ'] = 0.1 + i * 0.01
        hdu.header['BMIN'] = 0.1 + i * 0.01
        hdu.header['BPA'] = 45.0
        filename = tmp_path / f'test_{i}.fits'
        hdu.writeto(filename)
        filenames.append(str(filename))
    return filenames


def test_common_beam_files(fits_files):
    common_beam_files(fits_files)

    for f in fits_files:
        smoothed_file = f.replace('.fits', '_smooth.fits')
        assert os.path.exists(smoothed_file)

        original_cube = SpectralCube.read(f)
        smoothed_cube = SpectralCube.read(smoothed_file)

        assert smoothed_cube.shape == original_cube.shape
        assert smoothed_cube.beam == smoothed_cube.beam.common_beam()
