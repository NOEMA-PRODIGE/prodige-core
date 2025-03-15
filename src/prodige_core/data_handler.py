import radio_beam
from astropy.io import fits
from astropy import units as u
from spectral_cube import SpectralCube
import glob


def check_fits_files(fits_files: list) -> None:
    """
    Check that all files in the list are FITS files and are present in the path.

    Parameters:
    fits_files (list of str): List of paths to FITS files.

    Returns:
        None
    """
    for f in fits_files:
        if not f.endswith('.fits'):
            raise ValueError('All files must be FITS files.')
        if not glob.glob(f):
            raise FileNotFoundError(f'File not found: {f}')

    return


def common_beam_files(fits_files: list, suffix: str = '_smooth') -> None:
    """
    Load a list of FITS files, find the common beam using radio_beam. 
    The smoothed cubes are saved with the suffix '_smooth' added to the original 
    file name and before the '.fits' extension.

    Parameters:
    fits_files (list of str): List of paths to FITS files.
    suffix (str): Suffix to add to the output file names. Default is '_smooth'.

    Returns:
        None
    """
    # check that all files are FITS files and are present in the path
    check_fits_files(fits_files)

    # load all cubes and find the common beam
    cubes = [SpectralCube.read(f) for f in fits_files]
    # make list of BMAJ, BMIN, and BPA to create a Beams object
    bmaj_list, bmin_list, bpa_list = [], [], []
    for cube in cubes:
        bmaj_list.append(cube.header['BMAJ'])
        bmin_list.append(cube.header['BMIN'])
        bpa_list.append(cube.header['BPA'])

    my_beams = radio_beam.Beams(
        major=bmaj_list * u.deg, minor=bmin_list * u.deg, pa=bpa_list * u.deg)
    common_beam = my_beams.common_beam()

    convolved_cubes = [cube.convolve_to(common_beam) for cube in cubes]
    # write out the smoothed cubes
    for i, new_cube in enumerate(convolved_cubes):
        new_cube.write(fits_files[i].replace(
            '.fits', f'{suffix}.fits'), overwrite=True)

    return


def regrid_cubes_from_files(fits_files: list, template_file: str, suffix: str = '_regrid') -> None:
    """
    Load a list of FITS files, regrid them onto a common grid defined by a template file.
    The regridded cubes are saved with the suffix '_regrid' added to the original 
    file name and before the '.fits' extension.

    Parameters:
    fits_files (list of str): List of paths to FITS files.
    template_file (str): Path to the template FITS file defining the common grid.
    suffix (str): Suffix to add to the output file names. Default is '_regrid'.

    Returns:
        None
    """
    # check that all files are FITS files and are present in the path
    check_fits_files(fits_files)
    check_fits_files([template_file])

    # load the template cube
    template_hd = fits.getheader(template_file)
    list_key = ['CRVAL1', 'CRVAL2', 'CDELT1', 'CDELT2', 'CRPIX1', 'CRPIX2',
                'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2', 'NAXIS1', 'NAXIS2']
    for key in list_key:
        if key not in template_hd:
            raise KeyError(f'Header key {key} not found in the template file.')
    # load all cubes and regrid them
    for f in fits_files:
        cube = SpectralCube.read(f)
        target_hd = cube.header
        for key in list_key:
            target_hd[key] = template_hd[key]
        regridded_cube = cube.reproject(target_hd)
        regridded_cube.write(
            f.replace('.fits', f'{suffix}.fits'), overwrite=True)

    return
