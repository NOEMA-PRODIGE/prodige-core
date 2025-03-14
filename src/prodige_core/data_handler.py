import radio_beam
# from astropy.io import fits
from astropy import units as u
from spectral_cube import SpectralCube
import glob


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
    for f in fits_files:
        if not glob.glob(f):
            raise FileNotFoundError(f'File not found: {f}')
        if not f.endswith('.fits'):
            raise ValueError('All files must be FITS files.')

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
