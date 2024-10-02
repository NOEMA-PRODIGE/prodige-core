import numpy as np

from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

from importlib.resources import files
from .config import source_filename, distance

# name of the region
# region = 'L1448N'
# distance to the region
# distance = 288.0  # pc


def load_sources_table():
    """
    Load the source table containing the sources within the field of view.
    It uses the source_filename variable from the config.py file.
    """
    # load table containing sources within FoV
    # sources within FOV
    data_file = files('prodige_core').joinpath(source_filename)
    # data_file = files('prodige_core').joinpath('sources.dat')  # .read_text()
    # sources_tab = np.loadtxt('sources.dat',dtype='U')
    sources_tab = np.loadtxt(data_file, dtype='U')

    # source name
    name_list = sources_tab[:, 0]
    # source RA
    RA_list = sources_tab[:, 1]
    # source DEC
    Dec_list = sources_tab[:, 2]
    # source color
    color_list = sources_tab[:, 3]
    # source vlsr
    vlsr_list = sources_tab[:, 4].astype(float)
    # source outflow PA
    outflowPA = sources_tab[:, 5].astype(float)

    return name_list, RA_list, Dec_list, color_list, vlsr_list, outflowPA


# regions dictionary storing the region name, the central RA and DEC (FK5), and
# the size of the region (height and width wiht units)
# these are supposed to be used for plotting purposes
region_dic = {
    'L1448N': {'RA0': '3:25:36.44', 'Dec0': '30:45:18.3',
               'height': 33*u.arcsec, 'width': 30*u.arcsec,
               'fig_width': 6, 'fig_height': 6},
    'B5-IRS1': {'RA0': '3:25:36.44', 'Dec0': '30:45:18.3',
                'height': 33*u.arcsec, 'width': 30*u.arcsec,
                'fig_width': 6, 'fig_height': 6},
}


def load_cutout(file_in, source='L1448N', is_hdu=False):
    """
    Convenience function to load a FITS file, or an existing HDU,
    and to generate a cutout following the requested center and sizes.

    params:

    file_in : FITS file name to be loaded (if is_hdu=False). If is_hdu=True,
              the file_in represents a valid HDU to be used.
    source : String parameter. It represents the source name to be used to
             define the center of the cutout. This will be loaded from the
             region dictionary.
    is_hdu : Boolean parameter. It controls what is passed to the function.

    Return:
    It returns the HDU of the cutout defined by the global variables: position and cutout_size.
    """
    if source not in region_dic.keys():
        raise ValueError('Source not in region dictionary')
    else:
        # get the region center
        position = SkyCoord(region_dic[source]['RA0'] + ' ' +
                            region_dic[source]['Dec0'], unit=(u.hourangle, u.deg))
        cutout_size = u.Quantity((region_dic[source]['height'],
                                  region_dic[source]['width']))

    if is_hdu == False:
        hdu = fits.open(file_in)[0]
    else:
        hdu = file_in.copy()
    # Make the cutout, including the WCS
    wcs_original = WCS(hdu.header).dropaxis(2)
    cutout = Cutout2D(np.squeeze(hdu.data), position=position,
                      size=cutout_size, wcs=wcs_original)
    # update HDU with new data and updated header information
    hdu.data = cutout.data
    hdu.header.update(cutout.wcs.to_header())
    # list of keys to be removed (if they exist)
    key_list = ['CTYPE3', 'CRVAL3', 'CRPIX3',
                'CDELT3', 'CUNIT3', 'NAXIS3', 'CROTA3', 'SPECSYS', 'VELREF']
    for key in key_list:
        if key in hdu.header:
            del hdu.header[key]
    return hdu


def get_figsize(source):
    """
    Convenience function to get the figure size for a given source.
    """
    if source not in region_dic.keys():
        raise ValueError('Source not in region dictionary')
    else:
        fig_width = region_dic[source]['fig_width']
        fig_height = region_dic[source]['fig_height']
    return fig_width, fig_height
