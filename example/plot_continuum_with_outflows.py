import os
# import numpy as np

# from astropy.coordinates import SkyCoord
# from astropy import units as u
# from astropy.io import fits
# from astropy.wcs import WCS

# import matplotlib.pyplot as plt
# from matplotlib.patches import Arrow
# from matplotlib.ticker import MultipleLocator
# from matplotlib import ticker

import prodige_core as pcore
from prodige_core.config import pyplot_params

# NOEMA data directory
data_directory = os.getcwd() + '/'

# name of the region
region = 'L1448N'

# distance to the region
# distance = 288.0  # pc

# continuum baseband
bb = 'li'

pcore.plot_continuum(region, bb, data_directory,
                     color_nan='0.9', do_marker=True)
