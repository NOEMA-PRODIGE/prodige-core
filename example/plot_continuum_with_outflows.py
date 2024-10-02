import os
import prodige_core as pcore
# from prodige_core.config import pyplot_params

# NOEMA data directory
data_directory = os.getcwd() + '/'
# name of the region
region = 'L1448N'
# continuum baseband
bb = 'li'

pcore.plot_continuum(region, bb, data_directory, cmap='inferno',
                     color_nan='0.9', do_marker=True, do_outflow=True,)
