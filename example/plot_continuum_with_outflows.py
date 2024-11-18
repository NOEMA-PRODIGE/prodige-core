import os
import prodige_core as pcore
# from prodige_core.config import pyplot_params

# NOEMA data directory
data_directory = os.getcwd() + '/'

# # name of the region
region = 'L1448N'
# continuum baseband
bb = 'li'
print('>>>>> ' + region + ' continuum of Mosaic')
pcore.plot_continuum(region, bb, data_directory, cmap='inferno', mosaic=True,
                     color_nan='0.9', do_marker=True, do_outflow=True,)

# name of the region
region = 'B5-IRS1'
# continuum baseband
bb = 'lo'
print('>>>>> ' + region + ' continuum of single pointing')
pcore.plot_continuum(region, bb, data_directory, cmap='inferno', mosaic=False,
                     color_nan='0.9', do_marker=True, do_outflow=True,)


# # name of the region
region = 'HH211'
# # linename
linename = 'N2Dp_K'
bb = 'li'
print('>>>>> ' + region + ' Integrated Intensity')
pcore.plot_line_mom0(region, linename, bb, data_directory, mosaic=False,
                     color_nan='0.9', do_marker=True, do_outflow=True,)

print('>>>>> ' + region + ' Velocity map')
pcore.plot_line_vlsr(region, linename, data_directory, mosaic=False,
                     color_nan='0.9', do_marker=True, do_outflow=True, do_offsets=True)
