# Licensed under a MIT style license - see LICENSE

"""
This is a package to handle PRODIGE data. This package relies in Astropy.
"""

from .data_display import plot_continuum, plot_line_mom0, plot_line_vlsr, get_contour_params, load_continuum_data, determine_noise_map, annotate_sources, annotate_outflow
from .source_catalogue import load_sources_table

from .config import pyplot_params
