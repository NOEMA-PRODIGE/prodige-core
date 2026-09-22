# Licensed under a MIT style license - see LICENSE

"""
This is a package to handle PRODIGE data. This package relies in Astropy.
"""

from ._version import __version__
from .config import pyplot_params
from .data_display import (
    annotate_outflow,
    annotate_sources,
    determine_noise_map,
    get_contour_params,
    load_continuum_data,
    plot_continuum,
    plot_line_mom0,
    plot_line_vlsr,
)
from .source_catalogue import load_sources_table
