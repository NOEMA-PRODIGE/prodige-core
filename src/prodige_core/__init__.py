# Licensed under a MIT style license - see LICENSE

"""
This is a package to handle PRODIGE data. This package relies in Astropy.
"""

import importlib.metadata

from .config import pyplot_params
from .data_display import (
    annotate_outflow,
    annotate_sources,
    determine_noise_map,
    get_contour_params,
    LinePanelSpec,
    load_continuum_data,
    load_line_vlsr,
    plot_continuum,
    plot_line_mom0,
    plot_line_vlsr,
    plot_lines_for_source,
    plot_source_line_grid,
    plot_source_line_pairs,
    plot_sources_for_line,
)
from .source_catalogue import load_sources_table

__all__ = [
    "pyplot_params",
    "annotate_outflow",
    "annotate_sources",
    "determine_noise_map",
    "get_contour_params",
    "LinePanelSpec",
    "load_continuum_data",
    "load_line_vlsr",
    "plot_continuum",
    "plot_line_mom0",
    "plot_line_vlsr",
    "plot_lines_for_source",
    "plot_source_line_grid",
    "plot_source_line_pairs",
    "plot_sources_for_line",
    "load_sources_table",
]
try:
    from ._version import __version__
except ImportError:  # pragma: no cover
    __version__ = importlib.metadata.version("prodige-core")  # pragma: no cover
