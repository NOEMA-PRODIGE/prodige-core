from __future__ import annotations

from typing import cast

import matplotlib
import matplotlib.patheffects as PathEffects
from matplotlib.colors import Colormap

# try:
#     _ = matplotlib.get_backend()
# except Exception:
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.stats import sigma_clipped_stats
from astropy.visualization.wcsaxes import SphericalCircle, add_beam, add_scalebar
from astropy.wcs import WCS
from matplotlib import ticker
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.image import AxesImage
from numpy.typing import NDArray

from .config import (
    cmap_default,
    cmap_mom0_default,
    cmap_vlsr_default,
    distance,
    pyplot_params,
)
from .source_catalogue import (
    get_figsize,
    get_outflow_information,
    get_region_center,
    get_region_vlsr,
    load_cutout,
    load_sources_table,
)

# name of the region
label_col = "black"
label_col_back = "white"


def determine_noise_map(data_2d: NDArray[np.float64]) -> float:
    """
    Determine the noise in the continuum data.
    """
    # Ensure input is numeric before passing to Astropy
    if not isinstance(data_2d, np.ndarray) or not np.issubdtype(
        data_2d.dtype, np.number
    ):
        raise ValueError("Input data must be a numeric numpy array.")
    # compute noise in continuum data
    _, _, noise_2dmap = sigma_clipped_stats(data_2d, sigma=3.0)
    return noise_2dmap


def get_contour_params(
    maximum: float, noise: float
) -> tuple[NDArray[np.float64], list[str], bool]:
    """
    Compute the contour levels for the continuum data.
    maximum: maximum value to be shown in the plot
    noise: noise in the data
    """
    # compute contour levels at -5,5,10,20,40,80x... sigma
    # determines the number of contours to be plotted
    steps = int(np.log(maximum / (5.0 * noise)) // np.log(2.0)) + 1
    if steps < 1:
        return np.array([0.0], dtype=np.float64), ["solid"], False
    steps_arr = np.logspace(
        start=0,
        stop=steps,
        num=steps,
        endpoint=False,
        base=2.0,
        dtype=np.float64,
        axis=0,
    )
    # append -5 sigma to the array and multiply by step size
    steps_arr = np.append(-steps_arr[0], steps_arr) * 5.0 * noise
    line_styles = ["dotted"] + ["solid"] * steps
    return steps_arr, line_styles, True


def default_minmax(
    data: NDArray[np.float64],
    noise: float,
    vmin: float | None,
    vmax: float | None,
    vmax_scale: float = 1.0,
) -> tuple[float, float]:
    """
    Fill in vmin/vmax defaults shared by the continuum and moment-0 plots:
    vmin defaults to -5*noise, vmax defaults to vmax_scale*max(data).
    """
    if vmin is None:
        vmin = -5.0 * noise
    if vmax is None:
        vmax = vmax_scale * np.nanmax(data)
    return vmin, vmax


def plot_data_contours(
    ax: Axes,
    data: NDArray[np.float64],
    wcs: WCS,
    noise: float,
    linewidth_white: float = 0.75,
    linewidth_black: float = 0.35,
) -> bool:
    """
    Draw the white+black sigma contour levels shared by all panel plots.
    Returns whether any contour level was valid (and thus drawn).
    """
    cont_levels, style_levels, valid_contour = get_contour_params(
        np.nanmax(data), noise
    )
    if valid_contour:
        ax.contour(
            data,
            colors="white",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=linewidth_white,
            transform=ax.get_transform(wcs),
        )
        ax.contour(
            data,
            colors="black",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=linewidth_black,
            transform=ax.get_transform(wcs),
        )
    return valid_contour


def filename_continuum(region: str, bb: str, mosaic: bool = False) -> str:
    """Function to return the filename of the continuum data.
    It follows the naming convention of PRODIGE.
    Parameters:
    region: name of the region
    bb: baseband of the data (lo, li, ui, uo)
    mosaic: if True, mosaic data is used. This changes the filename format of the data.
    """
    if mosaic:
        datafile = region + "_CD_" + bb + "_cont_rob1-selfcal.fits"
    else:
        datafile = region + "_CD_" + bb + "_cont_rob1-selfcal-pbcor.fits"
    return datafile


def filename_line_TdV(region: str, linename: str, mosaic: bool = False) -> str:
    """Function to return the filename of the line integrated intensity data.
    It follows the naming convention of PRODIGE.
    Parameters:
    region: name of the region
    linename: line name (e.g., 'H2CO_l21')
    mosaic: if True, mosaic data is used. This changes the filename format of the data.
    """
    if mosaic:
        datafile = region + "_CD_" + linename + "_TdV.fits"
    else:
        datafile = region + "_CD_" + linename + "_TdV.fits"
    return datafile


def filename_line_vlsr(region: str, linename: str, mosaic: bool = False) -> str:
    """Function to return the filename of the line centroid velocity.
    It follows the naming convention of PRODIGE.
    Parameters:
    region: name of the region
    linename: line name (e.g., 'H2CO_l21')
    mosaic: if True, mosaic data is used. This changes the filename format of the data.
    """
    if mosaic:
        datafile = region + "_CD_" + linename + "_Vlsr.fits"
    else:
        datafile = region + "_CD_" + linename + "_Vlsr.fits"
    return datafile


def load_continuum_data(
    datafile: str,
    region: str,
) -> tuple[NDArray[np.float64], float, Header]:
    """
    Function to load the continuum data and return the cutout specified in the dictionary.
    It return the data, estimated noise, and the FITS header.
    Parameters:
    datafile: fileneame of the data to load
    region: name of the region

    Returns:
    data_cont: continuum data in mJy/beam
    noise_cont: estimated noise in the continuum data
    """

    # loads the cutout of the region. It uses the region dictionary to set the cutout size.
    hdu_cont = load_cutout(datafile, source=region, is_hdu=False)
    # set empty pixels (0.0) to NaN

    header = cast(Header, hdu_cont.header)
    data_cont = cast(NDArray, hdu_cont.data)
    data_cont[data_cont == 0.0] = np.nan
    # Update the header with the updated WCS from the cutout, as well as the data in mJy/beam.
    # header = hdu_cont.header
    if str(header["BUNIT"]).casefold() == "JY/BEAM".casefold():
        data_cont = np.squeeze(data_cont) * 1e3
        header["BUNIT"] = "mJy/beam"
    else:
        data_cont = np.squeeze(data_cont)
    # compute noise
    noise_cont = determine_noise_map(data_cont)
    return data_cont, noise_cont, header


def load_line_TdV(
    datafile: str,
    region: str,
) -> tuple[NDArray[np.float64], float, Header]:
    """
    Function to load the integrated intensity map and return the cutout specified in the dictionary.
    It return the data, estimated noise, and the FITS header.
    Parameters:
    datafile: fileneame of the data to load
    region: name of the region

    Returns:
    data_TdV: integrated intensity map in mJy/beam km/s or K km/s
    noise_map: estimated noise from the data
    """
    # loads the cutout of the region. It uses the region dictionary to set the cutout size.
    hdu = load_cutout(datafile, source=region, is_hdu=False)
    # set empty pixels (0.0) to NaN
    data = cast(NDArray, hdu.data)
    data = np.squeeze(data)
    data[data == 0.0] = np.nan
    # Update the header with the updated WCS from the cutout, as well as the data in mJy/beam.
    header = hdu.header
    unit_string = str(header["BUNIT"]).casefold()
    if unit_string == "JY/BEAM KM/S".casefold():
        data = data * 1e3
        header["BUNIT"] = "mJy/beam km/s"
    elif unit_string == "mJY/BEAM KM/S".casefold():
        header["BUNIT"] = "mJy/beam km/s"
    else:
        header["BUNIT"] = "K km/s"
    # compute noise
    noise_map = determine_noise_map(data)
    return data, noise_map, header


def get_frequency(header: Header) -> float:
    """
    Function to get the frequency from the header and convert it to GHz.
    Parameters:
    header: FITS header

    Returns:
    frequency: frequency in GHz
    """
    if "RESTFREQ" not in header:
        raise ValueError("RESTFREQ not found in header.")

    restfreq = float(header["RESTFREQ"])  # type: ignore[arg-type]
    return restfreq * 1e-9


def get_wavelength(header: Header) -> float:
    """
    Function to get the wavelength from the header.
    Parameters:
    header: FITS header

    Returns:
    wavelength: wavelength in mm
    """
    restfreq_ghz = get_frequency(header)
    # Perform unit conversion
    freq_quantity = u.Quantity(restfreq_ghz, u.GHz)
    wavelength_quantity: u.Quantity = freq_quantity.to(u.mm, equivalencies=u.spectral())  # type: ignore[reportUnknownMemberType]
    return round(float(wavelength_quantity.value), 1)  # type: ignore[reportUnknownMemberType]


def prodige_style(
    ax: Axes, do_offsets: bool = False, center_coord: SkyCoord | None = None
) -> None:
    """
    Setting a common style for the plots. This includes axis labels, tick labels, and minor ticks.
    Pararameters:
    ax: axis object.
    """
    # plot properties
    if do_offsets == False:
        RA = ax.coords[0]
        DEC = ax.coords[1]
        RA.set_axislabel(r"$\alpha$ (J2000)", minpad=0.7)
        DEC.set_major_formatter("dd:mm:ss")
        RA.set_major_formatter("hh:mm:ss.s")
        DEC.set_axislabel(r"$\delta$ (J2000)", minpad=0.8)
        DEC.set_ticklabel(rotation=90.0, color="black", exclude_overlapping=True)
        RA.set_ticklabel(color="black", exclude_overlapping=True)
        DEC.set_ticks(spacing=10 * u.arcsec, color="black")  # type: ignore
        RA.set_ticks(spacing=1.0 * 15 * u.arcsec, color="black")  # type: ignore
        RA.display_minor_ticks(True)
        DEC.display_minor_ticks(True)
        DEC.set_minor_frequency(5)
        RA.set_minor_frequency(5)
    else:
        if center_coord is None:
            raise ValueError("Center coordinate is not defined.")
        # Using implementation from
        # https://community.openastronomy.org/t/maps-in-relative-coordinates-with-wcsaxes/186/4
        RA = ax.coords[0]
        DEC = ax.coords[1]
        RA.set_ticks_visible(False)
        RA.set_ticklabel_visible(False)
        DEC.set_ticks_visible(False)
        DEC.set_ticklabel_visible(False)
        RA.set_axislabel("")
        DEC.set_axislabel("")

        off_frame = center_coord.skyoffset_frame()
        overlay_coord = ax.get_coords_overlay(off_frame)
        ra_offset = overlay_coord["lon"]
        dec_offset = overlay_coord["lat"]
        ra_offset.set_axislabel("R.A. offset")
        dec_offset.set_axislabel("Dec. offset")
        ra_offset.set_major_formatter("s")
        dec_offset.set_major_formatter("s")
        ra_offset.set_ticks_position("bt")
        ra_offset.set_ticklabel_position("b")
        dec_offset.set_ticks_position("lr")
        dec_offset.set_ticklabel_position("l")
        ra_offset.set_axislabel_position("b")
        dec_offset.set_axislabel_position("l")
        ra_offset.coord_wrap = 180 * u.deg  # avoid wrapping # type: ignore
        ra_offset.display_minor_ticks(True)
        dec_offset.display_minor_ticks(True)
        dec_offset.set_minor_frequency(5)
        ra_offset.set_minor_frequency(5)
        dec_offset.set_ticks(spacing=15 * u.arcsec, color="black")  # type: ignore
        ra_offset.set_ticks(spacing=15 * u.arcsec, color="black")  # type: ignore
        # remember the overlay coords, since they (not ax.coords) carry the visible labels
        ax._offset_coords = (ra_offset, dec_offset)


# @u.quantity_input
def annotate_sources(
    ax: Axes,
    wcs: WCS,
    color: str = "cornflowerblue",
    color_back: str = "black",
    marker: bool = False,
    label: bool = True,
    connect_line: bool = False,
    fontsize: int = 10,
    label_offset: u.Quantity = 1.0 * u.arcsec,  # type: ignore
) -> None:
    """
    Convenience function to annotate sources in the field of view.
    Parameters:
    ax: axis object
    wcs: WCS object
    color: color of the text
    color_back: color of the edge around the text (for better visibility)
    marker: if True, a marker is added to the source position usign the
    coordinates from the dictionary
    label: if True, the source name is added to the plot
    fontsize: fontsize of the text
    label_offset: offset of the labels
    """
    # load table containing sources within the region
    sources_name, sources_RA, sources_Dec, _, _, _, label_offsetPA = (
        load_sources_table()
    )
    # loop over all labels
    for source_i, RA_i, Dec_i, offset_PA_i in zip(
        sources_name, sources_RA, sources_Dec, label_offsetPA
    ):
        c = SkyCoord(ra=RA_i, dec=Dec_i, unit=(u.hourangle, u.deg))  # type: ignore
        # Check if source is within the field of view
        if not wcs.footprint_contains(c):  # type: ignore[reportUnknownMemberType]
            continue
        if marker:
            ax.scatter(
                c.ra,  # type: ignore
                c.dec,  # type: ignore
                marker="*",
                c=color,
                edgecolor="black",
                linewidth=0.5,
                s=20,
                transform=ax.get_transform("world"),
                zorder=40,
            )

        if label:
            c_label = c.directional_offset_by(offset_PA_i * u.deg, label_offset)  # type: ignore
            label_text = ax.text(
                c_label.ra.degree,  # type: ignore
                c_label.dec.degree,  # type: ignore
                r"\textbf{" + str(source_i) + r"}",
                transform=ax.get_transform("world"),
                color=color,
                fontsize=fontsize,
                verticalalignment="center",
                horizontalalignment="center",
            )
            label_text.set_path_effects(
                [PathEffects.withStroke(linewidth=1.0, foreground=color_back)]
            )

        if connect_line:
            c_line_start = c.directional_offset_by(
                offset_PA_i * u.deg,
                0.2 * label_offset,  # type: ignore
            )
            c_line_end = c.directional_offset_by(
                offset_PA_i * u.deg,
                0.5 * label_offset,  # type: ignore
            )
            ax.plot(
                [c_line_start.ra.degree, c_line_end.ra.degree],  # type: ignore
                [c_line_start.dec.degree, c_line_end.dec.degree],  # type: ignore
                color=color_back,
                lw=1.5,
                alpha=0.7,
                zorder=10,
                transform=ax.get_transform("world"),
            )
            ax.plot(
                [c_line_start.ra.degree, c_line_end.ra.degree],  # type: ignore
                [c_line_start.dec.degree, c_line_end.dec.degree],  # type: ignore
                color=color,
                lw=1.0,
                alpha=0.7,
                zorder=10,
                transform=ax.get_transform("world"),
            )


# @u.quantity_input
# @u.quantity_input(arrow_length="angle", arrow_offset="angle")
def annotate_outflow(
    ax: Axes,
    wcs: WCS,
    arrow_width: float = 1.0,
    arrow_length: u.Quantity = 3 * u.arcsec,  # type: ignore
    arrow_offset: u.Quantity = 0.05 * u.arcsec,  # type: ignore
) -> None:
    """
    Function to add outflow orientations to the plot.
    Parameters:
    ax: axis object
    wcs: WCS object
    arrow_width: width of the arrows
    arrow_length: length of the arrows
    arrow_offset: offset of the arrows
    """
    # add outflow orientation angle
    default_width = 0.000025
    default_head_width = 0.000075
    # load table containing sources within the region
    _, sources_RA, sources_Dec, sources_outflowPA, _ = get_outflow_information()
    # loop over all cores
    for RA_i, Dec_i, source_outflowPA_i in zip(
        sources_RA, sources_Dec, sources_outflowPA
    ):
        # for k in range(sources_name.size):
        # source coordinate
        c = SkyCoord(ra=RA_i, dec=Dec_i, unit=(u.hourangle, u.deg))  # type: ignore
        #  sources_Dec[k], unit=(u.hourangle, u.deg))
        # check if source is within the field of view
        # and if the outflow orientation is defined
        if not (wcs.footprint_contains(c) & np.isfinite(source_outflowPA_i)):
            continue
        c_blue_start = c.directional_offset_by(source_outflowPA_i * u.deg, arrow_offset)  # type: ignore
        c_blue_end = c.directional_offset_by(source_outflowPA_i * u.deg, arrow_length)  # type: ignore
        c_red_start = c.directional_offset_by(
            (180 + source_outflowPA_i) * u.deg,
            arrow_offset,  # type: ignore
        )
        c_red_end = c.directional_offset_by(
            (180 + source_outflowPA_i) * u.deg,
            arrow_length,  # type: ignore
        )
        # calculate the offset for the arrows
        dx_blue = c_blue_end.ra.degree - c_blue_start.ra.degree  # type: ignore
        dy_blue = c_blue_end.dec.degree - c_blue_start.dec.degree  # type: ignore
        dx_red = c_red_end.ra.degree - c_red_start.ra.degree  # type: ignore
        dy_red = c_red_end.dec.degree - c_red_start.dec.degree  # type: ignore
        # add blue and redshifted arrow
        plt.arrow(
            c_blue_start.ra.degree,  # type: ignore
            c_blue_start.dec.degree,  # type: ignore
            dx_blue,
            dy_blue,
            lw=1,
            fc="dodgerblue",
            ec="k",
            width=default_width * arrow_width,
            head_width=default_head_width * arrow_width,
            alpha=0.7,
            transform=ax.get_transform("fk5"),
            zorder=20,
        )
        plt.arrow(
            c_red_start.ra.degree,  # type: ignore
            c_red_start.dec.degree,  # type: ignore
            dx_red,
            dy_red,
            lw=1,
            fc="crimson",
            ec="k",
            width=default_width * arrow_width,
            head_width=default_head_width * arrow_width,
            alpha=0.7,
            transform=ax.get_transform("fk5"),
            zorder=21,
        )


def annotate_panel(
    ax: Axes,
    wcs: WCS,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
) -> None:
    """
    Add source names/markers and outflow orientations to a panel, using the
    styling shared by all panel plots.
    """
    ax.autoscale(enable=False)
    if do_annotation:
        annotate_sources(
            ax,
            wcs,
            color="white",
            color_back="black",
            fontsize=10,
            marker=do_marker,
            label=True,
            label_offset=4.0 * u.arcsec,  # type: ignore
            connect_line=True,
        )
    if do_outflow:
        annotate_outflow(ax, wcs, arrow_width=2.0)


def validate_frequency(frequency: u.Quantity[u.Hz]) -> bool:
    """
    Function to validate the frequency.
    Parameters:
    frequency: frequency in units of Hz (e.g., 1*u.GHz)

    Returns:
    True if the frequency is valid.
    """
    frequency.to(u.Hz)  # type: ignore
    return True


def pb_telecope(frequency: u.Quantity[u.Hz], telescope: str = "NOEMA") -> u.Quantity:
    """
    Function to compute the primary beam of the NOEMA telescope.
    Parameters:
    frequency: frequency in Hz
    telescope: name of the telescope

    Returns:
    primary beam in degrees
    """
    validate_frequency(frequency)
    if telescope == "NOEMA":
        # NOEMA primary beam
        pb = (64.1 * u.arcsec * 72.78382 * u.GHz / frequency).decompose()  # type: ignore
    elif telescope == "ALMA":
        pb = (19.0 * u.arcsec * 300 * u.GHz / frequency).decompose()  # type: ignore
    elif telescope == "SMA":
        pb = (36.0 * u.arcsec * 345 * u.GHz / frequency).decompose()  # type: ignore
    elif telescope == "VLA":
        pb = (45.0 * u.arcmin * 1 * u.GHz / frequency).decompose()  # type: ignore
    elif telescope == "30m":
        pb = (2460 * u.arcsec * 1 * u.Ghz / frequency).decompose()  # type: ignore
    else:
        raise ValueError(
            "Telescope not supported. Please choose NOEMA, ALMA, SMA, VLA, or 30m."
        )
    return pb.to(u.degree)  # type: ignore


def plot_PB(
    ax: Axes,
    header: Header,
    ra0: u.Quantity,
    dec0: u.Quantity,
    color: str = "white",
    lw: float = 1.0,
) -> None:
    frequency = get_frequency(header) * u.GHz
    pb_noema = pb_telecope(frequency, telescope="NOEMA")
    circ = SphericalCircle(
        (ra0, dec0),
        pb_noema / 2.0,
        ls=(0, (5, 10)),
        lw=lw,
        edgecolor=color,
        facecolor="none",
        transform=ax.get_transform("fk5"),
    )
    ax.add_patch(circ)


def prepare_color_map(cmap: Colormap | str, color_nan: str = "0.1") -> Colormap:
    """Build a colormap instance with NaN pixels colored, shared by all panel plots."""
    base_cmap = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    return base_cmap.with_extremes(bad=color_nan)


def add_scalebar_and_beam(
    ax: Axes,
    header: Header,
    label_col: str = "black",
    bkgrd_col: str = "white",
    show_beam: bool = True,
    with_stroke: bool = False,
) -> None:
    """Add the standard 1000 au scale bar and synthesized beam to a panel."""
    length = (1e3 * u.au / (distance * u.pc)).to(u.deg, u.dimensionless_angles())
    add_scalebar(ax, length, label=r"1\,000 au", color=label_col, corner="bottom right")
    if with_stroke:
        scalebar = ax.artists[-1]  # get the last added artist, which is the scalebar
        scalebar.txt_label._text.set_path_effects(
            [PathEffects.withStroke(linewidth=1.0, foreground=bkgrd_col)]
        )
        scalebar.size_bar.get_children()[0].set_path_effects(
            [PathEffects.withStroke(linewidth=2.0, foreground=bkgrd_col)]
        )
    if show_beam:
        add_beam(
            ax,
            header=header,
            frame=False,
            pad=0.2,
            color=label_col,
            corner="bottom left",
        )


def add_side_colorbar(
    fig: Figure,
    ax: Axes,
    im: AxesImage,
    label: str | None = None,
    label_fontsize: float | None = None,
    nbins: int = 5,
    fmt: str = "{x:.0f}",
    pad: float = 0.005,
    width: float = 0.025,
) -> None:
    """Add a colorbar to the right of ax, styled consistently across the panel plots."""
    cax = fig.add_axes(
        [
            ax.get_position().x1 + pad,
            ax.get_position().y0,
            width,
            ax.get_position().height,
        ]
    )
    cb = fig.colorbar(im, cax=cax)
    if label is not None:
        cb.set_label(label, fontsize=label_fontsize)
    cb.ax.yaxis.set_tick_params(color="black", labelcolor="black", direction="out")
    cb.ax.locator_params(nbins=nbins)
    cb.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter(fmt))


def plot_continuum_panel(
    data_cont: NDArray[np.float64],
    header: Header,
    wcs: WCS,
    ax: Axes,
    noise_cont: float,
    color_map: Colormap | str = "inferno",
    vmin: float | None = None,
    vmax: float | None = None,
    ra0: u.Quantity | None = None,
    dec0: u.Quantity | None = None,
    show_pb: bool = True,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    show_beam: bool = True,
    label_col: str = "black",
    bkgrd_col: str = "white",
    do_offsets: bool = False,
) -> AxesImage:
    """
    Draw a single continuum panel (image, contours, PB, annotations, scalebar, beam)
    into an existing axis. Does not create a figure, colorbar, or save to disk, so it
    can be reused in a loop to build multipanel figures.
    Parameters:
    data_cont: continuum data (e.g., from load_continuum_data)
    header: FITS header matching data_cont
    wcs: WCS matching data_cont, used for the transforms and PB/contour overlays
    ax: axis to draw into (must already be created with projection=wcs)
    color_map: colormap instance to use for the image
    noise_cont: estimated noise in data_cont, used for contour levels
    vmin, vmax: color scale limits. If None, defaults to -5*noise_cont and 0.3*max(data_cont)
    ra0, dec0: region center in degrees, required if show_pb is True or do_offsets is True
    show_pb: if True, the primary beam circle is drawn (requires ra0, dec0)
    do_marker: if True, markers are added to the source positions
    do_outflow: if True, outflow orientations are added to the plot
    do_annotation: if True, source names are added to the plot
    show_beam: if True, the synthesized beam is drawn
    label_col: color used for the scalebar and beam
    do_offsets: if True, the axes are displayed as offsets from (ra0, dec0)

    Returns:
    im: the AxesImage returned by imshow, for use with fig.colorbar()
    """
    vmin, vmax = default_minmax(data_cont, noise_cont, vmin, vmax, vmax_scale=0.3)

    im = ax.imshow(
        data_cont,
        origin="lower",
        interpolation="None",
        cmap=color_map,
        alpha=1.0,
        transform=ax.get_transform(wcs),
        vmin=vmin,
        vmax=vmax,
    )
    if show_pb:
        if ra0 is None or dec0 is None:
            raise ValueError("ra0 and dec0 are required when show_pb is True.")
        plot_PB(ax, header, ra0, dec0, color="white", lw=1.0)
        plot_PB(ax, header, ra0, dec0, color="black", lw=0.5)

    plot_data_contours(ax, data_cont, wcs, noise_cont)
    annotate_panel(
        ax, wcs, do_marker=do_marker, do_outflow=do_outflow, do_annotation=do_annotation
    )
    add_scalebar_and_beam(
        ax,
        header,
        label_col=label_col,
        bkgrd_col=bkgrd_col,
        show_beam=show_beam,
        with_stroke=True,
    )

    if do_offsets:
        if ra0 is None or dec0 is None:
            raise ValueError("ra0 and dec0 are required when do_offsets is True.")
        center_coord = SkyCoord(ra=ra0, dec=dec0, unit=(u.deg, u.deg))
    else:
        center_coord = None
    prodige_style(ax, do_offsets=do_offsets, center_coord=center_coord)
    return im


def hide_axis_labels(ax: Axes, hide_x: bool = True, hide_y: bool = True) -> None:
    # use the offset overlay coords (set by prodige_style with do_offsets=True) if present,
    # since those - not ax.coords - carry the visible tick/axis labels in that mode
    lon, lat = getattr(ax, "_offset_coords", (ax.coords[0], ax.coords[1]))  # type: ignore[reportUnknownMemberType]
    if hide_x:
        lon.set_ticklabel_visible(False)
        lon.set_axislabel("")
    if hide_y:
        lat.set_ticklabel_visible(False)
        lat.set_axislabel("")


def prepare_continuum_panel(
    region: str,
    bb: str,
    data_directory: str,
    mosaic: bool = False,
) -> tuple[NDArray[np.float64], Header, WCS, float, u.Quantity, u.Quantity]:
    """
    Load and preprocess the continuum data needed to draw one panel with
    plot_continuum_panel. The WCS is returned separately since it is needed
    to create the axis (projection=wcs) before the panel can be drawn.
    Parameters:
    region: name of the region
    bb: baseband of the data (lo, li, ui, uo)
    data_directory: directory where the data is stored
    mosaic: if True, mosaic data is used. This changes the filename format of the data.

    Returns:
    data_cont, header, wcs_cont, noise_cont, ra0, dec0
    """
    file_name = filename_continuum(region, bb, mosaic)
    data_cont, noise_cont, header = load_continuum_data(
        data_directory + file_name, region
    )
    wcs_cont = WCS(header)
    ra0, dec0 = get_region_center(region)
    return data_cont, header, wcs_cont, noise_cont, ra0, dec0


def plot_continuum_grid(
    panels: list[tuple[str, str]],
    data_directory: str,
    fig_directory: str = "./",
    fig_name: str = "continuum_grid.pdf",
    ncols: int = 3,
    panel_size: tuple[float, float] = (3.0, 2.7),
    cmap: Colormap | str = cmap_default,
    color_nan: str = "0.1",
    vmin: float | None = None,
    vmax: float | None = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    do_offsets: bool = False,
    labels: list[str] | None = None,
    label_col: str = "black",
    bkgrd_col: str = "white",
    show_colorbar: bool = True,
    save_fig: bool = True,
) -> tuple[Figure, list[Axes]]:
    """
    Plot a grid of continuum panels, one per (region, bb) pair in `panels`.
    Each panel keeps its own WCS projection, since different regions/basebands
    are not guaranteed to share a common pointing or pixel grid.
    Parameters:
    panels: list of (region, bb) tuples, one entry per panel
    data_directory: directory where the data is stored
    fig_directory: directory where the figure will be stored
    fig_name: filename of the saved figure
    ncols: number of columns in the grid
    panel_size: (width, height) in inches of a single panel
    cmap: colormap for the plot (default is the one listed in config.py)
    color_nan: color for NaN values
    vmin, vmax: shared color scale limits. If None, computed per panel as in plot_continuum
    mosaic: if True, mosaic data is used for all panels
    do_marker: if True, markers are added to the source positions
    do_outflow: if True, outflow orientations are added to the plot
    do_annotation: if True, source names are added to the plot
    do_offsets: if True, the axes are displayed as offsets from the region center
    labels: optional list of text labels, one per panel
    label_col: color used for the scale bar, beam, and panel labels
    show_colorbar: if True, a colorbar is added to each panel
    save_fig: if True, the figure is saved to disk

    Returns:
    fig, axs: the created figure and the list of axes (one per panel)
    """
    plt.rcParams.update(pyplot_params)
    color_map = prepare_color_map(cmap, color_nan)

    n_panels = len(panels)
    nrows = int(np.ceil(n_panels / ncols))
    panel_width, panel_height = panel_size
    fig = plt.figure(figsize=(panel_width * ncols + 1, panel_height * nrows + 1))  # type: ignore[reportUnknownMemberType]
    gs = fig.add_gridspec(nrows, ncols, hspace=0.0, wspace=0.0)  # type: ignore[reportUnknownMemberType]

    axs: list[Axes] = []
    for index, (region, bb) in enumerate(panels):
        row, col = divmod(index, ncols)
        data_cont, header, wcs_cont, noise_cont, ra0, dec0 = prepare_continuum_panel(
            region, bb, data_directory, mosaic
        )

        ax: Axes = fig.add_subplot(gs[row, col], projection=wcs_cont)
        axs.append(ax)
        im = plot_continuum_panel(
            data_cont,
            header,
            wcs_cont,
            ax,
            noise_cont,
            color_map=color_map,
            vmin=vmin,
            vmax=vmax,
            ra0=ra0,
            dec0=dec0,
            show_pb=not mosaic,
            do_marker=do_marker,
            do_outflow=do_outflow,
            do_annotation=do_annotation,
            show_beam=True,
            label_col=label_col,
            bkgrd_col=bkgrd_col,
            do_offsets=do_offsets,
        )
        hide_axis_labels(ax, hide_x=(row == 0) | (col != 0), hide_y=(col != 0))

        if labels is not None:
            label_text = ax.text(
                0.05,
                0.95,
                labels[index],
                transform=ax.transAxes,
                color=label_col,
                verticalalignment="top",
            )  # type: ignore[reportUnknownMemberType]
            label_text.set_path_effects(
                [PathEffects.withStroke(linewidth=1.0, foreground=label_col_back)]
            )  # type: ignore[reportUnknownMemberType]

        if show_colorbar:
            wavelength = get_wavelength(header)
            add_side_colorbar(
                fig,
                ax,
                im,
                label=r"$I_{" + str(wavelength) + "\\, \\rm mm}$ (mJy\\,beam$^{-1}$)",
                label_fontsize=8,
                nbins=4,
                width=0.015,
            )
    fig.tight_layout()
    if save_fig:
        fig.savefig(
            fig_directory + fig_name,
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )  # type: ignore[reportUnknownMemberType]
        plt.close(fig)

    return fig, axs


def plot_continuum(
    region: str,
    bb: str,
    data_directory: str,
    fig_directory: str = "./",
    cmap: Colormap | str = cmap_default,
    color_nan: str = "0.1",
    vmin: float | None = None,
    vmax: float | None = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    do_offsets: bool = False,
    label_col: str = "black",
    bkgrd_col: str = "white",
    save_fig: bool = True,
) -> None:
    """
    Function to plot the continuum data with the sources and outflow orientations.
    Labels and annotations are added to the plot.
    Parameters:
    region: name of the region
    bb: baseband of the data (lo, li, ui, uo)
    data_directory: directory where the data is stored
    fig_directory: directory where the figure will be stored
    cmap: colormap for the plot (default is the one listed in config.py)
    color_nan: color for NaN values
    vmin: minimum value for the color scale. If None, it is set to -5*noise
    vmax: maximum value for the color scale. If None, it is set to 0.3*max(data)
    mosaic: if True, mosaic data is used. This changes the filename format of the data.
    do_marker: if True, markers are added to the source positions
    do_outflow: if True, outflow orientations are added to the plot
    do_annotation: if True, source names are added to the plot
    do_offsets: if True, the axes are displayed as offsets from the region center

    Returns:
    A PDF file with the continuum plot is saved on disk with the following name
             'continuum_' + region + '_' + bb + '.pdf'
    """
    # plot continuum in color and contours, add source names, add outflow directions
    # use general plot parameters
    plt.rcParams.update(pyplot_params)  # type: ignore[reportUnknownMemberType]
    color_map = prepare_color_map(cmap, color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    ra0, dec0 = get_region_center(region)
    # load continuum data
    file_name = filename_continuum(region, bb, mosaic)
    data_cont, noise_cont, hd_cont = load_continuum_data(
        data_directory + file_name, region
    )
    vmin, vmax = default_minmax(data_cont, noise_cont, vmin, vmax, vmax_scale=0.3)

    wavelength = get_wavelength(hd_cont)
    wcs_cont = WCS(hd_cont)

    # create figure
    fig = plt.figure(1, figsize=(fig_width, fig_height))  # type: ignore[reportUnknownMemberType]
    ax: Axes = plt.subplot(1, 1, 1, projection=wcs_cont)  # type: ignore[reportUnknownMemberType]

    im = plot_continuum_panel(
        data_cont,
        hd_cont,
        wcs_cont,
        ax,
        noise_cont,
        color_map=color_map,
        vmin=vmin,
        vmax=vmax,
        ra0=ra0,
        dec0=dec0,
        show_pb=not mosaic,
        do_marker=do_marker,
        do_outflow=do_outflow,
        do_annotation=do_annotation,
        show_beam=True,
        label_col=label_col,
        bkgrd_col=bkgrd_col,
        do_offsets=do_offsets,
    )
    # Get coordinates for colorbar
    add_side_colorbar(
        fig,
        ax,
        im,
        label=r"$I_{" + str(wavelength) + "\\, \\rm mm}$ (mJy\\,beam$^{-1}$)",
    )
    # save plot
    if save_fig:
        fig.savefig(
            fig_directory + "continuum_" + region + "_" + bb + ".pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close(fig)


def plot_line_mom0(
    region: str,
    linename: str,
    bb: str,
    data_directory: str,
    fig_directory: str = "./",
    cmap: Colormap | str = cmap_mom0_default,
    color_nan: str = "0.1",
    vmin: float | None = None,
    vmax: float | None = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    save_fig: bool = True,
    label_col_TdV: str = "white",
    bkgrd_col_TdV: str = "black",
) -> None:
    # use general plot parameters
    plt.rcParams.update(pyplot_params)  # type: ignore[reportUnknownMemberType]
    color_map = prepare_color_map(cmap, color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    ra0, dec0 = get_region_center(region)
    # load integrated intensity data
    file_name = filename_line_TdV(region, linename, mosaic)
    data, noise_map, hd_TdV = load_line_TdV(data_directory + file_name, region)
    vmin, vmax = default_minmax(data, noise_map, vmin, vmax, vmax_scale=1.0)

    wcs_TdV = WCS(hd_TdV)

    # create figure
    fig: Figure = plt.figure(1, figsize=(fig_width, fig_height))  # type: ignore[reportUnknownMemberType]
    ax: Axes = plt.subplot(1, 1, 1, projection=wcs_TdV)  # type: ignore[reportUnknownMemberType]
    # plot continuum in color
    ax.imshow(
        data,
        origin="lower",
        interpolation="None",
        cmap=color_map,
        alpha=1.0,
        vmin=vmin,
        vmax=vmax,
    )
    if mosaic == False:
        plot_PB(ax, hd_TdV, ra0, dec0)
    plot_data_contours(ax, data, wcs_TdV, noise_map)
    annotate_panel(
        ax,
        wcs_TdV,
        do_marker=do_marker,
        do_outflow=do_outflow,
        do_annotation=do_annotation,
    )
    prodige_style(ax)

    add_scalebar_and_beam(
        ax,
        hd_TdV,
        label_col=label_col_TdV,
        bkgrd_col=bkgrd_col_TdV,
        show_beam=True,
        with_stroke=True,
    )
    # save plot
    if save_fig:
        fig.savefig(
            fig_directory + region + "_" + linename + "_TdV.pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close(fig)


def plot_line_vlsr(
    region: str,
    linename: str,
    data_directory: str,
    fig_directory: str = "./",
    cmap: Colormap | str = cmap_vlsr_default,
    color_nan: str = "0.1",
    vmin: float | None = None,
    vmax: float | None = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    do_offsets: bool = False,
    save_fig: bool = True,
    label_col_Vlsr: str = "black",
    bkgrd_col_Vlsr: str = "white",
) -> None:
    """
    Function to plot the line centroid velocity data with the sources and outflow orientations.
    Labels and annotations are added to the plot.
    Parameters:
    region: name of the region
    linename: line name (e.g., 'H2CO_l21')
    data_directory: directory where the data is stored
    fig_directory: directory where the figure will be stored
    cmap: colormap for the plot (default is the one listed in config.py)
    color_nan: color for NaN values
    vmin: minimum value for the color scale. If None, it is set to minimum value of the data
    vmax: maximum value for the color scale. If None, it is set to maximum value of the data
    if vmin and vmax are not set, the color scale is symmetric around the line center, with a width estimated from the largest from minimum and maximum separation between 5 and 95 percentail and Vlsr value from the source catalogue.
    mosaic: if True, mosaic data is used. This changes the filename format of the data.
    do_marker: if True, markers are added to the source positions
    do_outflow: if True, outflow orientations are added to the plot
    do_annotation: if True, source names are added to the plot
    do_offsets: if True, the axes are displayed in offsets
    save_fig: if True, the figure is saved to disk
    """
    # use general plot parameters
    plt.rcParams.update(pyplot_params)  # type: ignore[reportUnknownMemberType]
    color_map = prepare_color_map(cmap, color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    ra0, dec0 = get_region_center(region)
    v_lsr = get_region_vlsr(region)
    # load integrated intensity data
    file_name = filename_line_vlsr(region, linename, mosaic)
    file_TdV = filename_line_TdV(region, linename, mosaic)

    # load velocity data
    hdu = load_cutout(data_directory + file_name, source=region, is_hdu=False)
    data: NDArray[np.float64] = hdu.data
    # load integrated intensity data
    data_TdV, noise_map, hd_TdV = load_line_TdV(data_directory + file_TdV, region)
    if vmin is None and vmax is None:
        vmin_tmp, vmax_tmp = np.nanpercentile(data, [5, 95])
        delta = np.max([np.abs(vmin_tmp - v_lsr), np.abs(vmax_tmp - v_lsr)])
        vmin = v_lsr - delta
        vmax = v_lsr + delta
    elif vmax is None:
        vmax = np.nanmax(data)
    elif vmin is None:
        vmin = np.nanmin(data)

    wcs_TdV = WCS(hd_TdV)
    wcs_Vlsr = WCS(hdu.header)

    # create figure
    fig: Figure = plt.figure(1, figsize=(fig_width, fig_height))
    ax: Axes = plt.subplot(1, 1, 1, projection=wcs_Vlsr)
    # plot continuum in color
    im = ax.imshow(
        data,
        origin="lower",
        interpolation="None",
        cmap=color_map,
        alpha=1.0,
        vmin=vmin,
        vmax=vmax,
    )
    if mosaic == False:
        plot_PB(ax, hd_TdV, ra0, dec0, color=label_col_Vlsr)
    plot_data_contours(ax, data_TdV, wcs_TdV, noise_map)

    # annotate source names
    annotate_panel(
        ax,
        wcs_TdV,
        do_marker=do_marker,
        do_outflow=do_outflow,
        do_annotation=do_annotation,
    )
    # style
    prodige_style(
        ax,
        do_offsets=do_offsets,
        center_coord=SkyCoord(ra=ra0, dec=dec0, unit=(u.deg, u.deg)),  # type: ignore
    )

    add_side_colorbar(fig, ax, im, nbins=5, fmt="{x:.1f}")
    add_scalebar_and_beam(
        ax,
        hd_TdV,
        label_col=label_col_Vlsr,
        bkgrd_col=bkgrd_col_Vlsr,
        show_beam=True,
        with_stroke=True,
    )
    # save plot
    if save_fig:
        fig.savefig(
            fig_directory + region + "_" + linename + "_Vlsr.pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close(fig)
