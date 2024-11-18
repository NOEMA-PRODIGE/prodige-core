from __future__ import annotations
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import SphericalCircle, add_beam, add_scalebar
from astropy.stats import sigma_clipped_stats

import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.patheffects as PathEffects

from .source_catalogue import (
    load_sources_table,
    load_cutout,
    get_figsize,
    get_region_center,
    get_outflow_information,
    get_region_vlsr
)

from .config import pyplot_params, distance, cmap_default, cmap_mom0_default, cmap_vlsr_default

# name of the region
label_col = "black"
label_col_back = "white"


def determine_noise_map(data_2d: np.ndarray) -> float:
    """
    Determine the noise in the continuum data.
    """
    # compute noise in continuum data
    if isinstance(data_2d, np.ndarray) == False:
        raise ValueError("Input data is not a numpy array.")
    noise_2dmap = sigma_clipped_stats(data_2d, sigma=3.0)[-1]
    return noise_2dmap


def get_contour_params(maximum: float, noise: float) -> tuple[float, float]:
    """
    Compute the contour levels for the continuum data.
    maximum: maximum value to be shown in the plot
    noise: noise in the data
    """
    # compute contour levels at -5,5,10,20,40,80x... sigma
    # determines the number of contours to be plotted
    steps = int(np.log(maximum / (5.0 * noise)) // np.log(2.0)) + 1
    if steps < 1:
        return [0], ['solid'], False
    steps_arr = np.logspace(
        start=0,
        stop=steps,
        num=steps,
        endpoint=False,
        base=2.0,
        dtype=None,
        axis=0,
    )
    # append -5 sigma to the array and multiply by step size
    steps_arr = np.append(-steps_arr[0], steps_arr) * 5.0 * noise
    line_styles = ["dotted"] + ["solid"] * steps
    return steps_arr, line_styles, True


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
) -> tuple[np.ndarray, float, fits.header.Header]:
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
    hdu_cont.data[hdu_cont.data == 0.0] = np.nan
    # Update the header with the updated WCS from the cutout, as well as the data in mJy/beam.
    header = hdu_cont.header
    if header["BUNIT"].casefold() == "JY/BEAM".casefold():
        data_cont = np.squeeze(hdu_cont.data) * 1e3
        header["BUNIT"] = "mJy/beam"
    else:
        data_cont = np.squeeze(hdu_cont.data)
    # compute noise
    noise_cont = determine_noise_map(data_cont)
    return data_cont, noise_cont, header


def load_line_TdV(
    datafile: str,
    region: str,
) -> tuple[np.ndarray, float, fits.header.Header]:
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
    hdu.data[hdu.data == 0.0] = np.nan
    # Update the header with the updated WCS from the cutout, as well as the data in mJy/beam.
    header = hdu.header
    if header["BUNIT"].casefold() == "JY/BEAM KM/S".casefold():
        data = np.squeeze(hdu.data) * 1e3
        header["BUNIT"] = "mJy/beam km/s"
    elif header["BUNIT"].casefold() == "mJY/BEAM KM/S".casefold():
        data = np.squeeze(hdu.data)
        header["BUNIT"] = "mJy/beam km/s"
    else:
        header["BUNIT"] = "K km/s"
        data = np.squeeze(hdu.data)
    # compute noise
    noise_map = determine_noise_map(data)
    return data, noise_map, header


def get_frequency(header: fits.header.Header) -> float:
    """
    Function to get the frequency from the header and convert it to GHz.
    Parameters:
    header: FITS header

    Returns:
    frequency: frequency in GHz
    """
    if "RESTFREQ" not in header:
        raise ValueError("RESTFREQ not found in header.")
    else:
        restfreq = header["RESTFREQ"]  # Hz
    return restfreq * 1e-9


def get_wavelength(header: fits.header.Header) -> float:
    """
    Function to get the wavelength from the header.
    Parameters:
    header: FITS header

    Returns:
    wavelength: wavelength in mm
    """
    restfreq = get_frequency(header)
    wavelength = (restfreq * u.GHz).to(u.mm, equivalencies=u.spectral())
    return np.around(wavelength, decimals=1)  # mm


def prodige_style(ax: plt.Axes, do_offsets: bool = False, center_coord=None) -> None:
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
        DEC.set_ticklabel(rotation=90.0, color="black",
                          exclude_overlapping=True)
        RA.set_ticklabel(color="black", exclude_overlapping=True)
        DEC.set_ticks(spacing=10 * u.arcsec, color="black")
        RA.set_ticks(spacing=1.0 * 15 * u.arcsec, color="black")
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
        ra_offset.coord_wrap = 180*u.deg  # avoid wrapping
        ra_offset.display_minor_ticks(True)
        dec_offset.display_minor_ticks(True)
        dec_offset.set_minor_frequency(5)
        ra_offset.set_minor_frequency(5)
        dec_offset.set_ticks(spacing=15 * u.arcsec, color="black")
        ra_offset.set_ticks(spacing=15 * u.arcsec, color="black")


def annotate_sources(
    ax: plt.Axes,
    wcs: WCS.wcs,
    color: str = "cornflowerblue",
    color_back: str = "black",
    marker: bool = False,
    label: bool = True,
    connect_line: bool = False,
    fontsize: int = 10,
    label_offset: u.Quantity = 1.0 * u.arcsec,
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
        c = SkyCoord(ra=RA_i, dec=Dec_i, unit=(u.hourangle, u.deg))
        # Check if source is within the field of view
        if wcs.footprint_contains(c) == False:
            continue
        if marker == True:
            ax.scatter(
                c.ra,
                c.dec,
                marker="*",
                c=color,
                edgecolor="black",
                linewidth=0.5,
                s=20,
                transform=ax.get_transform("world"),
                zorder=40,
            )

        if label == True:
            c_label = c.directional_offset_by(
                offset_PA_i * u.deg, label_offset)
            label_text = ax.text(
                c_label.ra.degree,
                c_label.dec.degree,
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

        if connect_line == True:
            c_line_start = c.directional_offset_by(
                offset_PA_i * u.deg, 0.2 * label_offset
            )
            c_line_end = c.directional_offset_by(
                offset_PA_i * u.deg, 0.5 * label_offset
            )
            ax.plot(
                [c_line_start.ra.degree, c_line_end.ra.degree],
                [c_line_start.dec.degree, c_line_end.dec.degree],
                color=color_back,
                lw=1.5,
                alpha=0.7,
                zorder=10,
                transform=ax.get_transform("world"),
            )
            ax.plot(
                [c_line_start.ra.degree, c_line_end.ra.degree],
                [c_line_start.dec.degree, c_line_end.dec.degree],
                color=color,
                lw=1.0,
                alpha=0.7,
                zorder=10,
                transform=ax.get_transform("world"),
            )


def annotate_outflow(
    ax: plt.Axes,
    wcs: WCS.wcs,
    arrow_width: float = 1.0,
    arrow_length: u.Quantity = 3 * u.arcsec,
    arrow_offset: u.Quantity = 0.05 * u.arcsec,
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
        c = SkyCoord(ra=RA_i, dec=Dec_i, unit=(u.hourangle, u.deg))
        #  sources_Dec[k], unit=(u.hourangle, u.deg))
        # check if source is within the field of view
        # and if the outflow orientation is defined
        if (wcs.footprint_contains(c) & np.isfinite(source_outflowPA_i)) == False:
            continue
        c_blue_start = c.directional_offset_by(
            source_outflowPA_i * u.deg, arrow_offset)
        c_blue_end = c.directional_offset_by(
            source_outflowPA_i * u.deg, arrow_length)
        c_red_start = c.directional_offset_by(
            (180 + source_outflowPA_i) * u.deg, arrow_offset
        )
        c_red_end = c.directional_offset_by(
            (180 + source_outflowPA_i) * u.deg, arrow_length
        )
        # calculate the offset for the arrows
        dx_blue = c_blue_end.ra.degree - c_blue_start.ra.degree
        dy_blue = c_blue_end.dec.degree - c_blue_start.dec.degree
        dx_red = c_red_end.ra.degree - c_red_start.ra.degree
        dy_red = c_red_end.dec.degree - c_red_start.dec.degree
        # add blue and redshifted arrow
        plt.arrow(
            c_blue_start.ra.degree,
            c_blue_start.dec.degree,
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
            c_red_start.ra.degree,
            c_red_start.dec.degree,
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


def validate_frequency(frequency: u.Hz) -> bool:
    """
    Function to validate the frequency.
    Parameters:
    frequency: frequency in units of Hz (e.g., 1*u.GHz)

    Returns:
    True if the frequency is valid.
    """
    frequency.to(u.Hz)
    return True


def pb_telecope(frequency: u.Hz, telescope: str = "NOEMA") -> u.degree:
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
        pb = (64.1 * u.arcsec * 72.78382 * u.GHz / frequency).decompose()
    elif telescope == "ALMA":
        pb = (19.0 * u.arcsec * 300 * u.GHz / frequency).decompose()
    elif telescope == "SMA":
        pb = (36.0 * u.arcsec * 345 * u.GHz / frequency).decompose()
    elif telescope == "VLA":
        pb = (45.0 * u.arcmin * 1 * u.GHz / frequency).decompose()
    else:
        raise ValueError(
            "Telescope not supported. Please choose NOEMA, ALMA, SMA, or VLA."
        )
    return pb.to(u.degree)


def plot_PB(ax: plt.Axes, header: fits.header.Header, ra0: float, dec0: float, color: str = 'white') -> None:
    frequeny = get_frequency(header) * u.GHz
    pb_noema = pb_telecope(frequeny, telescope="NOEMA")
    circ = SphericalCircle(
        (ra0 * u.deg, dec0 * u.deg),
        pb_noema / 2.0,
        ls=(0, (5, 10)),
        lw=1.0,
        edgecolor=color,
        facecolor="none",
        transform=ax.get_transform("fk5"),
    )
    ax.add_patch(circ)


def plot_continuum(
    region: str,
    bb: str,
    data_directory: str,
    fig_directory: str = "./",
    cmap: str | None = None,
    color_nan: str = "0.1",
    vmin: float = None,
    vmax: float = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
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

    Returns:
    A PDF file with the continuum plot is saved on disk with the following name
             'continuum_' + region + '_' + bb + '.pdf'
    """
    # plot continuum in color and contours, add source names, add outflow directions

    if cmap == None:
        cmap = cmap_default
    # use general plot parameters
    plt.rcParams.update(pyplot_params)
    color_map = plt.get_cmap(cmap).copy()
    color_map.set_bad(color=color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    ra0, dec0 = get_region_center(region)
    # load continuum data
    file_name = filename_continuum(region, bb, mosaic)
    data_cont, noise_cont, hd_cont = load_continuum_data(
        data_directory + file_name, region
    )
    if vmin == None:
        vmin = -5.0 * noise_cont
    if vmax == None:
        vmax = 0.3 * np.nanmax(data_cont)

    wavelength = get_wavelength(hd_cont)
    wcs_cont = WCS(hd_cont)

    # create figure
    fig = plt.figure(1, figsize=(fig_width, fig_height))
    ax = plt.subplot(1, 1, 1, projection=wcs_cont)

    # plot continuum in color
    im = ax.imshow(
        data_cont,
        origin="lower",
        interpolation="None",
        cmap=color_map,
        alpha=1.0,
        transform=ax.get_transform(wcs_cont),
        vmin=vmin,
        vmax=vmax,
    )
    if mosaic == False:
        plot_PB(ax, hd_cont, ra0, dec0)

    # add continuum contour levels
    cont_levels, style_levels, valid_contour = get_contour_params(
        np.nanmax(data_cont), noise_cont)

    if valid_contour:
        ax.contour(
            data_cont,
            colors="white",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=0.75,
            transform=ax.get_transform(wcs_cont),
        )

        ax.contour(
            data_cont,
            colors="black",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=0.35,
            transform=ax.get_transform(wcs_cont),
        )

    # annotate source names
    ax.autoscale(enable=False)
    if do_annotation == True:
        annotate_sources(
            ax,
            wcs_cont,
            color="white",
            color_back="black",
            fontsize=10,
            marker=do_marker,
            label=True,
            label_offset=4.0 * u.arcsec,
            connect_line=True,
        )

    # add outflow orientations
    if do_outflow:
        annotate_outflow(ax, wcs_cont, arrow_width=2.0)
    prodige_style(ax)
    # Get coordinates for colorbar
    cax = fig.add_axes(
        [
            ax.get_position().x1 + 0.005,
            ax.get_position().y0,
            0.025,
            ax.get_position().height,
        ]
    )
    # add colorbar
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$I_{" + str(wavelength.value) +
                 "\\, \\rm mm}$ (mJy\\,beam$^{-1}$)")
    cb.ax.yaxis.set_tick_params(
        color="black", labelcolor="black", direction="out")
    cb.ax.locator_params(nbins=5)

    # cb.locator = MultipleLocator(10.0)
    cb.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    # add linear scale bar (1000 au)
    length = (1e3 * u.au / (distance * u.pc)
              ).to(u.deg, u.dimensionless_angles())
    add_scalebar(ax, length, label=r"1\,000 au",
                 color=label_col, corner="bottom right")
    # add beam
    add_beam(
        ax, header=hd_cont, frame=False, pad=0.2, color=label_col, corner="bottom left"
    )
    # save plot
    if save_fig:
        plt.savefig(
            fig_directory + "continuum_" + region + "_" + bb + ".pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close()


def plot_line_mom0(
    region: str,
    linename: str,
    bb: str,
    data_directory: str,
    fig_directory: str = "./",
    cmap: str | None = None,
    color_nan: str = "0.1",
    vmin: float = None,
    vmax: float = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    save_fig: bool = True,
) -> None:
    label_col_TdV = 'white'
    if cmap == None:
        cmap = cmap_mom0_default
    # use general plot parameters
    plt.rcParams.update(pyplot_params)
    color_map = plt.get_cmap(cmap).copy()
    color_map.set_bad(color=color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    ra0, dec0 = get_region_center(region)
    # load integrated intensity data
    file_name = filename_line_TdV(region, linename, mosaic)
    data, noise_map, hd_TdV = load_line_TdV(
        data_directory + file_name, region
    )
    if vmin == None:
        vmin = -5.0 * noise_map
    if vmax == None:
        vmax = np.nanmax(data)

    wcs_TdV = WCS(hd_TdV)

    # create figure
    fig = plt.figure(1, figsize=(fig_width, fig_height))
    ax = plt.subplot(1, 1, 1, projection=wcs_TdV)
    # plot continuum in color
    im = ax.imshow(
        data,
        origin="lower",
        interpolation="None",
        cmap=color_map,
        alpha=1.0,
        transform=ax.get_transform(wcs_TdV),
        vmin=vmin,
        vmax=vmax,
    )
    if mosaic == False:
        plot_PB(ax, hd_TdV, ra0, dec0)
    cont_levels, style_levels, valid_contour = get_contour_params(
        np.nanmax(data), noise_map)

    if valid_contour:
        ax.contour(
            data,
            colors="white",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=0.75,
            transform=ax.get_transform(wcs_TdV),
        )

        ax.contour(
            data,
            colors="black",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=0.35,
            transform=ax.get_transform(wcs_TdV),
        )

    # annotate source names
    ax.autoscale(enable=False)
    if do_annotation == True:
        annotate_sources(
            ax,
            wcs_TdV,
            color="white",
            color_back="black",
            fontsize=10,
            marker=do_marker,
            label=True,
            label_offset=4.0 * u.arcsec,
            connect_line=True,
        )

    # add outflow orientations
    if do_outflow:
        annotate_outflow(ax, wcs_TdV, arrow_width=2.0)
    prodige_style(ax)

   # Get coordinates for colorbar
    # cax = fig.add_axes(
    #     [
    #         ax.get_position().x1 + 0.005,
    #         ax.get_position().y0,
    #         0.025,
    #         ax.get_position().height,
    #     ]
    # )
    # add colorbar
    # cb = fig.colorbar(im, cax=cax)
    # cb.set_label(r"$I_{" + str(wavelength.value) + "}$ mm (mJy\\,beam$^{-1}$)")
    # cb.ax.yaxis.set_tick_params(
    #     color="black", labelcolor="black", direction="out")
    # cb.ax.locator_params(nbins=5)

    # cb.locator = MultipleLocator(10.0)
    # cb.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))
    # add linear scale bar (1000 au)
    length = (1e3 * u.au / (distance * u.pc)
              ).to(u.deg, u.dimensionless_angles())
    add_scalebar(ax, length, label=r"1\,000 au",
                 color=label_col_TdV, corner="bottom right")
    # add beam
    add_beam(
        ax, header=hd_TdV, frame=False, pad=0.2, color=label_col_TdV, corner="bottom left"
    )
    # save plot
    if save_fig:
        plt.savefig(
            fig_directory + region + "_" + linename + "_TdV.pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close()


def plot_line_vlsr(
    region: str,
    linename: str,
    data_directory: str,
    fig_directory: str = "./",
    cmap: str | None = None,
    color_nan: str = "0.1",
    vmin: float = None,
    vmax: float = None,
    mosaic: bool = False,
    do_marker: bool = False,
    do_outflow: bool = False,
    do_annotation: bool = True,
    do_offsets: bool = False,
    save_fig: bool = True,
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
    label_col_Vlsr = 'black'
    if cmap == None:
        cmap = cmap_vlsr_default
    # use general plot parameters
    plt.rcParams.update(pyplot_params)
    color_map = plt.get_cmap(cmap).copy()
    color_map.set_bad(color=color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    ra0, dec0 = get_region_center(region)
    v_lsr = get_region_vlsr(region)
    # load integrated intensity data
    file_name = filename_line_vlsr(region, linename, mosaic)
    file_TdV = filename_line_TdV(region, linename, mosaic)

    # load velocity data
    hdu = load_cutout(data_directory + file_name, source=region, is_hdu=False)
    data = hdu.data
    # load integrated intensity data
    data_TdV, noise_map, hd_TdV = load_line_TdV(
        data_directory + file_TdV, region
    )
    if vmin == None and vmax == None:
        vmin, vmax = np.nanpercentile(data, [5, 95])
        delta = np.max([np.abs(vmin-v_lsr), np.abs(vmax-v_lsr)])
        vmin = v_lsr - delta
        vmax = v_lsr + delta
    elif vmax == None:
        vmax = np.nanmax(data)
    else:
        vmin = np.nanmin(data)

    wcs_TdV = WCS(hd_TdV)
    wcs_Vlsr = WCS(hdu.header)

    # create figure
    fig = plt.figure(1, figsize=(fig_width, fig_height))
    ax = plt.subplot(1, 1, 1, projection=wcs_Vlsr)
    # plot continuum in color
    im = ax.imshow(
        data,
        origin="lower",
        interpolation="None",
        cmap=color_map,
        alpha=1.0,
        transform=ax.get_transform(wcs_Vlsr),
        vmin=vmin,
        vmax=vmax,
    )
    if mosaic == False:
        plot_PB(ax, hd_TdV, ra0, dec0, color=label_col_Vlsr)
    #
    cont_levels, style_levels, valid_contour = get_contour_params(
        np.nanmax(data_TdV), noise_map)

    if valid_contour:
        ax.contour(
            data_TdV,
            colors="white",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=0.75,
            transform=ax.get_transform(wcs_TdV),
        )

        ax.contour(
            data_TdV,
            colors="black",
            alpha=1.0,
            levels=cont_levels,
            linestyles=style_levels,
            linewidths=0.35,
            transform=ax.get_transform(wcs_TdV),
        )

    # annotate source names
    ax.autoscale(enable=False)
    if do_annotation == True:
        annotate_sources(
            ax,
            wcs_TdV,
            color="white",
            color_back="black",
            fontsize=10,
            marker=do_marker,
            label=True,
            label_offset=4.0 * u.arcsec,
            connect_line=True,
        )

    # add outflow orientations
    if do_outflow:
        annotate_outflow(ax, wcs_TdV, arrow_width=2.0)
    # style
    prodige_style(ax, do_offsets=do_offsets, center_coord=SkyCoord(
        ra=ra0, dec=dec0, unit=(u.deg, u.deg)))

   # Get coordinates for colorbar
    cax = fig.add_axes(
        [
            ax.get_position().x1 + 0.005,
            ax.get_position().y0,
            0.025,
            ax.get_position().height,
        ]
    )
    # add colorbar
    cb = fig.colorbar(im, cax=cax)
    # cb.set_label(r"$V_{LSR}$ (km \\,s$^{-1}$)")
    cb.ax.yaxis.set_tick_params(
        color="black", labelcolor="black", direction="out")
    cb.ax.locator_params(nbins=5)

    # cb.locator = MultipleLocator(10.0)
    cb.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
    # add linear scale bar (1000 au)
    length = (1e3 * u.au / (distance * u.pc)
              ).to(u.deg, u.dimensionless_angles())
    add_scalebar(ax, length, label=r"1\,000 au",
                 color=label_col_Vlsr, corner="bottom right")
    # add beam
    add_beam(
        ax, header=hd_TdV, frame=False, pad=0.2, color=label_col_Vlsr, corner="bottom left"
    )
    # save plot
    if save_fig:
        fig.savefig(
            fig_directory + region + "_" + linename + "_Vlsr.pdf",
            format="pdf",
            bbox_inches="tight",
            pad_inches=0.01,
        )
        plt.close()
