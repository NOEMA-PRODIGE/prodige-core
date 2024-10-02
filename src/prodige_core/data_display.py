import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import add_beam, add_scalebar


import matplotlib.pyplot as plt
# from matplotlib.cm import cmap
from matplotlib.patches import Arrow
from matplotlib.ticker import MultipleLocator
from matplotlib import ticker
# from matplotlib.patches import Ellipse, Rectangle

from astropy.stats import sigma_clipped_stats

from .source_catalogue import load_sources_table, load_cutout, get_figsize
from .source_catalogue import region_dic, distance

from .config import pyplot_params, distance, cmap_default
# name of the region
# region = 'L1448N'
label_col = 'black'
label_col_back = 'white'

# distance to the region
# distance = 288.0  # pc


def determine_noise_map(data_2d):
    """
    Determine the noise in the continuum data.
    """
    # compute noise in continuum data

    # noise_continuum = np.nanstd(data_cont[170:210,100:130])
    noise_2dmap = sigma_clipped_stats(data_2d, sigma=3.0)[-1]
    return noise_2dmap


def get_contour_params(maximum, noise):
    """
    Compute the contour levels for the continuum data.
    maximum: maximum value to be shown in the plot
    noise: noise in the data
    """
    # compute contour levels at -5,5,10,20,40,80x... sigma
    # determines the number of contours to be plotted
    steps = int(np.log(maximum / (5.0 * noise)) / np.log(2.0)) + 1
    steps_arr = np.logspace(start=0, stop=steps+1, num=steps+1,
                            endpoint=False, base=2.0, dtype=None, axis=0)
    # append -5 sigma to the array and multiply by step size
    steps_arr = np.append(-steps_arr[0], steps_arr) * 5.0 * noise
    line_styles = ['dotted'] + ['solid'] * steps
    return steps_arr, line_styles


def load_continuum_data(data_directory, region, bb):
    """
    Function to load the continuum data and return the cutout specified in the dictionary.
    It return the data, estimated noise, and the FITS header.
    """
    # datafile
    datafile = region + '_CD_' + bb + '_cont_rob1-selfcal.fits'
    # header
    hdu_cont = load_cutout(data_directory + datafile,
                           source='L1448N', is_hdu=False)
    # world coordinate system
    header = hdu_cont.header  #
    # continuum data
    data_cont = np.squeeze(hdu_cont.data) * 1000.0  # mJy/beam
    # compute noise
    noise_cont = determine_noise_map(data_cont)
    return data_cont, noise_cont, header


def get_frequency(header):
    restfreq = header['RESTFREQ']  # Hz
    c_km_s = 299792.458  # km/s
    wavelength = (restfreq * u.Hz).to(u.mm, equivalencies=u.spectral())
    # wavelength =
    return np.around(wavelength, decimals=1)  # mm


def prodige_style(ax):
    # plot properties
    RA = ax.coords[0]
    DEC = ax.coords[1]
    RA.set_axislabel(r'$\alpha$ (J2000)', minpad=0.7)
    DEC.set_major_formatter('dd:mm:ss')
    RA.set_major_formatter('hh:mm:ss.s')
    DEC.set_axislabel(r'$\delta$ (J2000)', minpad=0.8)
    DEC.set_ticklabel(rotation=90., color='black', exclude_overlapping=True)
    RA.set_ticklabel(color='black', exclude_overlapping=True)
    DEC.set_ticks(spacing=10*u.arcsec, color='black')
    RA.set_ticks(spacing=1.0 * 15*u.arcsec, color='black')
    RA.display_minor_ticks(True)
    DEC.display_minor_ticks(True)
    DEC.set_minor_frequency(5)
    RA.set_minor_frequency(5)


def plot_continuum(region, bb, data_directory, cmap=None, color_nan='0.1', do_marker=False):
    # plot continuum in color and contours, add source names, add outflow directions

    if cmap == None:
        cmap = cmap_default
    # use general plot parameters
    plt.rcParams.update(pyplot_params)
    color_map = plt.get_cmap(cmap).copy()
    color_map.set_bad(color=color_nan)
    # figure size from dictionary
    fig_width, fig_height = get_figsize(region)
    # load continuum data
    data_cont, noise_cont, hd_cont = load_continuum_data(
        data_directory, region, bb)

    wavelength = get_frequency(hd_cont)
    wcs_cont = WCS(hd_cont)

    # create figure
    fig = plt.figure(1, figsize=(fig_width, fig_height))
    ax = plt.subplot(1, 1, 1, projection=wcs_cont)

    # plot continuum in color
    im = ax.imshow(data_cont, origin='lower', interpolation='None', cmap=color_map,
                   alpha=1.0, transform=ax.get_transform(
                       wcs_cont), vmin=-5.0*noise_cont, vmax=0.3*np.nanmax(data_cont))

    # add continuum contour levels
    cont_levels, style_levels = get_contour_params(
        np.nanmax(data_cont), noise_cont)

    ax.contour(data_cont, colors='white', alpha=1.0, levels=cont_levels,
               linestyles=style_levels, linewidths=1.0, transform=ax.get_transform(wcs_cont))

    ax.contour(data_cont, colors='black', alpha=1.0, levels=cont_levels,
               linestyles=style_levels, linewidths=0.5, transform=ax.get_transform(wcs_cont))

    # annotate source names
    ax.autoscale(enable=False)
    annotate_sources(ax, wcs_cont, color='white',
                     fontsize=10, marker=do_marker)

    # add outflow orientations
    # annotate_outflow(ax, wcs_cont, width=2.0)
    prodige_style(ax)
    # add colorbar
    cb = fig.colorbar(im, pad=0.0, shrink=0.855)
    cb.set_label(r'$I_{' + str(wavelength) +
                 '\mathrm{mm}}$ (mJy\,beam$^{-1}$)')
    cb.ax.yaxis.set_tick_params(
        color='black', labelcolor='black', direction='out')
    cb.locator = MultipleLocator(10.0)
    cb.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

    # add linear scale bar (1000 au)
    # Scalebar
    length = (1e3*u.au / (distance*u.pc)).to(u.deg, u.dimensionless_angles())
    add_scalebar(ax, length, label="1,000 au",
                 color=label_col, corner='bottom right')
    # add beam
    add_beam(ax, header=hd_cont, frame=False, pad=0.2,
             color=label_col, corner='top left')
    # save plot
    plt.savefig('continuum_' + region + '_' + bb + '.pdf',
                format='pdf', bbox_inches='tight', pad_inches=0.01)


def annotate_sources(ax, wcs, color='cornflowerblue', marker=False, label=True, fontsize=10):
    # annotate mm sources

    # load table containing sources within the region
    sources_name, sources_RA, sources_Dec, sources_color, sources_vlsr, sources_outflowPA = load_sources_table()

    # loop over all cores
    for k in range(sources_name.size):

        c = SkyCoord(sources_RA[k] + ' ' +
                     sources_Dec[k], unit=(u.hourangle, u.deg))

        sources_RA_pix, sources_Dec_pix = wcs.wcs_world2pix(c.ra, c.dec, 0)

        if (sources_name[k] == 'IRS3A'):
            off_x, off_y = -30, 0
        elif (sources_name[k] == 'IRS3B'):
            off_x, off_y = +40, -10
        elif (sources_name[k] == 'IRS3C'):
            off_x, off_y = +10, -10
        elif (sources_name[k] == 'mm'):
            off_x, off_y = -10, 0
        else:
            print('error')

        if marker == True:
            ax.scatter(c.ra, c.dec, marker='*', c=color, edgecolor='black', linewidth=0.5, s=20,
                       transform=ax.get_transform('world'))

        if label == True:
            ax.annotate(r'\textbf{'+str(sources_name[k])+r'}', xy=(sources_RA_pix, sources_Dec_pix), xytext=(sources_RA_pix+off_x, sources_Dec_pix+off_y), xycoords='data',
                        arrowprops=dict(color='white', arrowstyle='-', linestyle='-',
                                        linewidth=0.0, alpha=0.7, shrinkA=0, shrinkB=0),
                        color=color, fontsize=fontsize, ha='center', va='center', alpha=1.0, zorder=5)


def annotate_outflow(ax, wcs, width=5.0):
    # add outflow orientation angle

    # load table containing sources within the region
    sources_name, sources_RA, sources_Dec, sources_color, sources_vlsr, sources_outflowPA = load_sources_table()

    # loop over all cores
    for k in range(sources_name.size):

        # source coordinate
        c = SkyCoord(sources_RA[k] + ' ' +
                     sources_Dec[k], unit=(u.hourangle, u.deg))

        # converte coordinate to pixels
        sources_RA_pix, sources_Dec_pix = wcs.wcs_world2pix(c.ra, c.dec, 0)

        # length and start of arrow
        if sources_name[k] == 'IRS3A':
            arrow_length = 12.5  # pixel
            arrow_start = 0.4  # %
        elif sources_name[k] == 'IRS3B':
            arrow_length = 20.0  # pixel
            arrow_start = 0.4  # %
        elif sources_name[k] == 'IRS3C':
            arrow_length = 30.0  # pixel
            arrow_start = 0.1  # %
        elif sources_name[k] == 'mm':
            arrow_length = 410.0  # pixel
            arrow_start = 0.1  # %
        else:
            print('error')

        # minus because angle is defined counter clockwise
        dx = -1.0*np.sin(sources_outflowPA[k]
                         * np.pi/180.0)*arrow_length  # pixel
        dy = np.cos(sources_outflowPA[k]*np.pi/180.0)*arrow_length  # pixel

        # add blue and redshifted arrow
        plt.arrow(sources_RA_pix+arrow_start*dx, sources_Dec_pix+arrow_start*dy, dx, dy, lw=0.2,
                  fc='dodgerblue', ec='k', width=width, zorder=1, head_width=2.0*width, alpha=0.7)
        plt.arrow(sources_RA_pix-arrow_start*dx, sources_Dec_pix-arrow_start*dy, -1.0*dx, -1.0*dy,
                  lw=0.2, fc='crimson', ec='k', width=width, zorder=1, head_width=2.0*width, alpha=0.7)
