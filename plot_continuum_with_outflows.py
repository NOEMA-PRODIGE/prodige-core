import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS

import matplotlib.pyplot as plt
from matplotlib.patches import Arrow
from matplotlib.ticker import MultipleLocator
from matplotlib import ticker
from matplotlib.patches import Ellipse, Rectangle



#general plot parameters
params = { 'text.usetex': True,
    'font.size' : 12,
    'axes.labelsize' : 12,
    'axes.linewidth' : 1.5,
    'legend.fontsize': 12,
    'xtick.labelsize' : 12,
    'ytick.labelsize' : 12,
    'errorbar.capsize' : 1,
    'lines.linewidth'   : 1.0,
    'xtick.top' : True,
    'ytick.right' : True,
    'legend.fancybox' : False,
    'xtick.major.size' : 8.0 ,
    'xtick.minor.size' : 4.0,
    'ytick.major.size' : 8.0 ,
    'ytick.minor.size' : 4.0,
    'xtick.major.width' : 1.0,
    'xtick.minor.width' : 1.0,
    'ytick.major.width' : 1.0,
    'ytick.minor.width' : 1.0,
    'xtick.direction' : 'in',
    'ytick.direction' : 'in',
    'xtick.color' : 'black',
    'ytick.color' : 'black',
    'font.family' : 'serif'
}

def load_sources_table():
    #load table containing sources within FoV

    #sources within FOV
    sources_tab = np.loadtxt('sources.dat',dtype='U')

    #source name
    sources_name = sources_tab[:,0]
    #source RA
    sources_RA = sources_tab[:,1]
    #source DEC
    sources_Dec = sources_tab[:,2]
    #source color
    sources_color = sources_tab[:,3]
    #source vlsr
    sources_vlsr = sources_tab[:,4].astype(float)
    #source outflow PA
    sources_outflowPA = sources_tab[:,5].astype(float)

    return sources_name, sources_RA, sources_Dec, sources_color, sources_vlsr, sources_outflowPA


def compute_contour_levels(maximum,noise):
    ### compute contour levels at -5,5,10,20,40,80x... sigma

    #number of contour steps
    steps = int(np.log(maximum/(5.0*noise))/np.log(2.0))+1

    #first step: -5 sigma
    steps_arr_0 = np.array([-1.0])

    #other steps: 5,10,20,...sigma
    steps_arr_1 = np.logspace(start=0, stop=steps+1, num=steps+1, endpoint=False, base=2.0, dtype=None, axis=0)

    #concatenate both arrays
    steps_arr = np.append(steps_arr_0,steps_arr_1)

    #create contour level array
    cont_levels = 5.0*noise*steps_arr

    #return array containing contour levels
    return cont_levels


def load_continuum_data(data_directory, region, bb):

    #datafile
    datafile = region + '_CD_' + bb + '_cont_rob1-selfcal.fits'
    #header
    hdu_cont = fits.open(data_directory + datafile)[0]
    #world coordinate system
    wcs_cont = WCS(hdu_cont.header, naxis=['longitude', 'latitude'])
    #pixel scale
    delta = hdu_cont.header['CDELT2'] * 3600.0 #arcsec/pixel
    #synthesized beam
    bmin_cont = hdu_cont.header['BMIN'] * 3600.0 #arcsec
    bmaj_cont = hdu_cont.header['BMAJ']* 3600.0 #arcsec
    bpa_cont = hdu_cont.header['BPA'] #degree

    restfreq = hdu_cont.header['RESTFREQ'] #Hz
    c_km_s = 299792.458 #km/s
    wavelength = np.around(c_km_s / restfreq * 10.0**6.0,decimals=1) #mm

    #continuum data
    data_cont = hdu_cont.data[0,:,:]*1000.0 #mJy/beam

    #compute noise
    noise_cont = determine_noise_continuum(data_cont)

    #contour levels (-5,5,10,20,40,80,..sigma)
    cont_levels = compute_contour_levels(np.nanmax(data_cont),noise_cont)

    return data_cont, noise_cont, cont_levels, wcs_cont, delta, bmin_cont, bmaj_cont, bpa_cont, wavelength


def determine_noise_continuum(data_cont):
    ###compute noise in continuum data

    noise_continuum = np.nanstd(data_cont[170:210,100:130])

    return noise_continuum


def annotate_sources(ax,wcs,color='cornflowerblue',marker=False, label=True, fontsize=10):
    ### annotate mm sources

    #load table containing sources within the region
    sources_name, sources_RA, sources_Dec, sources_color, sources_vlsr, sources_outflowPA = load_sources_table()

    #loop over all cores
    for k in range(sources_name.size):

        c = SkyCoord(sources_RA[k] + ' ' + sources_Dec[k], unit=(u.hourangle, u.deg))

        sources_RA_pix, sources_Dec_pix =  wcs.wcs_world2pix(c.ra, c.dec, 0)

        if (sources_name[k] == 'IRS3A'):
            off_x, off_y = -30, 0
        elif (sources_name[k] == 'IRS3B'):
            off_x, off_y = +40, -10
        elif (sources_name[k] == 'IRS3C'):
            off_x, off_y = +10,-10
        elif (sources_name[k] == 'mm'):
            off_x, off_y = -10,0
        else:
            print('error')

        if marker==True:
            ax.plot(sources_RA_pix, sources_Dec_pix, '.', color=color, markersize=2)

        if label==True:
            ax.annotate(r'\textbf{'+str(sources_name[k])+r'}', xy=(sources_RA_pix, sources_Dec_pix),xytext=(sources_RA_pix+off_x, sources_Dec_pix+off_y), xycoords='data', arrowprops=dict(color='white', arrowstyle='-',linestyle='-',linewidth=0.0,alpha=0.7,shrinkA=0,shrinkB=0), color=color, fontsize=fontsize,ha='center',va='center',alpha=1.0, zorder=5)


def annotate_outflow(ax, wcs, width=5.0):
    #add outflow orientation angle

    #load table containing sources within the region
    sources_name, sources_RA, sources_Dec, sources_color, sources_vlsr, sources_outflowPA = load_sources_table()

    #loop over all cores
    for k in range(sources_name.size):

        #source coordinate
        c = SkyCoord(sources_RA[k] + ' ' + sources_Dec[k], unit=(u.hourangle, u.deg))

        #converte coordinate to pixels
        sources_RA_pix, sources_Dec_pix =  wcs.wcs_world2pix(c.ra, c.dec, 0)

        #length and start of arrow
        if sources_name[k] == 'IRS3A':
            arrow_length = 12.5 #pixel
            arrow_start = 0.4 #%
        elif sources_name[k] == 'IRS3B':
            arrow_length = 20.0 #pixel
            arrow_start = 0.4 #%
        elif sources_name[k] == 'IRS3C':
            arrow_length = 30.0 #pixel
            arrow_start = 0.1 #%
        elif sources_name[k] == 'mm':
            arrow_length = 410.0 #pixel
            arrow_start = 0.1 #%
        else:
            print('error')

        # minus because angle is defined counter clockwise
        dx = -1.0*np.sin(sources_outflowPA[k]*np.pi/180.0)*arrow_length #pixel
        dy = np.cos(sources_outflowPA[k]*np.pi/180.0)*arrow_length #pixel

        #add blue and redshifted arrow
        plt.arrow(sources_RA_pix+arrow_start*dx, sources_Dec_pix+arrow_start*dy, dx, dy, lw=0.2, fc='dodgerblue', ec='k', width=width,zorder=1, head_width=2.0*width,alpha=0.7)
        plt.arrow(sources_RA_pix-arrow_start*dx, sources_Dec_pix-arrow_start*dy, -1.0*dx, -1.0*dy, lw=0.2, fc='crimson', ec='k', width=width,zorder=1, head_width=2.0*width,alpha=0.7)


def add_continuum_contours(data_cont, cont_levels, color, wcs, ax):

    #positive contours: solid contours
    cont_ls =  ['solid' for x in range(cont_levels.size)]
    cont_lw = [0.3 for x in range(cont_levels.size)]
    #negative contour: dotted
    cont_ls[0] = 'dotted'
    cont_lw[0] = 0.5
    cont_lw[1] = 0.5

    #plot contours
    ax.contour(data_cont, colors=color, alpha=1.0, levels=cont_levels,linestyles=cont_ls,linewidths=cont_lw,transform=ax.get_transform(wcs))


def plot_continuum(region, bb, distance, data_directory):
### plot continuum in color and contours, add source names, add outflow directions

    #use general plot parameters
    plt.rcParams.update(params)

    #figure size
    fig_width = 3.4 #inch
    fig_height = 3.4 #inch

    #figure limits in pixel
    cut1_x = 60
    cut2_x = 260
    cut1_y = 92
    cut2_y = 292


    #load continuum data
    data_cont, noise_cont, cont_levels, wcs_cont, delta, bmin_cont, bmaj_cont, bpa_cont, wavelength = load_continuum_data(data_directory, region, bb)

    #create figure
    fig = plt.figure(1, figsize=(fig_width,fig_height))
    ax = plt.subplot(1, 1, 1, projection=wcs_cont)

    #plot continuum in color
    im = ax.imshow(data_cont, origin='lower', interpolation='nearest', cmap='viridis',alpha=1.0,transform=ax.get_transform(wcs_cont),vmin=-5.0*noise_cont,vmax=0.3*np.nanmax(data_cont[cut1_y:cut2_y,cut1_x:cut2_x]))

    #add continuum contour levels
    add_continuum_contours(data_cont, cont_levels, 'white', wcs_cont, ax)

    #annotate source names
    annotate_sources(ax,wcs_cont, color='white', fontsize=10)

    #add outflow orientations
    annotate_outflow(ax, wcs_cont, width=2.0)

    #plot properties
    RA = ax.coords[0]
    DEC = ax.coords[1]
    RA.set_axislabel(r'$\alpha$ (J2000)',minpad=0.5)
    DEC.set_major_formatter('dd:mm:ss')
    RA.set_major_formatter('hh:mm:ss.s')
    DEC.set_axislabel(r'$\delta$ (J2000)',minpad=0.3)
    DEC.set_ticklabel(rotation=90., color='black', exclude_overlapping=True)
    RA.set_ticklabel(color='black', exclude_overlapping=True)
    DEC.set_ticks(spacing=10*u.arcsec,color='grey')
    RA.set_ticks(spacing=7.5*u.arcsec,color='grey')
    RA.display_minor_ticks(True)
    DEC.display_minor_ticks(True)
    DEC.set_minor_frequency(5)
    RA.set_minor_frequency(5)

    #add colorbar
    cb=fig.colorbar(im, pad=0.0, shrink=0.855)
    cb.set_label(r'$I_{' + str(wavelength)+ '\mathrm{mm}}$ (mJy\,beam$^{-1}$)')
    cb.ax.yaxis.set_tick_params(color='black',labelcolor='black',direction='out')
    cb.locator = MultipleLocator(10.0)
    cb.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f}"))

    #add linear scale bar (1000 au)
    pix_to_au = distance*delta #au/pixel
    ax.add_artist(Rectangle((cut1_x+10,cut2_y-18), width=1000.0/pix_to_au, height=2, color='k'))
    plt.annotate(r'1000\,au', xy=(cut1_x+10,cut2_y-14), color='k',xycoords='data',va='bottom',ha='left',fontsize=6)

    #add beam ellipse (continuum data)
    ax.add_artist(Ellipse((cut1_x+15, cut1_y+15), width=bmaj_cont/delta, height=bmin_cont/delta, angle=-90.0+bpa_cont, alpha=1.0, edgecolor='black', facecolor='black',lw=0.5,zorder=3))

    #axis limits
    plt.xlim(cut1_x,cut2_x)
    plt.ylim(cut1_y,cut2_y)

    #save plot
    plt.savefig('continuum_' + region + '_' + bb + '.pdf', format='pdf', bbox_inches='tight',pad_inches=0.01)
    #plt.show()
    plt.close(1)


#NOEMA data directory
data_directory = './'

#name of the region
region = 'L1448N'

#distance to the region
distance = 288.0 #pc

#continuum baseband
bb = 'li'

plot_continuum(region, bb, distance, data_directory)
