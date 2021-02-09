import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits

def sp(x, y, gi=None, xl="sigma", yl="match_dist"):
        
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.008

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(figsize=(8, 8))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)
    ax_histx = plt.axes(rect_histx, sharex=ax_scatter)
    ax_histx.tick_params(direction='in', labelbottom=False)
    ax_histy = plt.axes(rect_histy, sharey=ax_scatter)
    ax_histy.tick_params(direction='in', labelleft=False)

    ax_scatter.hexbin(x,y, bins='log', mincnt=2, gridsize=50)
    if gi is not None:
        ax_scatter.scatter(x[gi], y[gi], s=3,marker='.', color='r')
    
    x0, x1 = np.quantile(x, 0.05), np.quantile(x, 0.95)
    y0, y1 = np.quantile(y, 0.05), np.quantile(y, 0.95)
    
    ax_scatter.set_xlim(0.8*x0, 1.2*x1)
    ax_scatter.set_ylim(0.8*y0, 1.2*y1)
    ax_scatter.set_xlabel(xl)
    ax_scatter.set_ylabel(yl)
    
    xlim = ax_scatter.get_xlim()
    ylim = ax_scatter.get_ylim()
    print("x range: ",xlim)
    print("y range: ",ylim)
    
    xbins = np.linspace(xlim[0], xlim[1], 20)            
    ax_histx.hist(x, bins=xbins, density=True, label="All", alpha=0.4)
    if gi is not None:
        ax_histx.hist(x[gi], bins=xbins, density=True, label="Real", alpha=0.4)
    ax_histx.legend()
    
    ybins = np.linspace(ylim[0], ylim[1], 20)
    ax_histy.hist(y, bins=ybins, orientation='horizontal', density=True, range=(0,90), alpha=0.4, label="All")
    if gi is not None:
        ax_histy.hist(y[gi], bins=ybins, orientation='horizontal', density=True, range=(0,90), alpha=0.4, label="Real")

    ax_histy.legend()
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_histx.set_ylabel("Normalized Density")
    ax_histy.set_xlabel("Normalized Density")
    
    return ax_scatter


fn = "../../ero_data/merged_eFEDS_EDR3.fits"
ff = pyfits.open(fn)
fd = ff[1].data

gi = np.where(fd["NN"] == 3)[0]
sp(fd["RADEC_ERR"][gi], fd["eligible_sky_density"][gi], xl="RADEC_ERR", yl="Sky Density")
plt.legend()
#plt.xlabel("Match Distance (arcsec)")
plt.show()
exit()

fn0 = "simu.dat"
#fn0 = "../offs2.dat"

#------------------------------
#sigout, pos_off, skdens, nth, cls

dd0 = np.genfromtxt(fn0, unpack=True)

gi0 = np.where((dd0[3]==1) & (dd0[4] ==0 ))[0] # dd[3] = NN, dd[4] = class (0=real)

sp(dd0[2], dd0[1], gi0, xl="Sky Density")

plt.legend()
#plt.xlabel("Match Distance (arcsec)")
plt.show()
