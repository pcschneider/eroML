import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from eroML.positions import calc_sigma_from_RADEC_ERR 

def sp(x, y, gi, xl="sigma", yl="match_dist"):
        
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

    # the scatter plot:
    #ax_scatter.scatter(x, y)
    
    #ax_scatter.hexbin(x,y, bins='log')
    print(x)
    print(y)
    print(gi)
    h = ax_scatter.hexbin(x,y, bins='log', mincnt=2)
    xy = h.get_offsets()
    xi, yi = xy[::,0], xy[::,1]
    zi = h.get_array()
    
    
    h = np.histogram2d(x,y, bins=20, range=((0, 12),(0,60)))
    xi = (h[1][1::] + h[1][0:-1])/2
    yi = (h[2][1::] + h[2][0:-1])/2
    zi = h[0]
    #ax_scatter.imshow(np.transpose(zi))
    ax_scatter.contour(xi, yi, np.transpose(zi), colors='w')
    #print("h",xi, yi, zi)
    #ax_scatter.tricontour(xi, yi, zi, levels=3)
             #norm=plt.Normalize(vmax=abs(zi).max(), vmin=-abs(zi).max()))

    #ax_scatter.hexbin(x[gi],y[gi], bins='log', mincnt=2)

    #ax_scatter.scatter(x[gi], y[gi],s=1)
    
    # now determine nice limits by hand:
    xbinwidth = (max(x)-min(x))/20
    ybinwidth = (max(y)-min(y))/20
    if xbinwidth < 1e-3: xbinwidth=1
    if ybinwidth < 1e-3: ybinwidth=1
    
    
    
    ax_scatter.scatter(x[gi], y[gi], s=3,marker='.', color='r', alpha=0.2)
    
    ax_scatter.set_xlim(min(x), max(x))
    ax_scatter.set_ylim(min(y), max(y))
    ax_scatter.set_xlabel(xl)
    ax_scatter.set_ylabel(yl)

    
    
    xlim = ax_scatter.get_xlim()
    ylim = ax_scatter.get_ylim()

    xbins = np.linspace(xlim[0], xlim[1], 20)
    ybins = np.arange(ylim[0], ylim[1], ybinwidth)
    print(xbins)
    xbins=40
    ax_histx.hist(x, bins=xbins, density=True, label="All", alpha=0.4)
    #ax_histy.hist(y, bins=ybins, orientation='horizontal')

    ax_histx.hist(x[gi], bins=xbins, density=True, label="Real", alpha=0.4)
    #ax_histx.hist(ffd["RADEC_ERR"], bins=xbins, density=True, label="Measured", alpha=0.4)
    
    #gg = np.where(ffd["NN"] == 1)[0]

    
    #ax_histx.hist(ffd[mfn_xkey][gg], bins=xbins, density=True, label="Measured", alpha=0.4)
    ybins=61
    ax_histx.legend()
    
    #ax_histy.hist(ffd[mfn_ykey][gg], orientation='horizontal', density=True, range=(0,90), bins=ybins, label="Measured")

    ax_histy.hist(y, bins=ybins, orientation='horizontal', range=(0,90), alpha=0.4, label="All", density=True)
    #ax_histy.hist(y[gi], bins=ybins, orientation='horizontal', range=(0,90), alpha=0.4, label="Real")
    xx = np.linspace(0,20,100)
    #ax_histx.plot(xx,65*xx**2)

    
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_histx.set_ylabel("Normalized Density")
    ax_histy.set_xlabel("Normalized Density")


    fn = "../../ero_data/merged_major_eFEDS_EDR3.fits"
    ff = pyfits.open(fn)
    fd = ff[1].data
    gi = np.where((fd["NN"] == 1) & (fd["match_dist"]<60))[0]
    print("Number of nearest neighbours within 60 arcsec:", len(gi))
    x, y = fd["RADEC_ERR"][gi], fd["match_dist"][gi]
    x = calc_sigma_from_RADEC_ERR(x)
    h = np.histogram2d(x,y, bins=20, range=((0, 12),(0,60)))
    xi = (h[1][1::] + h[1][0:-1])/2
    yi = (h[2][1::] + h[2][0:-1])/2
    zi = h[0]
    #ax_scatter.imshow(np.transpose(zi))
    ax_scatter.contour(xi, yi, np.transpose(zi), colors='r', linestyles=':')
    
    ax_histy.hist(y, bins=ybins, density=True, zorder=0, range=(0,90), orientation='horizontal', label="Measured")
    ax_histy.legend(loc='upper left')
    return ax_scatter


fn0 = "simu.dat"
fn0 = "../offs2.dat"
fn0 = "../offs5.dat"

#------------------------------
#sigout, pos_off, skdens, nth, cls

dd0 = np.genfromtxt(fn0, unpack=True)

giNN = np.where(dd0[3]==1)[0]
print("Number of nearest neighbours: ",len(giNN))
gi0 = np.where(dd0[4][giNN] ==0 )[0] # dd[3] = NN, dd[4] = class (0=real)
#gi0 = np.where(dd0[3]==1)[0]

#sp(dd0[2][giNN], dd0[1][giNN], gi0, xl="Sky Density")
ax = sp(dd0[0][giNN], dd0[1][giNN], gi0, xl="sigma")
ax.set_xlim(0.9,6.2)

plt.legend()
#plt.xlabel("Match Distance (arcsec)")
plt.show()
