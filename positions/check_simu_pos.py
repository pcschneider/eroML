import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits

def sp(x, y, gi, mfn_xkey="RADEC_ERR", mfn_ykey="match_dist"):
        
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
    ax_scatter.hexbin(x[gi],y[gi], bins='log', mincnt=2)
    #ax_scatter.scatter(x[gi], y[gi],s=1)
    
    # now determine nice limits by hand:
    xbinwidth = (max(x)-min(x))/20
    ybinwidth = (max(y)-min(y))/20
    if xbinwidth < 1e-3: xbinwidth=1
    if ybinwidth < 1e-3: ybinwidth=1
    
    
    
    ax_scatter.set_xlim(min(x), max(x))
    ax_scatter.set_ylim(min(y), max(y))
    xlim = ax_scatter.get_xlim()
    ylim = ax_scatter.get_ylim()

    xbins = np.arange(xlim[0], xlim[1], xbinwidth)
    ybins = np.arange(ylim[0], ylim[1], ybinwidth)
    ax_histx.hist(x, bins=xbins, density=True, label="Simulated")
    #ax_histy.hist(y, bins=ybins, orientation='horizontal')

    #ax_histx.hist(x[gi], bins=xbins, density=True, label="?")
    #ax_histx.hist(ffd["RADEC_ERR"], bins=xbins, density=True, label="Measured", alpha=0.4)
    
    gg = np.where(ffd["NN"] == 1)[0]

    
    ax_histx.hist(ffd[mfn_xkey][gg], bins=xbins, density=True, label="Measured", alpha=0.4)
    ybins=30
    ax_histx.legend()
    
    ax_histy.hist(ffd[mfn_ykey][gg], orientation='horizontal', density=True, range=(0,90), bins=ybins, label="Measured")

    ax_histy.hist(y[gi], bins=ybins, orientation='horizontal', density=True, range=(0,90), alpha=0.4, label="Simulated")
    xx = np.linspace(0,20,100)
    #ax_histx.plot(xx,65*xx**2)

    ax_histy.legend()
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_histx.set_ylabel("Normalized Density")
    ax_histy.set_xlabel("Normalized Density")
    
    return ax_scatter

#oo = np.transpose([sigout, pos_off, skdens, nth, cls])

mfn = "../../ero_data/merged_eFEDS.fits"
sfn = "../offs2.dat"
print("Reading sfn=",sfn," for comparison with mfn=",mfn)

ff = pyfits.open(mfn)
ffd = ff[1].data
dd = np.genfromtxt(sfn, unpack=True)

gi = np.where((dd[3]==1) & (dd[4] >= 0))[0]
ax = sp(dd[0], dd[1], gi)


ax.set_xlabel("RADEC_ERR (arcsec)")
ax.set_ylabel("match dist (arcsec)")
ax.set_xlim(0,19.5)
ax.set_ylim(0,97)
plt.show()

ax = sp(dd[2], dd[1],gi, mfn_xkey="eligible_sky_density")
#ax.scatter(dd[1][gi], dd[0][gi],s=1)



ax.set_xlabel("sk_dens")
ax.set_ylabel("pos_off")
ax.set_xlim(0.1,1.5)
ax.set_ylim(0,100)
plt.show()



exit()
sk = np.unique(dd[1])
print(len(sk))
for s in sk:
    gi = np.where(dd[1] == s)[0]
    plt.hist(dd[0][gi], bins=60, range=(0,30),label=str(s))
    
    plt.hist(dd[0][gi], bins=60, range=(0,30),label=str(s))

plt.legend()
plt.show()
