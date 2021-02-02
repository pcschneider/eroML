import numpy as np
import matplotlib.pyplot as plt



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
    ax_scatter.hexbin(x,y, bins='log', mincnt=2)
    #ax_scatter.hexbin(x[gi],y[gi], bins='log', mincnt=2)

    #ax_scatter.scatter(x[gi], y[gi],s=1)
    
    # now determine nice limits by hand:
    xbinwidth = (max(x)-min(x))/20
    ybinwidth = (max(y)-min(y))/20
    if xbinwidth < 1e-3: xbinwidth=1
    if ybinwidth < 1e-3: ybinwidth=1
    
    
    
    ax_scatter.scatter(x[gi], y[gi], s=3,marker='.', color='r')
    
    ax_scatter.set_xlim(min(x), max(x))
    ax_scatter.set_ylim(min(y), max(y))
    ax_scatter.set_xlabel(xl)
    ax_scatter.set_ylabel(yl)

    
    
    xlim = ax_scatter.get_xlim()
    ylim = ax_scatter.get_ylim()

    xbins = np.linspace(xlim[0], xlim[1], 20)
    ybins = np.arange(ylim[0], ylim[1], ybinwidth)
    print(xbins)
    ax_histx.hist(x, bins=xbins, density=True, label="All", alpha=0.4)
    #ax_histy.hist(y, bins=ybins, orientation='horizontal')

    ax_histx.hist(x[gi], bins=xbins, density=True, label="Real", alpha=0.4)
    #ax_histx.hist(ffd["RADEC_ERR"], bins=xbins, density=True, label="Measured", alpha=0.4)
    
    #gg = np.where(ffd["NN"] == 1)[0]

    
    #ax_histx.hist(ffd[mfn_xkey][gg], bins=xbins, density=True, label="Measured", alpha=0.4)
    ybins=30
    ax_histx.legend()
    
    #ax_histy.hist(ffd[mfn_ykey][gg], orientation='horizontal', density=True, range=(0,90), bins=ybins, label="Measured")

    ax_histy.hist(y, bins=ybins, orientation='horizontal', density=True, range=(0,90), alpha=0.4, label="All")
    ax_histy.hist(y[gi], bins=ybins, orientation='horizontal', density=True, range=(0,90), alpha=0.4, label="Real")
    xx = np.linspace(0,20,100)
    #ax_histx.plot(xx,65*xx**2)

    ax_histy.legend()
    ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_histx.set_ylabel("Normalized Density")
    ax_histy.set_xlabel("Normalized Density")
    
    return ax_scatter


fn0 = "simu.dat"
fn0 = "../offs2.dat"

#------------------------------
#sigout, pos_off, skdens, nth, cls

dd0 = np.genfromtxt(fn0, unpack=True)

gi0 = np.where((dd0[3]==1) & (dd0[4] ==0 ))[0] # dd[3] = NN, dd[4] = class (0=real)

sp(dd0[2], dd0[1], gi0, xl="Sky Density")

plt.legend()
#plt.xlabel("Match Distance (arcsec)")
plt.show()
