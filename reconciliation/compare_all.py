import numpy as np 
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations
import copy


def shrinked_fits(gi, ofn):
    ocols = []
    for c in ff0[1].columns:
        #print(c.name)
        ocols.append(
            pyfits.Column(name=c.name, format=c.format, array=ff0[1].data[c.name][gi])
        )
    new_cols = pyfits.ColDefs(ocols)
    hdu = pyfits.BinTableHDU.from_columns(new_cols)
    hdu.writeto(ofn, overwrite=True)
        

def activity_plot(hdu, gi, ax=None, ec='r', ll="", ps=2):
    x = hdu["BP_RP"]
    Fg = 10**(-0.4*hdu["Gmag"])*3.660e-08*720
    y = hdu["Fx"] / Fg
    #ax.scatter(x, np.log10(y), color='0.8', s=3)
    ax.scatter(x[gi], np.log10(y[gi]), edgecolor = ec, label=ll, s=ps, facecolor='none')
    


ff0 = pyfits.open("ero_master.fits")
ff = ff0[1].data

#catalogs = list(fnames.keys())#[::-1]
p_lim = {"SVM":0.01, "Bayes":0.57,"NWAY":0.01}



def one_panel(x,y=None, bg_points=None, ax=None, anno="", color="r",annox=0.5, **kwargs):
    if y is not None and y[0] is not None:
        return one_panel_scatter(x, y, bg_points=bg_points, ax=ax, anno=anno, color=color, annox=annox, **kwargs)
    else:
        return one_panel_hist(x, bg_points=bg_points, ax=ax, anno=anno, color=color, annox=annox, **kwargs)
    #if bg_points is not None: ax.scatter(bg_points[0], bg_points[1], color='0.8', s=3)
    #sc = ax.scatter(x, y, c=color, s=10, vmin=vmin, vmax=vmax)
    #ax.annotate(anno, xy=(annox, 0.1), xycoords="axes fraction")
    #ax.set_xlim(0,4.7)
    #ax.set_ylim(-8,-1)
    #return sc


def one_panel_hist(x,bg_points=None, ax=None, anno="", color="r",annox=0.2, range=(0,1), bins=10, **kwargs):
    if bg_points is not None: ax.hist(bg_points[0], density=True, bins=bins, range=range)#, color='0.8', **kwargs)
    sc = ax.hist(x, density=True, bins=bins, range=range, alpha=0.5)#, color=color)
    ax.annotate(anno, xy=(annox, 0.1), xycoords="axes fraction")
    print("XXX")
    #ax.set_xlim(0,4.7)
    #ax.set_ylim(-8,-1)
    return sc

def one_panel_scatter(x,y=None, bg_points=None, ax=None, anno="", color="r",annox=0.2, xlim=(0,1), ylim=(0,1), **kwargs):
    if bg_points is not None: ax.scatter(bg_points[0], bg_points[1], color='0.8', s=3)
    sc = ax.scatter(x, y, c=color, s=10, vmin=vmin, vmax=vmax)
    ax.annotate(anno, xy=(annox, 0.1), xycoords="axes fraction")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return sc



def one_page(x, y=None, color=None, xlabel="BP-RP", ylabel="log Fx/Fg", **kwargs):
    if y is None:
        if color is None:
            color=len(x)*[color]
        elif len(x) != len(color):
            color=len(x)*[color]
        y = len(x)*[y]
    y = np.array(y)
    color = np.array(color)
    print(len(x), len(color))
        
    gs = fig.add_gridspec(3, 3)

    ax0 = fig.add_subplot(gs[0, :])
    cond0 = (ff["SVM"] > p_lim["SVM"]) & (ff["SVM_ij"] > 0.5)
    cond1 = (ff["Bayes"] > p_lim["Bayes"]) & (ff["Bayes_ij"] > 0.5)
    cond2 = ff["NWAY"] > p_lim["NWAY"]          
    gi_all = np.where(cond0 & cond1 & cond2)[0]
    x_all = x[gi_all]
    if y is not None:
        y_all = y[gi_all]
    else:
        y_all = None
        
    print("all: ",len(gi_all))
    one_panel(x_all, y_all, color=color[gi_all], ax=ax0, annox=0.8, anno="SVM, Bayes, NWAY (%i)" % len(gi_all), **kwargs)
    
    #plt.scatter(x[gi_all], y[gi_all], c=color[gi_all], s=10, vmin=vmin, vmax=vmax)
    #plt.annotate("JS, SF, pcs (%i)" % len(gi_all), xy=(0.8, 0.1), xycoords="axes fraction")
    for k, cat in enumerate(["SVM", "Bayes", "NWAY"]):
        ll = len(np.where((ff[cat] > p_lim[cat]) & (ff[cat+"_ij"] > 0.5))[0])
        print(cat,ll)
        plt.annotate(cat+ ": %i" % ll, xy=(0.8, 0.18+k*.05), xycoords="axes fraction")
        
    #plt.xlim(0,4.7)
    #plt.ylim(-8,-1)

    ax11 = fig.add_subplot(gs[1, 0])
    cond0 = (ff["SVM"] > p_lim["SVM"]) & (ff["SVM_ij"] > 0.5)
    cond1 = (ff["Bayes"] > p_lim["Bayes"]) & (ff["Bayes_ij"] > 0.5)
    cond2 = (ff["NWAY"] < p_lim["NWAY"]) | ((ff["NWAY"]>p_lim["NWAY"]) & (ff["NWAY_ij"] <= 0.5))
    gi = np.where(cond0 & cond1 & cond2)[0]
    sc = one_panel(x[gi], y[gi], bg_points=(x_all, y_all), ax=ax11, anno="SVM, Bayes, !NWAY (%i)" % len(gi), color=color[gi], **kwargs)
    if ylabel=="log Fx/Fg": ax11.plot([0, 2.9], [-3.8, -1], color='r')
    shrinked_fits(gi, ofn="missing.fits")

    if y[0]!=None:
        cax = plt.axes([0.91, 0.08, 0.03, 0.9])
        plt.colorbar(sc, cax=cax)

    ax12 = fig.add_subplot(gs[1, 1])
    cond0 = (ff["SVM"] > p_lim["SVM"]) & (ff["SVM_ij"] > 0.5)
    cond1 = (ff["Bayes"] <= p_lim["Bayes"]) | ((ff["Bayes"]>p_lim["Bayes"]) & (ff["Bayes_ij"] <= 0.5) )
    cond2 = (ff["NWAY"] > p_lim["NWAY"]) & (ff["NWAY_ij"] > 0.5) 
    gi = np.where(cond0 & cond1 & cond2)[0]
    sc = one_panel(x[gi], y[gi], bg_points=(x_all, y_all), ax=ax12, anno="SVM, !Bayes, NWAY (%i)" % len(gi), color=color[gi], **kwargs)

    ax13 = fig.add_subplot(gs[1, 2])
    cond0 = (ff["SVM"] <= p_lim["SVM"]) | ((ff["SVM"]>p_lim["SVM"]) & (ff["SVM_ij"] <= 0.5) )
    cond1 = (ff["Bayes"] > p_lim["Bayes"]) & (ff["Bayes_ij"] > 0.5)
    cond2 = (ff["NWAY"] > p_lim["NWAY"]) & (ff["NWAY_ij"] > 0.5) 
    gi = np.where(cond0 & cond1 & cond2)[0]
    sc = one_panel(x[gi], y[gi], bg_points=(x_all, y_all), ax=ax13, anno="!SVM, Bayes, NWAY (%i)" % len(gi), color=color[gi], **kwargs)

    ax21 = fig.add_subplot(gs[2, 0])
    cond0 = (ff["SVM"] > p_lim["SVM"]) & (ff["SVM_ij"] > 0.5)
    cond1 = (ff["Bayes"] <= p_lim["Bayes"]) | ((ff["Bayes"]>p_lim["Bayes"]) & (ff["Bayes_ij"] <= 0.5) )
    cond2 = (ff["NWAY"] <= p_lim["NWAY"]) | ((ff["NWAY"]>p_lim["NWAY"]) & (ff["NWAY_ij"] <= 0.5) )
    gi = np.where(cond0 & cond1 & cond2)[0]
    sc = one_panel(x[gi], y[gi], bg_points=(x_all, y_all), ax=ax21, anno="SVM, !Bayes, !NWAY (%i)" % len(gi), color=color[gi], **kwargs)
    if ylabel=="log Fx/Fg": ax21.plot([0, 2.9], [-3.8, -1], color='r')

    ax22 = fig.add_subplot(gs[2, 1])
    cond0 = (ff["SVM"] <= p_lim["SVM"]) | ((ff["SVM"]>p_lim["SVM"]) & (ff["SVM_ij"] <= 0.5) )
    cond1 = (ff["Bayes"] <= p_lim["Bayes"]) | ((ff["Bayes"]>p_lim["Bayes"]) & (ff["Bayes_ij"] <= 0.5) )
    cond2 = ff["NWAY"]>p_lim["NWAY"]
    gi = np.where(cond0 & cond1 & cond2)[0]
    sc = one_panel(x[gi], y[gi], bg_points=(x_all, y_all), ax=ax22, anno="!SVM, !Bayes, NWAY (%i)" % len(gi), color=color[gi], **kwargs)
    #shrinked_fits(gi, ofn="pcs_only.fits")
 
    ax23 = fig.add_subplot(gs[2, 2])
    cond0 = (ff["SVM"] <= p_lim["SVM"]) | ((ff["SVM"]>p_lim["SVM"]) & (ff["SVM_ij"] <= 0.5) )
    cond1 = (ff["Bayes"] > p_lim["Bayes"]) & (ff["Bayes_ij"] > 0.5)
    cond2 = (ff["NWAY"] <= p_lim["NWAY"]) | ((ff["NWAY"]>p_lim["NWAY"]) & (ff["NWAY_ij"] <= 0.5) )
    gi = np.where(cond0 & cond1 & cond2)[0]
    sc = one_panel(x[gi], y[gi], bg_points=(x_all, y_all), ax=ax23, anno="!SVM, Bayes, !NWAY (%i)" % len(gi), color=color[gi], **kwargs)

    ax0.set_ylabel(ylabel)
    ax11.set_ylabel(ylabel)
    ax21.set_ylabel(ylabel)

    ax21.set_xlabel(xlabel)
    ax22.set_xlabel(xlabel)
    ax23.set_xlabel(xlabel)


fig = plt.figure(figsize=(12,10))#constrained_layout=True)
fig.subplots_adjust(top=0.98, right=0.90, left=0.08, bottom=0.08)


x = ff["BP_RP"]
Fg = 10**(-0.4*ff["Gmag"])*3.660e-08*720
y = np.log10(ff["Fx"] / Fg)
c = ff["match_dist"]
c = ff["Gmag"]
vmin, vmax= 8, 18
c = np.log10(ff["plx"])
vmin, vmax= -1, 1
one_page(x, y=y, color=c, xlim=(0, 4.8), ylim=(-8,-1))
plt.show()


fig = plt.figure(figsize=(12,10))#constrained_layout=True)
fig.subplots_adjust(top=0.98, right=0.90, left=0.08, bottom=0.08)
x = ff["match_dist"]
y = ff["plx"]
c = ff["Gmag"]
vmin, vmax= 8, 18
one_page(x, y=y, color=c, xlabel="match_dist", ylabel="plx", xlim=(0, 15), ylim=(0,80))
plt.show()


fig = plt.figure(figsize=(12,10))#constrained_layout=True)
fig.subplots_adjust(top=0.98, right=0.90, left=0.08, bottom=0.08)
x = ff["match_dist"]
one_page(x, xlabel="match_dist", ylabel="N", bins=30, range=(0,30))
plt.show()



fig = plt.figure(figsize=(12,10))#constrained_layout=True)
fig.subplots_adjust(top=0.98, right=0.90, left=0.08, bottom=0.08)
x = ff["plx"]
one_page(x, xlabel="plx", ylabel="N", bins=30, range=(0,50))


plt.show()



