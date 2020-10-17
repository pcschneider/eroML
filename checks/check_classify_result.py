import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u

fn = "../classify/major_proba.fits"
fn = "../classify/random_proba.fits"
#fn = "../../ero_data/training_eFEDS.fits"
ff = pyfits.open(fn)

if "category" in ff[1].data.columns.names:
    gi = np.where((ff[1].data["category"]==0) & (ff[1].data["NN"]<2) )[0]
    plt.hist(ff[1].data["pos"][gi])
    plt.xlabel("pos")
    plt.show()             
    gi = np.where((ff[1].data["category"]==0) & (ff[1].data["NN"]<2)  & (ff[1].data["pos"]<110))[0]
else:
    gi = np.where((ff[1].data["predicted"]==0) & (ff[1].data["NN"]<2)  & (ff[1].data["pos"]<9.5))[0]
    
print(len(gi))
dist = 1000/10**ff[1].data["log_plx"][gi]
iii = np.where(dist>700)[0]
print(len(iii))
nbins=20
logbins = np.logspace(np.log10(10),np.log10(5000),nbins)
plt.hist(dist, bins=logbins)
plt.xlabel("d (pc)")
plt.ylabel("N")
plt.xscale("log")
plt.show()


plt.scatter(ff[1].data["bp_rp"], ff[1].data["logFxFg"])
plt.scatter(ff[1].data["bp_rp"][gi], ff[1].data["logFxFg"][gi])
plt.xlabel("bp_rp")
plt.ylabel("logFxFg")
plt.show()

plt.scatter(ff[1].data["pos"], ff[1].data["logFxFg"])
plt.scatter(ff[1].data["pos"][gi], ff[1].data["logFxFg"][gi])
plt.xlabel("pos")
plt.ylabel("logFxFg")
plt.show()


plt.scatter(ff[1].data["pos"], ff[1].data["log_plx"])
plt.scatter(ff[1].data["pos"][gi], ff[1].data["log_plx"][gi])
plt.xlabel("pos")
plt.ylabel("log_plx")
plt.show()

