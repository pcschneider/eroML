from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u


def gen_random_pos_offset(N, dens=1.):
    rnd = np.random.rand(N)
    rndx = np.sqrt(-np.log(1-rnd)/(np.pi*dens))
    return rndx


fn = "../../ero_data/merged_eFEDS.fits"
ff = pyfits.open(fn)
print("Using \'%s\'" % fn)
mn = np.mean(ff[1].data["eligible_sky_density"])
print(r"mean sky density: %6.4e (1/arcmin^2)" % mn)
#print(np.mean(ff[1].data["eligible_sky_density"]))
ff.close()

fn = "../../ero_data/Gaia_merged_eFEDS.fits"
ff = pyfits.open(fn)
gi = np.where(ff[1].data["eligible_Gaia"] == 1)[0]
print("#eligible sources %i" % len(gi))
coord = SkyCoord(ra=ff[1].data["RA"][gi], dec=ff[1].data["Dec"][gi], unit=("degree", "degree"))

around=10
idx, d2d, d3d = match_coordinates_sky(coord, coord, nthneighbor=2)

rnd_dist = gen_random_pos_offset(len(coord), dens=mn/3600)


#uni, cnt = np.unique(idxc, return_counts=True)
print(np.mean(d2d.arcsec), np.mean(rnd_dist))
plt.hist(d2d.arcsec, bins=30, range=(0, 100))
plt.hist(rnd_dist, bins=30, range=(0, 100), alpha=0.5)
plt.show()
