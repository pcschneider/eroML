from astropy.io import fits as pyfits
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

fn = "merged_training.fits"
fn = "merged_random.fits"
fn = "../ero_data/merged_major.fits"
ff = pyfits.open(fn)
NSIDE=32
NPIX = hp.nside2npix(NSIDE)
#print(NPIX)
#hps = np.unique(ff[1].data["healpix"])
hps = ff[1].data["healpix"]

uni, cnt = np.unique(hps, return_counts=True)
print(len(hps), len(uni))

m = np.ones(NPIX)
m[uni] = cnt

hp.mollview(m, norm="log")

plt.figure()

ra, dec = ff[1].data["RA"], ff[1].data["Dec"]
c = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")
theta, phi = c.galactic.l.degree, c.galactic.b.degree



theta, phi = ff[1].data["RA"], ff[1].data["Dec"]
theta[theta<0]+=360
phi+=90
hb = plt.gca().hexbin(theta, phi, gridsize=70, mincnt=2, bins='log')#, cmap=plt.cm.BuGn_r)

plt.figure()
gi = np.random.choice(len(theta), 10000)
plt.scatter(theta[gi], phi[gi], c= 1./ff[1].data["sky_density"][gi], alpha=0.5, vmin=0, vmax=200)

plt.show()
