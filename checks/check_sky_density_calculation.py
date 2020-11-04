from eroML.ensemble import fake_ensemble
from eroML.utils import sky_density
import matplotlib.pyplot as plt
import numpy as np

def gen_random_pos_offset(N, dens=1.):

    rnd = np.random.rand(N)
    rndx = np.sqrt(-np.log(1-rnd)/(np.pi*dens))
    return rndx

e = fake_ensemble(N=500, center=(3,0), width=(1,2), random_pos=True)
c = e.skyCoords()

ra, dec = c.ra.degree, c.dec.degree
#plt.scatter(ra,dec)
#plt.show()

idx, d2d, d3d = c.match_to_catalog_sky(c, nthneighbor=2)

M = 1000

dens = sky_density(e, filter_prop=None).to_array("eligible_sky_density", array_type="array")
eta = np.mean(dens) / 3600
mean_dist = (2*np.pi*eta)**-0.5
print("mean dist: ", mean_dist)
x = np.linspace(0,mean_dist*3, 1000)
ppf = 2*np.pi*x*eta * np.exp(-np.pi*eta*x**2)
plt.hist(d2d.arcsec, range=(0, mean_dist*3), bins=50, density=True, label="random")

rnd = gen_random_pos_offset(M, dens=eta)
plt.hist(rnd, range=(0, mean_dist*3), bins=50, density=True, alpha=0.2, label="generated")

plt.axvline(x=mean_dist, color='r')
plt.plot(x, ppf, lw=3, label="pdf")
plt.legend()
plt.show()
