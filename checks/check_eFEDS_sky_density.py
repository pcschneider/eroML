from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u

fn = "../../ero_data/merged_eFEDS.fits"
ff = pyfits.open(fn)
print("Using \'%s\'" % fn)
gi = np.where(ff[1].data["NN"] == 1)[0]
coord = SkyCoord(ra=ff[1].data["RA"][gi], dec=ff[1].data["Dec"][gi], unit=("degree", "degree"))

around=10
idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, around*u.arcmin)
uni, cnt = np.unique(idxc, return_counts=True)

print(len(gi), len(uni), np.mean(cnt), np.median(cnt))
print(r"mean sky density: %6.4e (1/arcmin^2)" % (np.mean(cnt)/(np.pi*(around)**2)))
print(np.mean(ff[1].data["eligible_sky_density"]))
plt.hist(cnt)
plt.show()
