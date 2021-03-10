from eroML.utils import major_set

fn0, fn1 =  "../ero_data/ero_eFEDS_EDR3_nside32_6827.fits", "../ero_data/Gaia_EDR3_nside32_6827.fits"

x = major_set(fn0, fn1, "x.fits", overwrite=True)
print(len(x), len(x)/3)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits


fn = "../../ero_data/major_eFEDS_EDR3_nside32_6827.fits"
fn = "../../ero_data/Gaia_EDR3_nside32_6827.fits"
fn = "x.fits"
ff = pyfits.open(fn)
fd = ff[1].data
#plt.hist(np.log10(fd["FxFg"]))

plt.scatter(fd["bp_rp"], np.log10(fd["FxFg"]), c=fd["match_dist"], alpha=0.1)
gi = np.where(fd["match_dist"] < 10)[0]
plt.scatter(fd["bp_rp"][gi], np.log10(fd["FxFg"][gi]), c=fd["match_dist"][gi])
plt.colorbar()

activity_poly = [-3.22, 3.6/5.5]
x = np.linspace(0,4, 100)
ys = activity_poly[1]
y = x*ys + activity_poly[0]
plt.plot(x, y, color='k')

activity_poly = [-3.22+0.5, 3.6/5.5]
x = np.linspace(0,4, 100)
ys = activity_poly[1]
y = x*ys + activity_poly[0]
plt.plot(x, y, color='k')


plt.show()
