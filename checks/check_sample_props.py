import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits

fn = "random_eFEDS_classified.fits"
#fn = "major_eFEDS_classified.fits"
ff = pyfits.open(fn)
fd = ff[1].data

color = fd["bp_rp"]
FxFg = fd["FxFg"]
log_plx = fd["log_plx"]
md = fd["match_dist"]
sig = fd["RADEC_sigma"]
gi = np.where(fd["category"] == 0)[0]
print("#stars: ",len(gi))
plt.scatter(color, FxFg, label="All", c=log_plx, vmin=-1, vmax=2.)
plt.scatter(color[gi], FxFg[gi], label="Stars", fc="None", ec="r", s=24)
cb = plt.colorbar()
cb.set_label("log plx")
#plt.yscale("log")
#plt.ylim(2e-7, 1e-1)
plt.legend()
plt.xlabel("BP-RP (mag)")
plt.ylabel("Fx/Fbol")
plt.show()

plt.hist(md/sig)
plt.hist(md[gi]/sig[gi])
x = np.linspace(0,5)
plt.plot(x, len(gi)*0.5*x*np.exp(-0.5*x**2))
plt.show()

plt.hist(log_plx)
plt.hist(log_plx[gi])
#x = np.linspace(0,5)
#plt.plot(x, len(gi)*0.5*x*np.exp(-0.5*x**2))
plt.show()
