from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

fn = "../ero_data/cat2rxs.fits"
ff = pyfits.open(fn)

#plt.hist(ff[1].data["Xerr"]*45, label="xerr")
#plt.hist(ff[1].data["Yerr"]*45, label="yerr")

pos_err = np.sqrt(ff[1].data["Xerr"]**2 + ff[1].data["Yerr"])
print(np.mean(pos_err)*45/2, np.median(pos_err)*45/2)
#plt.show()
