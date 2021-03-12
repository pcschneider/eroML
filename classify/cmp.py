import numpy as np
from astropy.io import fits as pyfits

fn  = "major_eFEDS_classified.fits"
ff = pyfits.open(fn)
fd = ff[1].data
gi = np.where(fd["category"] == 0)[0]
s_names = fd["original_srcID"][gi]
print(s_names, len(s_names))


fn1 = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"
ff1 = pyfits.open(fn1)
fd1 = ff1[1].data

gi = np.where(fd1["p_stellar"] > 0.735)[0]
print(len(gi))
#b_names = fd1["ero_NAME"][gi]
b_names = np.array([str("ML%05i" % int(d)) for d in fd1["ero_ID"][gi]])
print(b_names)

ovl = np.in1d(s_names, b_names) # s_name in b_names
print(ovl, len(ovl), np.sum(ovl))

