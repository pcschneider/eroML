from eroML.tile import file4
from astropy.io import fits as pyfits
import numpy as np

fn = file4("major", cconfig="eFEDS_EDR3.ini")
print(fn)
ff = pyfits.open(fn)
fd = ff[1].data

gi = np.where(fd["NN"] == 1)[0]
print("Found nearest neighbour to ",len(gi))
mm = np.unique(fd["original_srcID"])
mi = np.array([str("%i" % int(x[2:])) for x in mm])
print("   ",len(mm))
print(mi)
print()

#fn = '../ero_data/eFEDS_c001_hs.fits'
fn = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"

ff = pyfits.open(fn)
fd = ff[1].data
try:
    ll = np.unique(fd["ero_NAME"])
    print("unique ero_name: ",len(ll))
    print(np.setdiff1d(mm,ll))
except:
    pass

print(len(fd))
ee = ["fd[\"RADEC_ERR\"] > 0", "fd[\"DET_LIKE_0\"]>6", "fd[\"EXT_LIKE\"]<6"]
a = ""
for e in ee:
    print(e)
    gi = np.where(eval(e))[0]
    print(len(gi))
    print()

a = " & ".join(["("+x+")" for x in ee])
print(a)
gi = np.where(eval(a))[0]
print(len(gi))
    
    
