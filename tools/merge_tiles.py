from astropy.io import fits as pyfits
import numpy as np
from eroML.ensemble import from_fits, to_fits

ifn = "../../ero_data/merged_eFEDS.fits"
prefix="ID3"

ff = pyfits.open(ifn)
hpx = np.unique(ff[1].data["healpix"])
print(hpx)
ff.close()

fn = "../../ero_data/Gaia_ID2_nside32_"+str(hpx[0])+".fits"
e0 = from_fits(fn, maxN=100)
print(len(e0))

for hp in hpx[1:]:
    fn = "../../ero_data/Gaia_ID2_nside32_"+str(hp)+".fits"
    e = from_fits(fn)
    e0.append(e, duplicates='ignore')
    print(len(e0))

to_fits(e0, ofn="x.fits", overwrite=True)
