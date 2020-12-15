from astropy.io import fits as pyfits
import numpy as np
from eroML.ensemble import from_fits, to_fits
import matplotlib.pyplot as plt

ifn = "../../ero_data/merged_eFEDS.fits"
prefix="ID3"

ff = pyfits.open(ifn)
hpx = np.unique(ff[1].data["healpix"])
print(hpx)
ff.close()

fn = "../../ero_data/Gaia_ID2_nside32_"+str(hpx[0])+".fits"
#fn = "../../ero_data/Gaia_EDR3_nside32_"+str(hpx[0])+".fits"

e0 = from_fits(fn)
print(len(np.unique(e0.srcIDs())))
print(len(e0))
o = len(e0)
#exit()

srcIDs = pyfits.open(fn)[1].data["srcID"]

for hp in hpx[1:]:
    fn = "../../ero_data/Gaia_EDR3_nside32_"+str(hp)+".fits"
    fn = "../../ero_data/Gaia_ID2_nside32_"+str(hp)+".fits"
    e = from_fits(fn)
    ee = len(e)
    ra, dec = e.to_array("RA", array_type="array"), e.to_array("Dec", array_type="array")
    e0.append(e, duplicates='ignore')
    print(o,"+",ee, " = ", len(e0), "(",o+ee,")")
    o+=len(e)
    srcIDs = np.concatenate((srcIDs, pyfits.open(fn)[1].data["srcID"]))
    print(len(np.unique(srcIDs)))
    print()
    
    plt.scatter(ra, dec, label=str(hp), alpha=0.3)



print(len(srcIDs), len(np.unique(srcIDs)))  
ra, dec = e0.to_array("RA", array_type="array"), e0.to_array("Dec", array_type="array")
plt.scatter(ra, dec, label="merged", s=5)

plt.legend()
#exit()
#for si in np.unique(srcIDs):
    #print(si, si in e0.srcIDs())
    
    

to_fits(e0, ofn="x.fits", overwrite=True)
plt.show()
