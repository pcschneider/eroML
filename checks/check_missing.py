from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

fn0 = "../../ero_data/SrcCat_V2T_hp.fits"
fn1 = "../../ero_data/major_eFEDS.fits"

ff0 = pyfits.open(fn0)
ff1 = pyfits.open(fn1)

srcIDs0 = ff0[1].data["srcID"]
srcIDs1 = ff1[1].data["srcID"]

print(len(srcIDs0), len(srcIDs1))

idx = np.intersect1d(srcIDs0, srcIDs1)
print(len(idx))
missing = []
cnt = 0
for si in srcIDs0:
    if si not in idx:
        missing.append(si)
        cnt+=1
print(cnt, len(missing))        

#np.savetxt("missing_srcIDs.txt", np.array(missing), fmt='%s')
