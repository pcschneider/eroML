from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt

fn0 = "../ero_data/eFEDS_c001_main_ctp_star_1.0.fits"
#fn1 = "major_eFEDS_classified_HamStar.fits"
fn1 = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"

ff0 = pyfits.open(fn0)
ff1 = pyfits.open(fn1)

fd0 = ff0[1].data
fd1 = ff1[1].data

srcIDs0 = fd0["ero_ID"]
srcIDs1 = fd1["ero_ID"]

idx, i0, i1 = np.intersect1d(srcIDs0, srcIDs1, return_indices=True)
print(len(np.unique(idx)))
#for i in range(len(i0)):
    #print(i, " - ",idx[i], srcIDs0[i0[i]], srcIDs1[i1[i]], fd0["p_stellar"][i0[i]],  fd1["p_stellar"][i1[i]])
#print(len(idx))

#plt.scatter(fd0["p_stellar"][i0],  fd1["p_stellar"][i1])
#plt.plot([0.5, 1.0], [0.5,1], color='r')
#plt.show()

#plt.scatter(fd0["p_ij"][i0],  fd1["p_ij"][i1])
#plt.plot([0.5, 1.0], [0.5,1], color='r')
#plt.show()


ll = 0.5882
gi0 = np.where((fd0["p_stellar"] > ll) & (fd0["p_ij"]>0.05))[0]
gi1 = np.where((fd1["p_stellar"] > ll) & (fd1["p_ij"]>0.05))[0]
print(len(gi0), len(gi1))
print(len(np.unique(srcIDs0[gi0])))
print(len(np.unique(srcIDs1[gi1])))

print(np.sum(1-fd0["p_stellar"][gi0]))
print(np.sum(1-fd1["p_stellar"][gi1]))
print(np.sum(1-fd1["p_ij"][gi1]))
