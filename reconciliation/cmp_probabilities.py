import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
plt.rcParams.update({'font.size': 16})

#fn0 = "major4classify_eFEDS_HamStar.fits"
fn0 = "major_eFEDS_classified_HamStar.fits"
ff0 = pyfits.open(fn0)
fd0 = ff0[1].data
srcIDs0 = fd0["original_srcID"]

fn1 = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"
ff1 = pyfits.open(fn1)
fd1 = ff1[1].data

c = "ero_ID"
srcIDs1 = np.array([str("ML%05i" % int(d)) for d in fd1[c]])
print(srcIDs0, srcIDs1)

idx, i0, i1 = np.intersect1d(srcIDs0, srcIDs1, return_indices=True)
print(srcIDs0[i0], srcIDs1[i1])

plt.scatter(fd0["svm_prob"][i0], fd1["p_ij"][i1])
plt.show()

pp = [1.98372460, -4.53855281,  3.53929735, -4.44042056e-3]


fig = plt.figure()
fig.subplots_adjust(top=0.98, right=0.94, left=0.12, bottom=0.13)


svm = np.polyval(pp, fd0["svm_prob"][i0])
hb = plt.hexbin(svm, fd1["p_ij"][i1], mincnt=3, bins="log", gridsize=20)
cb = plt.colorbar(hb)
cb.set_label("Number per bin")

plt.xlabel("SVM")
plt.ylabel("Bayes")
plt.show()

si0 = np.argsort(fd0["svm_prob"][i0])[::-1]#[0:100]
si00 = np.argsort(svm)[::-1]#[0:100]
si1 = np.argsort(fd1["p_ij"][i1])[::-1]#[0:100]
plt.scatter(np.cumsum(1-fd0["svm_prob"][i0][si0]), np.nancumsum(1-fd1["p_stellar"][i1][si1]), label="1")
plt.scatter(np.cumsum(1-svm[si00]), np.nancumsum(1-fd1["p_stellar"][i1][si1]), label="2")
plt.legend()
plt.show()
