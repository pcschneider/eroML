import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits

fn = "major_eFEDS_classified_HamStar.fits"
ff = pyfits.open(fn)
fd = ff[1].data

rfn = "random_eFEDS_classified_HamStar.fits"
rff = pyfits.open(rfn)
rfd = rff[1].data

bfn = "../ero_data/eFEDS_c001_main_ctp_star_1.0.fits"
#bfn = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"

bff = pyfits.open(bfn)
bfd = bff[1].data

x = np.linspace(0.02, 0.95, 1000)
uids = np.unique(bfd["ero_ID"]).tolist()
#print(uids)
ui = [uids.index(a) for a in uids]
print(ui)
#exit()


num_ident = []
num_missed, num_spurious = [], []
for xx in x:
    gi = np.where(bfd["p_stellar"][ui]>xx)[0]
    gi2 = np.where(bfd["p_stellar"][ui]<xx)[0]
    num_ident.append(len(gi))
    #np.nansum(p_stellar[np.logical_and(no_mult, np.logical_not(temp))]/100)
    num_missed.append(np.nansum(bfd["p_stellar"][gi2]))
    #np.nansum(1-p_stellar[np.logical_and(no_mult, temp)]/100)
    num_spurious.append(np.nansum(1 - bfd["p_stellar"][gi]))
    if num_ident[-1] >=2020:
        print(xx, len(gi), num_missed[-1], num_spurious[-1])

num_ident = np.array(num_ident)

num_missed = np.array(num_missed)
num_spurious = np.array(num_spurious)

completeness = num_ident/(num_ident+num_missed)
completeness = (num_ident-num_spurious)/2060
reliability = (num_ident-num_spurious)/num_ident
plt.plot(x, num_ident, label="above limit")
plt.plot(x, num_missed, label="missed")
plt.plot(x, num_spurious, label="spurious")
plt.plot(x, completeness, label="Completeness")
plt.plot(x, reliability, label="Reliability")
plt.legend()
plt.show()
#exit()

a, b = fd["svm_prob"], fd["category"]
si = np.where(b==0)[0]
mn = min(a[si])


y, z = [],   []
for xx in x:
    y.append(len(np.where(fd["svm_prob"]>xx)[0]))
    z.append(len(np.where(rfd["svm_prob"]>xx)[0]))
    
y = np.array(y)#/1.05
z = np.array(z)

plt.plot(x, y/2060, label="predicted")
plt.plot(x, z/10./2060, label="random")
plt.plot(x, (y-z/10.)/2060, label="predicted - random")
plt.plot([0,0.5],[1,1])
plt.legend()
plt.show()

xx = (y-z/10)/2060 + (x-0.235)*0.4
yy = (y-z/10.)/y
plt.plot(x, (y-z/10)/2060 + (x-0.235)*0.4, label="completeness")
plt.plot(x, (y-z/10.)/y, label="reliability")
plt.plot([0,1],[0.9,0.9])
plt.legend()
plt.show()

plt.plot(yy, xx, label="SVC")
plt.plot(reliability, completeness, label="Bayesian")
plt.scatter([0.896],[0.896], marker='x', color='r')
plt.xlim(0.85, 0.99)
plt.legend()
plt.xlabel("Reliability")
plt.ylabel("Completeness")
plt.show()

    
print(mn, len(np.where(a>=mn)[0]))
print(len(np.where(a>0)[0]))
print()

a, b = rfd["svm_prob"], rfd["category"]
si = np.where(b==0)[0]
mn = min(a[si])
print(mn, len(np.where(a>=mn)[0]))
print(len(np.where(a>0)[0]))


plt.scatter(a, b)
plt.show()
