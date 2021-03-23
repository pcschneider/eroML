import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
from PyAstronomy.pyasl import binningx0dt
fn = "major_eFEDS_classified_HamStar.fits"
ff = pyfits.open(fn)
fd = ff[1].data

rfn = "random_eFEDS_classified_HamStar.fits"
rff = pyfits.open(rfn)
rfd = rff[1].data


a, b = fd["svm_prob"], fd["category"]
si = np.where(b==0)[0]
mn = min(a[si])

last = 0
pps = []
si = np.argsort(a)[::-1]
for i in si:
    lgreal = np.where(a>=a[i])[0]
    lgrand = np.where(rfd["svm_prob"]>a[i])[0]
    pred = np.sum(1-a[lgreal])
    Nrand = len(lgrand)/10.
    delta_rand = Nrand - last
    preal =  1-delta_rand
    #print(a[i], len(lgreal),  pred, Nrand, " - ",preal)
    last = Nrand
    pps.append(preal)
    
plt.plot(a[si], pps)
x,dt = binningx0dt(a[si], pps, reduceBy=50)
plt.plot(x[::,0], x[::,1])
pp = np.polyfit(x[::,0], x[::,1], 3)
print("poly",pp)
plt.plot(a[si], np.polyval(pp, a[si]))
plt.show()

x = np.linspace(0.02, 0.95, 1000)


y, z, e, n = [], [], [],   []
for xx in x:
    tmp = np.where(fd["svm_prob"]>xx)[0]
    y.append(len(tmp))
    z.append(len(np.where(rfd["svm_prob"]>xx)[0])/10.)
    e.append(np.sum(1 - fd["svm_prob"][tmp]))
    n.append(np.sum(1 - np.polyval(pp, fd["svm_prob"][tmp])))
y = np.array(y)#/1.05
z = np.array(z)
e = np.array(e)
n = np.array(n)

plt.plot(z, e)
plt.plot(z, n, label="n")
plt.plot([0, 1200],[0, 1200])
plt.legend()
plt.show()

plt.plot(x, y/2060, label="predicted")
plt.plot(x, z/2060, label="random")
plt.plot(x, (y-z)/2060, label="predicted - random")
plt.plot([0,0.5],[1,1])
plt.legend()
plt.show()
