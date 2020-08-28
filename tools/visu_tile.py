from astropy.io import fits as pyfits
import matplotlib.pyplot as plt

directory="../../ero_data/"
fn0 = directory+"ero_rID1_nside32_9027.fits"
fn1 = directory+"Gaia_rID1_nside32_9027.fits"
fn2 = directory+"major_ID2_nside32_9027.fits"

ff0 = pyfits.open(fn0)
ff1 = pyfits.open(fn1)
ff2 = pyfits.open(fn2)

plt.scatter(ff0[1].data["RA_CORR"], ff0[1].data["DEC_CORR"], color='r', label="ero")
plt.scatter(ff1[1].data["RA"], ff1[1].data["DEC"], color='b', label="Gaia")
plt.scatter(ff2[1].data["RA"], ff2[1].data["DEC"], color='g', marker='s', fc="none", label="matched")


plt.legend()
plt.show()
