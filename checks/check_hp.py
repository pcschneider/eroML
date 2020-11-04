import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

nside=32
hpix = 100

Npix = hp.nside2npix(32)
print("Npix", Npix)

i = 10
pix_center = hp.pix2ang(nside, i, nest=True)
a, b = pix_center[1]/np.pi*180, 90. - pix_center[0]/np.pi*180
print("RA", a, " Dec", b)
phi = 90-b
theta = a
c = hp.pixelfunc.ang2pix(nside, phi/180*np.pi, theta/180*np.pi, nest=True)
print(i, c)

ra, dec = [], []

resol = hp.nside2resol(nside, arcmin=True)
for i in range(Npix):
    pix_center = hp.pix2ang(nside, i, nest=True)
    a, b = pix_center[1]/np.pi*180, 90. - pix_center[0]/np.pi*180
    ra.append(a)
    dec.append(b)

ra = np.array(ra)
dec = np.array(dec)

ra[ra<0]+=360
phi = 90-dec
theta = ra
    
x = hp.pixelfunc.ang2pix(nside, phi/180*np.pi, theta/180*np.pi, nest=True)
print(len(np.unique(x)))
plt.scatter(ra,dec)

#plt.figure()
#hp.mollview(np.array(len(a)*[1]))
hp.mollview(x)
plt.show()

#print(ra, dec)
