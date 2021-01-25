import numpy as np
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt

N0 = 300000
ps = 0.1


ff = pyfits.open("../../ero_data/merged_eFEDS.fits")
ffd = ff[1].data
gi = np.where(ffd["NN"] == 1)[0]
sig = (0.6 * ffd["RADEC_ERR"][gi]) /60/60 * np.pi/180
md = ffd["match_dist"][gi] /60/60 * np.pi/180

dd = np.genfromtxt("../offs.dat", unpack=True)
    #oo = np.transpose([sigout, pos_off, skdens, nth, cls])
gi = np.where(dd[3] == 1)[0]
sig = (0.6 * dd[0][gi]) /60/60 * np.pi/180
md = dd[1][gi] /60/60 * np.pi/180




tmp1 = 2/sig**2*np.exp(-md**2 / (2* sig**2))                     
tmp2 = ps/(1-ps) * 1/N0 * tmp1            

p_stellar = tmp2 / (1+tmp2)

si = np.argsort(p_stellar)[::-1]

for j in range(350, 800):
    Nmax = j
    star_indices = si[0:Nmax ]
    plt.hist(dd[4][star_indices])
    ns = np.sum(dd[4][star_indices])
    print(j," -> Non stars: ", ns, ns/Nmax)
plt.show()

plt.hist(p_stellar[si][0:Nmax])
NN = ffd["NN"][gi][si][0:Nmax]
#print(NN)
plt.show()
plt.hist(NN)
plt.show()
                    
