import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from eroML.positions import gen_random_pos_offset, gen_real_pos_offset
from scipy.stats import uniform

ofn = "simu.dat"
NN=3
key=10

N1, sig0, sig1 = 6000, 1, 15    # Real matches
N2, dens0, dens1 = 56000, 0.2, 1.2 # Random matches

#-------------------------------------------------------------

print("Simulating ",N1+3*N2, " sources, with a real fraction of ",N1/N2)
sigs = uniform.rvs(size=N1)*(sig1-sig0)+sig0



offs1 = gen_real_pos_offset(N=N1, sigma=sigs)

lc = (dens0/dens1)**2


dens = uniform.rvs(size=N2+N1, loc=lc, scale=1-lc)**0.5 * dens1 # Scaling proportional to density
#dens = uniform.rvs(size=N2+N1, loc=dens0, scale=(dens1-dens0)) # Uniform scaling in density
dens_real = dens[0:N1]
dens_rand = dens[N1::]

dens_real = uniform.rvs(size=N1, loc=dens0, scale=(dens1-dens0))
dens_rand = uniform.rvs(size=N2, loc=lc, scale=1-lc)**0.5 * dens1


print(len(dens_rand), N2, len(dens))

#plt.hist(dens, density=True)
#plt.plot([dens0, dens1], [dens0/4.5, dens1/4.5])
#plt.show()


offs2, dens2, group, nth  = gen_random_pos_offset(N=N2, dens=dens_rand/3600)

plt.hist(offs1, label="real", alpha=0.5, density=True, range=(0, 30))

plt.hist(offs2, label="random",alpha=0.5, density=True, range=(0, 30))

sigs2 = uniform.rvs(size=N2)*(sig1-sig0)+sig0


s = np.concatenate((sigs, np.repeat(sigs2,NN))).flatten()
o = np.concatenate((offs1, offs2)).flatten()
d = np.concatenate((dens_real, dens2*3600)).flatten()
n = np.concatenate((np.ones(N1), nth)).flatten()
c = np.concatenate((np.zeros(N1), np.repeat(np.ones(N2),NN)) ).flatten()
k = np.repeat([key], len(c))
for x in [s,o,d,n,c,k]:
    print(np.shape(x))    
    
#sigout, pos_off, skdens, nth, cls
out = np.transpose([s, o, d, n, c, k])
print("Writing to ",ofn,np.shape(out))
np.savetxt(ofn, out)


