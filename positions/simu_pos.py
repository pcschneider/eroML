import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from eroML.positions import gen_random_pos_offset, gen_real_pos_offset, generate_random_densities
from scipy.stats import uniform, poisson
    
N = 30000
NN=N*5
max_dist = 60
dens0, dens1 = 0.95, 9.55
dens = uniform.rvs(size=N)*(dens1-dens0)+dens0

#offs = max_dist*np.random.rand(N,NN*2)**0.5
offs = max_dist*np.random.rand(N*20)**0.5

Nexp = dens/3600*max_dist**2*np.pi 
Nmeaus = poisson.rvs(Nexp, size=N)

#print(Nexp)

grp = np.arange(N)

skdens = np.repeat(dens, Nmeaus)
offs_idx = np.repeat(grp, Nmeaus)
pos_off = offs[0:len(offs_idx)]
sigout = [np.nan]*len(offs_idx)
nth = np.zeros(len(offs_idx))
j = 0
print("calculating 'nth'")
for i, g in enumerate(grp):
    gi = np.where(offs_idx == g)[0]
    n = np.argsort(np.argsort(pos_off[gi]))+1
    nth[j:j+len(gi)] = n
    #print(g, len(gi), Nmeaus[i], pos_off[gi], n)
    j+=len(gi)
print("                done.")    
    #nth.append(n)
#print(nth)    
#exit()    
#nth = [np.nan]*len(offs_idx)


np.savetxt("xxx.dat", np.transpose([sigout, pos_off, skdens, nth, offs_idx]))



#plt.hist(pos_off)
#plt.show()


plt.hexbin(skdens, pos_off)
plt.show()

#for i in [1,2,3,4]:
    #gi = np.where(nth == i)[0]
    #plt.hist(pos_off[gi], label=str(i))
    #print(i, len(gi), np.mean(pos_off[gi]), np.median(pos_off[gi]))
#plt.show()    



delta = (dens1-dens0)/20
print(delta)
for d in np.linspace(dens0, dens1, 11):
    gi = np.where((abs(skdens-d)<delta) & (nth == 1))[0]
    #gi = np.where(abs(skdens-d)<delta)[0]
    plt.hist(pos_off[gi], label=str(d))
    print(d, len(gi), np.mean(pos_off[gi]), np.median(pos_off[gi]))
plt.show()    

exit()

#print(offs_meaus)
print(len(pos_offs))
plt.hist(offs_meaus)
plt.show()
exit()

plt.hexbin(Nexp, Nmeaus, gridsize=10)
plt.show()
exit()
plt.hist(offs)
plt.show()
print(np.shape(offs))
exit()
    
    
def gen_random(dens, max_dist=60):
    """
    Parameters
    ----------
    max_dist : float (in arcsec)
    """
    #max_dist = np.repeat(np.sqrt(NN/np.pi/dens), NN*2).reshape((N,NN*2))
    #print(np.shape(max_dist), max_dist)
    dd = np.repeat(dens,NN*2).reshape((N,NN*2))
    gg = np.repeat(np.arange(len(dens)),NN*2).reshape((N,NN*2))
    #print(N, np.shape(max_dist), max(max_dist))
    NNsimu = np.random.poisson(lam=np.array(len(dens)*[NN]))
    #print(np.shape(NNsimu), NNsimu)
    offs = max_dist*np.random.rand(N,NN*2)**0.5
    nth = np.tile(np.arange(2*NN)+1, N)
    #print("offs: ",np.shape(offs))
    for j, n in enumerate(NNsimu):
        offs[j][n::] = np.nan
        #print(offs[j])
    offs = np.sort(offs, axis=1).flatten()
    

    
    
#ofn = "simu.dat"
NN=3
key=10

N1, sig0, sig1 = 5000, 1, 10    # Real matches
N2, dens0, dens1 = 73000, 0.2, 1.2 # Random matches
max_match_dist = 60 # arcsec


dens1 = uniform.rvs(size=N1)*(dens1-dens0)+dens0
xx = gen_random(dens1)
exit()

#-------------------------------------------------------------

sigs = uniform.rvs(size=N1)*(sig1-sig0)+sig0
offs1 = gen_real_pos_offset(N=N1, sigma=sigs)
dens1 = uniform.rvs(size=N1)*(dens1-dens0)+dens0
offs2, dens2, group, nth = gen_random_pos_offset(N=N2, dens=dens_rand/3600)

#plt.scatter(



dens_real, dens_rand = generate_random_densities(N1, N2, dens0=dens0, dens1=dens1, dens_scaling='uniform')
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


