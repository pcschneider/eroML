import numpy as np
#from eroML.utils import Dist_model
import matplotlib.pyplot as plt
from scipy.stats import norm

def gen_real_pos_offset(N, sig=1.):
        
    rnd = np.random.rand(N)
    rndx = np.sqrt(2*sig**2 * (-np.log(1-rnd)))

    return rndx


def gen_random_pos_offset(dens=1., NN=3):
    """
    Parameters
    -----------
    dens : array of float
        Density (arcsec-2) of stars

    Returns
    -------
    dens_simu, offs
    """
    if type(dens) != float:
        dens = np.array(dens)
    max_dist = np.repeat(np.sqrt(3/np.pi/dens), NN)
    dd = np.repeat(dens,NN)    
    N = len(dens)
    print(N, np.shape(max_dist))
    offs = max_dist*np.random.rand(N*NN)**0.5
    print(np.shape(offs))
    return offs



def gen_random_pos_offset2(dens=1., max_dist=70, sigma_in=None):
    """
    Parameters
    -----------
    dens : array of float
        Density (arcsec-2) of stars

    Returns
    -------
    dens_simu, offs
    """
    if type(dens) != float:
        dens = np.array(dens)
    NN = 3*len(dens)
    N_simu = np.random.poisson(lam=NN)    
    dens_simu = np.repeat(dens, N_simu)
    sigma_simu = np.repeat(sigma_in, N_simu)
    offs = max_dist*np.random.rand(np.sum(N_simu))**0.5
    return sigma_simu, dens_simu, offs


#sk_array = [0.0003, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1]
Nreal=2000
Nrnd = 5000
SIG = 4
SIG = np.random.rand(Nreal)*20+1
#sk = np.random.rand(Nreal)**2*0.001+0.0001 # (mean eFEDS: 0.00014)
sk = np.random.rand(Nreal)*0.001+0.0001 # (mean eFEDS: 0.00014)
#sk = Nreal*[0.004]

#rnd_offs = gen_random_pos_offset(dens=sk)
rnd_offs = gen_random_pos_offset(dens=sk)
sk_simu = np.repeat(sk,3)
sig_simu = np.repeat(SIG,3)



Nrnd = len(rnd_offs)
print("Simulating %i and %i real and random sources, respectively." % (Nreal, Nrnd))
print(len(sk_simu))


#pos_off = np.zeros((Nreal+Nrnd)
#skdens = np.zeros((Nreal+Nrnd)*len(sk_array))
cls = np.zeros(Nreal+Nrnd)
cls[0:Nreal] = 0
cls[Nreal::] = 1

tmp0 = gen_real_pos_offset(Nreal, sig=SIG)
#tmp1 = gen_random_pos_offset(Nrnd, dens=sk_simu)
pos_off = np.concatenate((tmp0, rnd_offs))
skdens = np.concatenate((sk, sk_simu))
sigout = np.concatenate((SIG, sig_simu))   
    
np.savetxt("offs.dat", np.transpose([sigout, pos_off, skdens*10000, cls]))
exit()
print(offs)
plt.hist(offs)

dm = Dist_model()

#"fraction","sig","dens","N", "err_scaling"]


dm["fraction"] = 1.0
dm["N"] = 1
dm["sig"] = SIG
dm["err_scaling"] = 1.

x = np.linspace(0, 5*SIG, 1000)
y = dm.evaluate(x)
plt.plot(x,y*N/2.4)
plt.show()
#plt.plot(x,np.cumsum(y)/np.sum(y))
#plt.show()
