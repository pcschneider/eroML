import numpy as np
#from eroML.utils import Dist_model
import matplotlib.pyplot as plt
from scipy.stats import norm

def gen_real_pos_offset(N, sig=1.):
        
    rnd = np.random.rand(N)
    rndx = np.sqrt(2*sig**2 * (-np.log(1-rnd)))

    return rndx


def gen_random_pos_offset(dens=1., max_dist=70, sigma_in=None):
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
    NN = np.pi*dens*max_dist**2
    N_simu = np.random.poisson(lam=NN)    
    dens_simu = np.repeat(dens, N_simu)
    sigma_simu = np.repeat(sigma_in, N_simu)
    offs = max_dist*np.random.rand(np.sum(N_simu))**0.5
    return sigma_simu, dens_simu, offs



#dens = np.linspace(0.1, 0.3, 3)
#xx = gen_random_pos_offset(12, dens=dens)
#print(xx)
#exit()
#Nreal=int(1e7)
#Nrnd = int(1e7)
#SIG = 4
#real = gen_real_pos_offset(Nreal, sig=SIG)
#fake = gen_random_pos_offset(Nrnd, dens=1/32/np.pi)

#plt.hist(real, alpha=0.5, label="Real", bins=40, range=(0, 20))
#plt.hist(fake, alpha=0.5, label="Fake", bins=40, range=(0, 20))
#plt.legend()
#plt.show()
#exit()



#sk_array = [0.0003, 0.0001, 0.001, 0.003, 0.01, 0.03, 0.1]
Nreal=50000
Nrnd = 5000
SIG = 4
SIG = np.random.rand(Nreal)*20+1
sk = np.random.rand(Nreal)**2*0.001+0.0001
#sk = Nreal*[0.004]

sig_simu, sk_simu, rnd_offs = gen_random_pos_offset(dens=sk, sigma_in=SIG)
Nrnd = len(rnd_offs)
print("Simulating %i and %i real and random sources, respectively." % (Nreal, Nrnd))

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
    
np.savetxt("offs.dat", np.transpose([sigout, pos_off, skdens*1000, cls]))
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
