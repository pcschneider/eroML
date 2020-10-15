import numpy as np
from eroML.utils import Dist_model
import matplotlib.pyplot as plt
from scipy.stats import norm

def gen_real_pos_offset(N, sig=1.):
        
    rnd = np.random.rand(N)
    rndx = np.sqrt(2*sig**2 * (-np.log(1-rnd)))

    return rndx

def gen_random_pos_offset(N, dens=1.):

    rnd = np.random.rand(N)
    rndx = np.sqrt(-np.log(1-rnd)/(np.pi*dens))
    return rndx


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
Nreal=5000
Nrnd = 5000
SIG = 4

#pos_off = np.zeros((Nreal+Nrnd)
#skdens = np.zeros((Nreal+Nrnd)*len(sk_array))
cls = np.zeros(Nreal+Nrnd)
cls[0:Nreal] = 0
cls[Nreal::] = 1

sk = np.random.rand(Nreal+Nrnd)**2*0.01+0.0001
tmp0 = gen_real_pos_offset(Nreal, sig=SIG)
tmp1 = gen_random_pos_offset(Nrnd, dens=sk[Nreal:Nreal+Nrnd])
pos_off = np.concatenate((tmp0, tmp1))
skdens = sk


#for j, sk in enumerate(sk_array):
    #SIG = 4
    
    #offs = gen_real_pos_offset(Nreal, sig=SIG)
    #i0=j*(Nreal+Nrnd)
    #i1=j*(Nreal+Nrnd)+Nreal
    #print(i0, i1)
    #pos_off[i0:i1] = offs
    #skdens[i0:i1] = sk
    #cls[i0:i1] = 0
    
    #gg = norm.rvs(size=Nrnd//20, loc=sk, scale=sk)
    #for g in gg:
        #offs = gen_random_pos_offset(Nrnd//len(gg), dens=abs(sk))
        #i0=i1
        #i1+=len(offs)
        #print(i0, i1, g)
        #pos_off[i0:i1] = offs
        #skdens[i0:i1] = abs(g)
        #cls[i0:i1] = 1    
    #print()
    
    
    
np.savetxt("offs.dat", np.transpose([pos_off, skdens, cls]))
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
