import matplotlib.pyplot as plt
import numpy as np
from scipy.special import factorial

np.random.seed(1)

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
    N = len(dens)
    M = NN
    NN*=3
    max_dist = np.repeat(np.sqrt(NN/np.pi/dens), NN*2).reshape((N,NN*2))
    #print(np.shape(max_dist), max_dist)
    dd = np.repeat(dens,NN*2).reshape((N,NN*2))
    gg = np.repeat(np.arange(len(dens)),NN*2).reshape((N,NN*2))
    print(np.shape(gg), np.shape(dd))
    #print(N, np.shape(max_dist), max(max_dist))
    NNsimu = np.random.poisson(lam=np.array(len(dens)*[NN]))
    print(np.shape(NNsimu), NNsimu)
    offs = max_dist*np.random.rand(N,NN*2)**0.5
    nth = np.tile(np.arange(2*NN)+1, N)
    print("offs: ",np.shape(offs))
    for j, n in enumerate(NNsimu):
        offs[j][n::] = np.nan
        #print(offs[j])
    offs = np.sort(offs, axis=1).flatten()
    
    step = 2*NN
    for i in range(step-M):
        j = M+i
        
        offs[j::step] = np.nan
        
    plt.hist(offs, range=(0, 2.3), bins=30)
    #print(offs[0,:])
    plt.show()
    #print(np.shape(offs))
    #nth = np.argsort(offs, axis=1)
    #offs = offs[nth]
    #nth = np.tile(np.arange(NN),N).flatten()[np.argsort(offs).flatten()]+1
    #print(nth[0,:])
    #print(offs[0,nth[0,:]-1])
    #offs = offs.flatten()
    #gi = np.where(nth.flatten() == 1)
    #print(offs)
    #print(nth)
    #print(offs[nth])
    #print(np.shape(offs))
    oo = offs.flatten()
    gi = np.where(np.isfinite(oo))[0]
    print(len(gi), len(dens)*NN)
    return np.array([offs.flatten()[gi], dd.flatten()[gi], gg.flatten()[gi], nth.flatten()[gi]])


NN=3
dens = NN/np.pi/(np.arange(5)+1)**2
dens = NN/(np.random.rand(100000)+1)
#print(dens)
#dens = 500000*[1]

offs = gen_random_pos_offset(dens=dens, NN=NN)
print(np.shape(offs))
#gg = np.where((offs[0] > 1.5) & (offs[1]==1))[0]
#print(gg)
#print(offs[3][gg], gg[0])
#cc = np.where(offs[3] == offs[3][gg[0]])[0]
##print(offs[0][cc], offs[1][cc], offs[3][cc])
print(offs[0][0::NN], offs[2][0::NN])
#exit()
#plt.hist(offs[0], range=(0, 1.1), bins=30)
px = np.linspace(0,2.3,100)

dens0 = 0.45
dens0 = np.mean(dens)
oo = np.zeros(30)
for i in [9,8,7,6,5,4,3,2,1,0]:
    nn = i+1
    gi = np.where(offs[3] == i+1)[0]
    print(i, "len",len(gi))
    bns = np.histogram(offs[0][gi], range=(0, 2.3), bins=30, density=True)
    bnx = (bns[1][1:] + bns[1][0:-1])/2
    #print(gi, np.shape(gi))
    oo+=bns[0]
    plt.bar(bnx, height=bns[0], width=0.05, label=str(i))
    plt.plot(px, 2*(np.pi*dens0)**nn * px**(2*nn-1) / factorial(nn-1) * np.exp(-np.pi*dens0*px**2))
    #oo+=
    print(i, np.mean(offs[0][i::NN]))
    print()
plt.legend()
plt.show()

plt.plot(bnx, oo, ls='steps-mid')
plt.show()
