import matplotlib.pyplot as plt
import numpy as np
from scipy.special import factorial
from eroML.positions import gen_random_pos_offset

np.random.seed(1)

NN=10
dens = NN/np.pi/(np.arange(5)+1)**2
dens = NN/(np.random.rand(100000)+1)
#print(dens)
dens = 500000*[1]

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
