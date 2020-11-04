import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn import mixture

dd = np.genfromtxt("offs.dat", unpack=True)
#dd[1]*=100
#X = np.transpose([dd[0], dd[1]])
X = dd[0:3].T

gi_real = np.where(dd[3]==0)[0]
gi_rand = np.where(dd[3]==1)[0]

real = KernelDensity(kernel='epanechnikov', bandwidth=5).fit(X[gi_real,:])
rand = KernelDensity(kernel='epanechnikov', bandwidth=5).fit(X[gi_rand,:])

z_real = real.score_samples(X[gi_rand,:])
z_rand = rand.score_samples(X[gi_rand,:])
#print(z)

plt.scatter(X[gi_rand,0], X[gi_rand,1], color='0.7', alpha=0.5, s=1)

plt.scatter(X[gi_rand,0], X[gi_rand,1], c=z_real-z_rand, vmin=-3, vmax=3, s=1)

print(np.shape(z_real), np.shape(z_rand), np.shape(X[gi_real]))

z_real = real.score_samples(X[gi_real,:])
z_rand = rand.score_samples(X[gi_real,:])

plt.scatter(X[gi_real,0], X[gi_real,1], c=z_real-z_rand, vmin=-3, vmax=3)





plt.colorbar()
plt.show()



