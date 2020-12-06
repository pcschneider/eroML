import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn import mixture

dd = np.genfromtxt("../offs.dat", unpack=True)
#dd[1]*=100
#X = np.transpose([dd[0], dd[1]])
X = dd[0:2].T

gi_real = np.where(dd[4]==0)[0]
gi_rand = np.where(dd[4]==1)[0]

real = KernelDensity(kernel='epanechnikov', bandwidth=1.5).fit(X[gi_real,:])
rand = KernelDensity(kernel='epanechnikov', bandwidth=1.5).fit(X[gi_rand,:])

z_real = real.score_samples(X[gi_rand,:])
z_rand = rand.score_samples(X[gi_rand,:])
#print(z)

#plt.scatter(X[gi_rand,0], X[gi_rand,1], color='0.7', alpha=0.5, s=1)

#plt.scatter(X[gi_rand,0], X[gi_rand,1], c=z_real-z_rand, vmin=-3, vmax=3, s=1)

print(np.shape(z_real), np.shape(z_rand), np.shape(X[gi_real]))

z_real = real.score_samples(X[gi_real,:])
z_rand = rand.score_samples(X[gi_real,:])


xlim = plt.xlim()
ylim = plt.ylim()

xlim = [0, 15]
ylim = [0, 30]

print(xlim, ylim)
xx = np.linspace(xlim[0], xlim[1], 100)
yy = np.linspace(ylim[0], ylim[1], 100)

YY, XX = np.meshgrid(yy, xx)
#print(np.shape(XX))
ZZ =  np.array([XX.ravel(), YY.ravel()]).T
#print("ZZ",np.shape(ZZ))
z_ratio = real.score_samples(ZZ) - rand.score_samples(ZZ)




score = real.score_samples(X) - rand.score_samples(X)
print(score)#si = np.argsort(score)
si = np.argsort(score)[::-1]
N = len(gi_real)
print("Nreal: ",N)
print("recovered: ",N - np.sum(dd[4][si[0:N]]))
print("ratio range: ",score[si[0]], score[si[N-1]])

#print(ZZ)
z_ratio = z_ratio.reshape((100,100))
print("z_ratio: ",np.shape(z_ratio))
mm = plt.imshow(z_ratio.T, extent=(min(xx), max(xx), min(yy), max(yy)), origin="lower", aspect='auto', vmin=-3, vmax=5)
con = plt.contour(XX, YY, z_ratio, colors='r', levels=[2, 3,4,5], alpha=1.0,
                linestyles=['-.', '-', '--',':'])

cb = plt.colorbar(mm)
cb.set_label("Log Density Ratio")

#print(np.shape(z_real), np.shape(z_rand), np.shape(X[gi_real]))
#plt.colorbar()

plt.scatter(X[:,0], X[:,1], color='0.7', alpha=0.5, s=2)
#plt.scatter(X[gi_real,0], X[gi_real,1], c=z_real-z_rand, vmin=-3, vmax=3)
plt.scatter(X[gi_real,0], X[gi_real,1], c='k', s=3, alpha=0.5)
plt.scatter(X[si[0:N],0], X[si[0:N],1], color='c', marker='.', s=1)
plt.plot([0,15],[0,1.5*15], color='k', ls=':')
plt.plot([1,10] ,[4.3,7], color='k', ls='--', lw=2)
#plt.scatter(X[gi_rand,0], X[gi_rand,1], color='0.7', alpha=0.5, s=1)
plt.xlabel("RADEC_ERR (arcsec)")
plt.ylabel("Match distance (arcsec)")
plt.xlim(1,11)
plt.ylim(0,26)
plt.show()

# Nreal:  2000
# recovered:  1587.0
# stars as stars: 1587
# random as stars: 413
# stars as random: 413


