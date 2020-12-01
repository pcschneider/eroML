from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from learn import get_props
from sklearn.neighbors import KernelDensity

fn0 = "../../ero_data/random_eFEDS.fits"
ff0 = pyfits.open(fn0)
fd0 = ff0[1].data
gi1 = np.where(fd0["NN"] == 1)[0]
gi2 = np.where(fd0["NN"] == 2)[0]
gi3 = np.where(fd0["NN"] == 3)[0]

fn = "../../ero_data/training_eFEDS_clean.fits"
ff = pyfits.open(fn)
fd = ff[1].data

gi = np.where(fd["category"]<1)[0]

#plt.hist(fd["offset_sig"][gi], bins=30, range=(0,10), density=True)
#plt.hist(fd0["offset_sig"], bins=30, range=(0,10), density=True)
#plt.hist(fd0["offset_sig"][gi1], bins=30, range=(0,10), density=True)
#plt.hist(fd0["offset_sig"][gi2], bins=30, range=(0,10), density=True)
#plt.hist(fd0["offset_sig"][gi3], bins=30, range=(0,10), density=True)
#plt.show()

#plt.scatter(fd["offset_sig"], fd["expected_rnd"])
plt.hexbin(fd0["offset_sig"], fd0["expected_rnd"], extent=(0,20,0,40), gridsize=30)
plt.scatter(fd["offset_sig"][gi], 10*4* fd["match_dist"][gi]**2 * np.pi*10**fd["log_sk"][gi]/3600, s=3, color='r')
plt.xlim(0,20)
plt.ylim(0,40)
plt.xlabel("offset_sig")
plt.ylabel("expected_rnd")

#ii = np.where((fd["offset_sig"][gi] < 4) & (fd["expected_rnd"][gi]>25))[0]
#print(gi[ii], fd["srcID"][gi[ii]])

plt.figure()
#plt.hexbin(fd0["offset_sig"], fd0["expected_rnd"], extent=(0,50,0,700))
plt.hexbin(fd0["offset_sig"], fd0["expected_rnd"],extent=(0,20,0,40), gridsize=30)
plt.xlim(0,20)
plt.ylim(0,40)
plt.title("random")

plt.figure()           
hb = plt.hexbin(fd["offset_sig"][gi], 10*4* fd["match_dist"][gi]**2 * np.pi*10**fd["log_sk"][gi]/3600, extent=(0,20,0,40), gridsize=30)



plt.xlim(0,20)
plt.ylim(0,40)
plt.title("training")

plt.show()


verts = hb.get_offsets()
for j, v in enumerate(verts):
    print(v, hb.get_array()[j])

plt.plot(hb.get_array())

plt.show()
x = hb.get_offsets()[:,0]
y = hb.get_offsets()[:,1]
plt.scatter(x, y, c=hb.get_array())
plt.show()


X = np.array([fd["offset_sig"], fd["expected_rnd"]]).T
print(np.shape(X))
gi_real = np.where(fd["category"]<1)[0]
gi_rand = np.where(fd["category"]>0)[0]
print("#real %i, #random: %i" % (len(gi_real), len(gi_rand)))
real = KernelDensity(kernel='epanechnikov', bandwidth=1).fit(X[gi_real,:])
rand = KernelDensity(kernel='epanechnikov', bandwidth=1).fit(X[gi_rand,:])

z_real = real.score_samples(X[gi_real,:])
z_rand = rand.score_samples(X[gi_rand,:])
#print(z)



#plt.scatter(X[gi_rand,0], X[gi_rand,1], color='0.7', alpha=0.5, s=1)

#plt.scatter(X[gi_real,0], X[gi_real,1], c=z_real, vmin=-3, vmax=3, s=1)


xlim = plt.xlim()
ylim = plt.ylim()

xlim = [0, 10]
ylim = [0, 100]

print(xlim, ylim)
xx = np.linspace(xlim[0], xlim[1], 100)
yy = np.linspace(ylim[0], ylim[1], 100)

YY, XX = np.meshgrid(yy, xx)
print(np.shape(XX))
ZZ =  np.array([XX.ravel(), YY.ravel()]).T
print("ZZ",np.shape(ZZ))

print(ZZ)
z_ratio = real.score_samples(ZZ) - rand.score_samples(ZZ)
z_ratio = z_ratio.reshape((100,100))
print("z_ratio: ",np.shape(z_ratio))
plt.imshow(z_ratio, origin="lower")
#print(np.shape(z_real), np.shape(z_rand), np.shape(X[gi_real]))
plt.colorbar()

plt.scatter(X[:,0], X[:,1], color='0.7', alpha=0.5, s=1)

#z_real = real.score_samples(X[gi_real,:])
#z_rand = rand.score_samples(X[gi_real,:])

#plt.scatter(X[gi_real,0], X[gi_real,1], c=z_real-z_rand, vmin=-3, vmax=3)

plt.xlabel("offset_sig")
plt.ylabel("expected_rnd")
plt.show()
