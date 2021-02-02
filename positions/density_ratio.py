import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn import mixture

def one_ratio(x0, y0, x1, y1, ax=None, scaling='number'):
    if ax is None: ax=plt.gca()
    N0, N1 = len(x0), len(x1)
    print("N0: %i, N1: %i " % (N0, N1))
    if type(scaling) == str:
        if scaling.lower()=="number":
            scaling = 1
        elif scaling.lower() == 'density':
            scaling = N0/N1
    
            
    
    x = np.concatenate((x0, x1))
    y = np.concatenate((y0, y1))
    
    X0 = np.transpose([x0, y0])
    X1 = np.transpose([x1, y1])
    
    d0 = KernelDensity(kernel='epanechnikov', bandwidth=1.5).fit(X0)
    d1 = KernelDensity(kernel='epanechnikov', bandwidth=1.5).fit(X1)

    #z0 = real.score_samples(X[gi_rand,:])
    #z1 = rand.score_samples(X[gi_rand,:])
    ##print(z)

    ##plt.scatter(X[gi_rand,0], X[gi_rand,1], color='0.7', alpha=0.5, s=1)

    ##plt.scatter(X[gi_rand,0], X[gi_rand,1], c=z_real-z_rand, vmin=-3, vmax=3, s=1)

    #print(np.shape(z_real), np.shape(z_rand), np.shape(X[gi_real]))

    #z_real = real.score_samples(X[gi_real,:])
    #z_rand = rand.score_samples(X[gi_real,:])


    xlim = (np.min(x), np.max(x)) # plt.xlim()
    ylim = (np.min(y), np.max(y)) # plt.ylim()

    xlim = [0, 15]
    ylim = [0, 30]

    print(xlim, ylim)
    xx = np.linspace(xlim[0], xlim[1], 100)
    yy = np.linspace(ylim[0], ylim[1], 100)

    YY, XX = np.meshgrid(yy, xx)
    #print(np.shape(XX))
    ZZ =  np.array([XX.ravel(), YY.ravel()]).T
    #print("ZZ",np.shape(ZZ))
    z_ratio = d0.score_samples(ZZ) - d1.score_samples(ZZ)

    z_ratio = z_ratio.reshape((100,100)) - np.log10(scaling)
    print("z_ratio: ",np.shape(z_ratio))
    mm = plt.imshow(z_ratio.T, extent=(min(xx), max(xx), min(yy), max(yy)), origin="lower", aspect='auto', vmin=-3, vmax=5)
    con = plt.contour(XX, YY, z_ratio, colors='r', levels=[2, 3,4,5], alpha=1.0,
                    linestyles=['-.', '-', '--',':'])

    Ncon = len(con.collections)
    for i in range(Ncon):
        #print(dir(con.collections[i]))
        Nislands = len(con.collections[i].get_paths())
        print(i, Nislands, con.levels[i])
        for j in range(Nislands):
            p = con.collections[i].get_paths()[j]
            #print(i, j, p)
            v = np.transpose(p.vertices)
            plt.plot(v[0], v[1], alpha=0.3, lw=4, color='k')
            #print(v, np.shape(np.transpose(v)))
        print()
        
    cb = plt.colorbar(mm)
    cb.set_label("Log Density Ratio")

    #print(np.shape(z_real), np.shape(z_rand), np.shape(X[gi_real]))
    #plt.colorbar()

    #plt.scatter(X[:,0], X[:,1], color='0.7', alpha=0.5, s=2)
    #plt.scatter(X[gi_real,0], X[gi_real,1], c=z_real-z_rand, vmin=-3, vmax=3)
    #plt.scatter(X[gi_real,0], X[gi_real,1], c='k', s=3, alpha=0.5)
    #plt.scatter(X[si[0:N],0], X[si[0:N],1], color='c', marker='.', s=1)
    plt.plot([0,15],[0,1.5*15], color='k', ls=':')
    plt.plot([1,10] ,[4.3,7], color='k', ls='--', lw=2)
    #plt.scatter(X[gi_rand,0], X[gi_rand,1], color='0.7', alpha=0.5, s=1)
    plt.xlabel("RADEC_ERR (arcsec)")
    plt.ylabel("Match distance (arcsec)")
    plt.xlim(1,11)
    plt.ylim(0,26)
    return ax
    #plt.show()

sfn = "../offs.dat"
#sfn = "simu.dat"
print("Reading sfn=",sfn)
dd = np.genfromtxt(sfn, unpack=True)
#dd[1]*=100
#X = np.transpose([dd[0], dd[1]])

gi0 = np.where(dd[4]==0)[0]
gi1 = np.where(dd[4]>0)[0]
ax = one_ratio(dd[0][gi0], dd[1][gi0], dd[0][gi1], dd[1][gi1], scaling='number')
plt.show()


dd1 = np.genfromtxt("simu.dat", unpack=True)
ax = one_ratio(dd[0][gi0], dd[1][gi0], dd1[0], dd1[1], scaling='density')
plt.show()
