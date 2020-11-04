from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob

directory = "../ero_data/"
glob_str = directory+"random_ID2_nside32_*_small.fits"
print(glob_str)


fnames = glob.glob(glob_str)


dens = []
number = []
number2 = []

for j, fn in enumerate(fnames):#[0:2]:
    print(fn, "(",j+1,"/",len(fnames),")")
    ff = pyfits.open(fn)
    sk = np.mean(ff[1].data["eligible_sky_density"])
    gi = np.where(ff[1].data["offset_sig"] < 2)[0]
    gi2 = np.where(ff[1].data["match_dist"] < 10)[0]
    dens.append(sk)
    number.append(len(gi) / len(ff[1].data["offset_sig"]))
    number2.append(len(gi2) / len(ff[1].data["offset_sig"]))
    N = len(ff[1].data["offset_sig"])
    
    plt.hist(ff[1].data["match_dist"], range=(0, 60), bins=30)
    
    x = np.linspace(0,60, 100)
    N  = np.mean(sk)*2 * np.pi* 9
    eta = (N-1)/np.pi/9/3600
    
    #plt.
    
    #eta = (-1)/3600
    
    plt.title(str(np.mean(eta)))
    
    y = x * np.pi * x * eta * np.exp(-np.pi*eta*x**2)
    plt.plot(x,y * N*2)
    plt.show()
    #print(len(gi))
    #print(sk)
    #print()

np.savetxt("sky.dat", np.transpose([dens, number, number2]))


