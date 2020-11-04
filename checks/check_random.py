from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

fn = "../ero_data/merged_random_test120.fits"
ff = pyfits.open(fn)
dd = ff[1].data

for i in [1,2,3]:
    gi = np.where(dd["NN"] == i)[0]
    print(i, len(gi))
    
    
for i in range(10):
    gi = np.where(dd["N_random"] == i)[0]
    print(i, len(gi))
    
   

fn = "../ero_data/merged_test120.fits"
ff1 = pyfits.open(fn)
dd2 = ff1[1].data

abs_dist = 3
 
plt.hist(dd["offset_sig"], range=(0.1,20), bins=100, label="Random")
gi = np.where( (dd["offset_sig"] < 1) & (dd["match_dist"] < abs_dist) )[0]
print(len(gi))

   
plt.hist(dd2["offset_sig"], range=(0.1,20), bins=100, label="Real")

gi = np.where( (dd["offset_sig"] < 1) & (dd["match_dist"] < abs_dist) )[0]
plt.hist(dd["offset_sig"][gi], color='b', alpha=0.5, label="good", range=(0.1,20), bins=100)


gi = np.where( (dd2["offset_sig"] < 1) & (dd2["match_dist"] < abs_dist) )[0]
#gi = np.where( dd2["offset_sig"] < 1)[0]

print(len(gi))


plt.hist(dd2["offset_sig"][gi], color='r', alpha=0.5, label="real good", range=(0.1,20), bins=100)


plt.legend()



#plt.yscale("log")
#plt.xscale("log")


plt.show()
