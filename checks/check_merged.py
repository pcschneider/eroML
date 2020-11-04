from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import PyAstronomy.funcFit as fuf
from scipy.stats import invgauss, chi2, lognorm, recipinvgauss,skewnorm,norminvgauss

fn = "../ero_data/merged_major.fits"
ff = pyfits.open(fn)
dd = ff[1].data

N = len(np.where(dd["NN"] == 1)[0])

for i in [1,2,3]:
    gi = np.where(ff[1].data["NN"] == i)[0]
    plt.hist(ff[1].data["match_dist"][gi], bins=61, range=(0, 40), label=str(i), alpha=0.3)
    
plt.legend()    
plt.show()


#print(dd["eligible"])
gi0 = np.where( (dd["eligible"]==0) & (dd["NN"] == 1) )[0]
gi1 = np.where( (dd["eligible"][gi0]==1) & (dd["NN"][gi0] == 2) )[0]
gi1 = gi0[gi1]
print(N, len(gi0), len(gi1))
