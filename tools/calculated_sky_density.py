import healpy as hp
from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import sys

nside=32
N = hp.nside2npix(nside)
directory = "../../ero_data/"

glob_str = directory+"Gaia_ID2_nside"+str(nside)+"_*.fits"
p = re.compile(directory+"Gaia_(\S*)_nside"+str(nside)+"_(\d+).fits")

#glob_str = directory+"Gaia_test30_nside"+str(nside)+"_*.fits"
#p = re.compile(directory+"Gaia_test30_nside"+str(nside)+"_(\d+).fits")



def progressbar(it, prefix="", size=10, file=sys.stdout):
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()        
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()

fnames = glob.glob(glob_str)
#fnames = ["../../ero_data/Gaia_ID2_nside32_1543.fits"]

print("Number of files: ",len(fnames), " (",glob_str,")")
cnt = np.ones(N)*0.01
for fn in progressbar(fnames, prefix="Scanning ", size=20):#[0:10]:
    ff = pyfits.open(fn)
    try:
        n = ff[1].header["NAXIS2"]# len(ff[1].data["srcID"])
        n = np.nanmean(ff[1].data["eligible_sky_density"])
    except Exception as ee:
        print("Cannot read ",fn)
        print(ee)
        continue
    mm = p.match(fn)
    if not mm: 
        print("No index found in ", fn)
        continue
    #print(mm)
    found = mm.group(2)
    #print(found)
    idx = int(found)
    cnt[idx] = n
    ff.close()
    
print("Mean sky density (1/arcmin^2): %5.2f" % np.mean(cnt)    )

hp.mollview(cnt, norm="log", nest=True)
plt.show()
    #idx = 
