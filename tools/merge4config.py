from eroML.tile import file4, calc_hpix
from eroML.config import *
from eroML.ensemble import from_fits, to_fits
import numpy as np

conf_fn = "eRASS1_datasets_only.ini"
conf_fn = "eFEDS_EDR3.ini"

def prep_one_tile(e):
    ra, dec = e.to_array("RA", array_type="array"), e.to_array("Dec", array_type="array")
    hpx = calc_hpix(ra, dec, nside=nside)
    this_hpx = int(np.median(hpx))
    gi = np.where(hpx==this_hpx)[0]
    print(len(gi), len(hpx))
    e.keep(np.array(e.srcIDs())[gi])

#tmp = read_config(conf_fn)
cc = read_config(conf_fn)
nside = int(config["Healpix"]["nside"])

fnames = file4("gaia_tiles", cconfig=conf_fn)
N = len(fnames)
print("Processing %i tiles (nside=%i)." % (N, nside))
e0 = from_fits(fnames[0])
prep_one_tile(e0)

for i, fn in enumerate(fnames[1:]):
    e = from_fits(fn)
    
    ee = len(e)
    prep_one_tile(e)
    
    e0.append(e, duplicates='ignore')
    print(fn,": ", len(e0), " (",i+1,"/",N,")")  
    print()

to_fits(e0, ofn="x.fits", overwrite=True)


import matplotlib.pyplot as plt
from eroML.ensemble import from_fits, to_fits
import numpy as np

e0 = from_fits("x.fits")
ra, dec = e0.to_array("RA", array_type="array"), e0.to_array("Dec", array_type="array")
skd = e0.to_array("eligible_sky_density", array_type="array")
N = len(ra)
print(N)

idx = np.random.choice(np.arange(N), N)
print(idx)
sc = plt.scatter(ra[idx], dec[idx], label="merged", s=10, c=skd[idx], alpha=0.5)
plt.xlabel("RA")
plt.ylabel("Dec")
cb = plt.colorbar(sc)
cb.set_label(r"#stars/arcmin$^2$")
#plt.legend()
plt.show()

