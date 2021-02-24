import glob
import numpy as np
from eroML.utils import sky_dens4coordinates, sky_density
from eroML.tile import file4
from eroML.config import *
import time
from astropy.io import fits as pyfits
#sky_dens4coordinates


conf_fn = "eRASS1_datasets_only.ini"
conf_fn = "eFEDS_EDR3.ini"
#tmp = read_config(conf_fn)

#print(config)

oo = open("dens_log.dat","w")
oo.write("# fn N dens\n")

fnames = file4("gaia_tiles", cconfig=conf_fn)
print("#Filenames: ",len(fnames))
print("    From ",fnames[0], " to ",fnames[-1])
for i, fn in enumerate(fnames):
    print("===============", time.asctime(), " ",i+1, "/",len(fnames), " - Working on: ",fn)
    sky_density(fn)
    print(" >>> Checking")
    ff = pyfits.open(fn)
    nm = np.nanmean(ff[1].data["eligible_sky_density"])
    gi = np.where(ff[1].data["eligible_Gaia"] == 1)[0]
    ostr = "%s %i %f\n" % (fn, len(gi), nm)
    oo.write(ostr)
    oo.flush()
    ff.close()
    print(ostr)
oo.close()    
