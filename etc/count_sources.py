from astropy.io import fits as pyfits
import numpy as np
from eroML.tile import file4


fn = file4("X_filename_hp", cconfig="eFEDS_EDR3.ini") # "../ero_data/Gaia_merged_eFEDS_EDR3.fits"
print("Reding ",fn)
ff = pyfits.open(fn)
fd = ff[1].data 
gi = np.where((fd["DET_LIKE_0"] > 6) & (fd["EXT_LIKE"]<=6))[0]
print("Number of suitable eROSITA sources: ",len(gi))
ff.close()

fn = "../ero_data/Gaia_merged_eFEDS_EDR3.fits"
print("Reding ",fn)
ff = pyfits.open(fn)
fd = ff[1].data 
gi = np.where(fd["eligible_Gaia"] ==1)[0]
print("Number of eligible Gaia sources: ",len(gi))
ff.close()





