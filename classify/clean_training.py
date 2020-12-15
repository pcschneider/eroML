from astropy.io import fits as pyfits
import numpy as np
#../ero_data/random_eFEDS.fits
fn = "../../ero_data/training_eFEDS.fits"
ofn = "../../ero_data/training_eFEDS_clean.fits"
ff = pyfits.open(fn)
fd = ff[1].data

ii = np.where((fd["srcID"] != "ML26649") & (fd["srcID"] != "ML23615"))[0]
print(len(fd["srcID"]), len(ii))
cols = []

for c in ff[1].columns:
    cols.append(
            pyfits.Column(name=c.name, format=c.format, unit=c.unit, array=fd[c.name][ii])
            )
              
new_cols = pyfits.ColDefs(cols)
hdu = pyfits.BinTableHDU.from_columns(new_cols)
hdu.writeto(ofn, overwrite=True)    
