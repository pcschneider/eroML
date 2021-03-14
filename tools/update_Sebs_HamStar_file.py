import numpy as np
from astropy.io import fits as pyfits

fn = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"
ofn = "../ero_data/efeds_c001_V3_main_HamStar_internal2.fits"
overwrite=True
verbose=1
ff = pyfits.open(fn)
ext = 1
fd = ff[ext].data
column_formats = {col.name:col.format for col in ff[ext].columns}
cols = []

for c in column_formats.keys():
    print(c)
    arr = fd[c]
    if c == "ero_ID":
        print("Updating")
        arr = np.array([str("ML%05i" % int(d)) for d in fd[c]])
    col = pyfits.Column(name=c, array=arr, format=column_formats[c])    
    cols.append(col)

      
hdu = pyfits.PrimaryHDU()    
cc = pyfits.ColDefs(cols)    
xx = pyfits.BinTableHDU.from_columns(cc)
dst = len(xx.data[cc[0].name])
hdul = pyfits.HDUList([hdu, xx])
print("XXXXXX")
if ofn is not None:
    hdul.writeto(ofn, overwrite=overwrite)        
    if verbose>0:
        print("datasets::prep_classify - Written ",dst," objects with ",len(cols)," properties to ",ofn)
    
        

