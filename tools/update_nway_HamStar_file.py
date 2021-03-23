import numpy as np
from astropy.io import fits as pyfits
from astropy.coordinates import SkyCoord
import astropy.units as u
from eroML.positions import calc_sigma_from_RADEC_ERR


fn = "../ero_data/nway.fits"
ofn = "../ero_data/nway2.fits"
overwrite=True
verbose=1
ff = pyfits.open(fn)
ext = 1
fd = ff[ext].data
column_formats = {col.name:col.format for col in ff[ext].columns}
cols = []

si = np.where((fd["Gaia_ID"] != "0") & (fd["type"] == "S"))[0]
print(len(si))

def ero_names(coords, prefix):
    names=[]
    for c in coords:
        #print(c.to_string("hmsdms"))
        ra = str("%02i%02i%05.2f " % (c.ra.hms.h, c.ra.hms.m, c.ra.hms.s))
        
        dec = str("%02i%02i%05.2f" % (abs(c.dec.dms.d), abs(c.dec.dms.m),  abs(c.dec.dms.s)))
        sign="+" if c.dec.degree>=0 else "-"
        
        name = prefix+" J"+ra[0:-2]+ sign+dec[0:-3]
        #print(c, name)
        names.append(name)
    return names

# Get names
print("Generating ero-names...")
cc0 = SkyCoord(fd["ero_RA"], fd["ero_DEC"], unit=(u.degree, u.degree))
names = ero_names(cc0[si], prefix="eFEDS")
col = pyfits.Column(name="ero_NAME" , array=names, format="32A")    
cols.append(col)
print("done.")

for c in column_formats.keys():
    print(c)
    arr = fd[c]
    
    if c == "ero_ID":
        arr = np.array([str("ML%05i" % int(d)) for d in fd[c]])
        col = pyfits.Column(name=c, array=arr[si], format=column_formats[c])    
    elif c == "Gaia_ID":    
        gi = np.where(arr == "0")[0]
        print("len",len(gi))
        arr[gi] = ""
        col = pyfits.Column(name="gaia_ID", array=arr[si], format=column_formats[c])
    elif c == "Gaia_RA":    
        col = pyfits.Column(name="gaia_RA", array=arr[si], format=column_formats[c])      
    elif c == "Gaia_Dec":    
        col = pyfits.Column(name="gaia_Dec", array=arr[si], format=column_formats[c])      
    elif c == "plx_err":    
        col = pyfits.Column(name="plx_error", array=arr[si], format=column_formats[c])
    else:
        col = pyfits.Column(name=c, array=arr[si], format=column_formats[c])    
        
    cols.append(col)


sigma_r = calc_sigma_from_RADEC_ERR(fd["RADEC_ERR"])
col = pyfits.Column(name="sigma_r", array=sigma_r[si], format="D")    
cols.append(col)


col = pyfits.Column(name="p_stellar", array=fd["p_any"][si], format="D")    
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
    
        

