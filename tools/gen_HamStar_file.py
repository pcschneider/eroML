import numpy as np
from astropy.io import fits as pyfits
from astropy.coordinates import SkyCoord
import astropy.units as u
import copy
from eroML.positions import calc_sigma_from_RADEC_ERR

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


fn = "major_eFEDS_classified.fits"
ofn = "eFEDS_HamStar.fits"
overwrite=True
verbose=1
ff = pyfits.open(fn)
ext = 1
fd = ff[ext].data
column_formats = {col.name:col.format for col in ff[ext].columns}
cols = []


# Get names
cc0 = SkyCoord(fd["RA"], fd["DEC"], unit=(u.degree, u.degree))
names = ero_names(cc0, prefix="eFEDS")
col = pyfits.Column(name="ero_NAME" , array=names, format="32A")    
cols.append(col)

mapper = {"match_dist":"match_dist",\
                "Gmag":"Gmag","plx":"plx","plx_error":"plx_error",\
                 "gaia_RA":"gaia_RA", "gaia_Dec":"gaia_Dec",\
                 "ero_RA":"RA", "ero_Dec":"Dec",\
             "sigma_r":"RADEC_sigma","gaia_ID":"srcID_NN","ero_ID":"original_srcID",\
                 "BP_RP":"bp_rp", "plx":"plx", "plx_error":"plx_error"}

for colname in mapper.keys():
    print(colname)
    arr = fd[mapper[colname]]
    col = pyfits.Column(name=colname , array=arr, format=column_formats[mapper[colname]])    
    cols.append(col)


col = pyfits.Column(name="Fx", array=10**(fd["Fx"]-13), format="D")    
cols.append(col)


#sigma_r = calc_sigma_from_RADEC_ERR(fd["RADEC_ERR"])
#col = pyfits.Column(name="sigma_r", array=sigma_r, format="D")    
#cols.append(col)


# Dummy columns:
for colname in ["rate","rate_error", "det_likeli", "Fx_error", "pm_RA", "pm_Dec"]:
    col = pyfits.Column(name=colname , array=np.zeros(len(fd["srcID"])), format="D")    
    cols.append(col)

col = pyfits.Column(name="p_stellar", array=1-fd["category"], format="D")    
cols.append(col)

col = pyfits.Column(name="p_ij", array=np.ones(len(fd["srcID"])), format="D")    
cols.append(col)
 
print(cols)
    
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
ff.close()            
#return hdul
    #col = pyfits.Column(name="sigma_r" , array=arr, format=column_formats[colname])    


