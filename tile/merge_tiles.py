from astropy.io import fits as pyfits
import glob
import numpy as np
import copy

ofn = "test_merged.fits"
overwrite=True
glob_str = "../../ero_data/*rID1_training.fits"
fnames = glob.glob(glob_str)[0:10]
print("Number of matching files: ",len(fnames))

# Count required entries
N = 0
for fn in fnames:
    ff = pyfits.open(fn)
    tmp = len(ff[1].data["srcID"])
    N+=tmp
    print(fn, " -> ", tmp, "(",N,")")
    
ff = pyfits.open(fnames[0])
cols = ff[1].columns
col_dct = {}
for c in cols:
    col_dct[c.name] = {"array":N*[0], "col":c}
    
i0, i1 = None, None

for fn in fnames:
    #print(fn)    
    ff = pyfits.open(fn)
    M = len(ff[1].data["srcID"])
    if i0 == None:
        i0, i1 = 0, M
    else:
        i0, i1 = i1, i1+M
        
    for c in cols:
        #col_dct[c.name]["array"][i0:i1] = copy.copy(ff[1].data[c.name])
        col_dct[c.name]["array"][i0:i1] = ff[1].data[c.name]
        col_dct[c.name]["col"] = c

#print(col_dct["srcID"])        
    
out_cols = []    
for c in cols:
    #print(c)
    if c.name == "Gaia_Quality":
        col = pyfits.Column(name=c.name, array=col_dct[c.name]["array"]=="True", format="J")
    else:
        try:
            col = pyfits.Column(name=c.name, array=col_dct[c.name]["array"], format=col_dct[c.name]["col"].format)
        except:
            print("XXX",c.name)
            col = pyfits.Column(name=c.name, array=col_dct[c.name]["array"]=="True", format=col_dct[c.name]["col"].format)
    out_cols.append(col)
        
hdu = pyfits.PrimaryHDU()    
cc = pyfits.ColDefs(out_cols)
xx = pyfits.BinTableHDU.from_columns(cc)
hdul = pyfits.HDUList([hdu, xx])

if ofn is not None:
    hdul.writeto(ofn, overwrite=overwrite)        
    print("Written ",N," objects with ",len(out_cols)," properties to ",ofn)
#return hdul
    
