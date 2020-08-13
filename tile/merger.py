from astropy.io import fits as pyfits
import glob
import numpy as np
import copy
#import eroML.config as conf
import sys

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
    
def merge_matching(glob_str, ofn=None, overwrite=True, verbose=1):
    """
    Merge files matching `glob_str` into one fits-file
    """
    fnames = glob.glob(glob_str)#[0:10]
    if verbose>0: print("tile.merger::merge_matching - Number of matching files: ",len(fnames))

    # Count required entries
    N = 0    
    for fn in progressbar(fnames, "tile.merger::merge_matching - Counting number of objects "):
        ff = pyfits.open(fn)
        tmp = len(ff[1].data["srcID"])
        N+=tmp
        if verbose>1: print(k+1,"/",len(fnames),"   ", fn, " -> ", tmp, "(",N,")")
        ff.close()
    
    
    # Extract relevant columns, prepare final array
    ff = pyfits.open(fnames[0])
    cols = ff[1].columns
    col_dct = {}
    for c in cols:
        col_dct[c.name] = {"array":N*[0], "col":c}
    ff.close()
    
    if verbose>0: print("tile.merger::merge_matching - Number of rows: ",N)

    i0, i1 = None, None
    #for k, fn in enumerate(fnames):
    for fn in progressbar(fnames, "tile.merger::merge_matching - Merging objects into one file "):        
        ff = pyfits.open(fn)
        #if verbose>2: print(k+1,"/",len(fnames),"  Adding sources from ", fn)
        M = len(ff[1].data["srcID"])
        if i0 == None:
            i0, i1 = 0, M
        else:
            i0, i1 = i1, i1+M
            
        for c in cols:
            col_dct[c.name]["array"][i0:i1] = ff[1].data[c.name]
            col_dct[c.name]["col"] = c
            
        ff.close()
        
    # Creating new fits-file and some streamlining of formats for output
    out_cols = []    
    for c in cols:
        #print(c)
        if c.name == "Gaia_quality" or c.name == "Gaia_Quality":
            ar = np.array(col_dct[c.name]["array"])=="True"
            #print(ar, col_dct[c.name]["array"])
            col = pyfits.Column(name=c.name, array=ar, format="J")
        else:
            try:
                col = pyfits.Column(name=c.name, array=col_dct[c.name]["array"], format=col_dct[c.name]["col"].format)
            except:
                #print("XXX",c.name)
                col = pyfits.Column(name=c.name, array=col_dct[c.name]["array"]=="True", format=col_dct[c.name]["col"].format)
        out_cols.append(col)
    #print(out_cols)        
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(out_cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])

    if ofn is not None:
        hdul.writeto(ofn, overwrite=overwrite)        
        print("Written ",N," objects with ",len(out_cols)," properties to ",ofn)
    
        
 
    
    
#ofn = "test_merged.fits"
#overwrite=True
#glob_str = "../../ero_data/*rID1_training.fits"
#fnames = glob.glob(glob_str)#[0:10]

#ofn = "random_merged.fits"
#overwrite=True
#glob_str = "../../ero_data/*rID1_random.fits"
#fnames = glob.glob(glob_str)[0:10]

