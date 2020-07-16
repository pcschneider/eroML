from .astro_object import Astro_Object
from .astro_ensemble import Ensemble
from astropy.io import fits as pyfits
import numpy as np

def from_fits(fn, mapper={}, verbose=1, extension=1, maxN=None):
    """
    Generate Astro_Ensemble from fits-file
    """

    col_mapper = lambda x:mapper[x] if x in mapper else x

    if verbose>0: print("ensemble.tools::from_fits - Reading ",fn)
    ff = pyfits.open(fn)    
    cols = ff[extension].columns
    if verbose>5:
        print("ensemble.tools::from_fits - ", fn, " contains ",cols)
    
    if maxN is None:
        maxN = len(ff[extension].data[cols[0].name])
        
    # Generate np.recarray for array    
    col_data = []
    names=[]
    for col in cols:
        col_data.append(ff[extension].data[col.name][0:maxN])
        if verbose>6: print("ensemble.tools::from_fits - Using col: ",col.name, " with format: ",col.format)
        names.append(col_mapper(col.name))
            
    names =     ",".join(names)    
    a = np.core.records.fromarrays(col_data,names=names)
    if verbose>1: print("ensemble.tools::from_fits - Read ",np.shape(a), " entries with ",len(names.split(","))," properties from ",fn)
    
    e = Ensemble().from_array(a)
    if verbose>0:
        print("ensemble.tools::from_fits - Generated Ensemble with",np.shape(a), " entries with ",len(names.split(","))," properties")
    return e
    
def to_fits(ensemble, ofn, overwrite=False, verbose=1, mapper={}, maxN=None):
    """
    Write Ensemble-data to fits-file
    """
    fmt_mapper = {"i":"I","u":"I","U":"32A","f":"D"}
    col_mapper = lambda x:mapper[x] if x in mapper else x

    if verbose>1:
        print("ensemble.tools::to_fits - Using cols: ",ensemble.known_cols)
    array = ensemble.array(colnames=ensemble.known_cols)
    if maxN is None:
        maxN = len(array)
    else:
        maxN = maxN if maxN<len(array) else len(array)
        
    cols = []
    for c in array.dtype.names:
        
        fmt = fmt_mapper[array[c].dtype.kind]
        if verbose>5: print("ensemble.tools::to_fits - Using column ",c, "with dtype=",array[c].dtype, " and fits-format=",fmt)
        col = pyfits.Column(name=col_mapper(c), array=array[c][0:maxN], format=fmt)
        cols.append(col)
        
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    hdul.writeto(ofn, overwrite=overwrite)        
    if verbose>0:
        print("ensemble.tools::to_fits - Written ",maxN," objects with ",len(cols)," properties to ",ofn)
