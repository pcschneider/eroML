"""
Collection of Tools for handling Ensembles.

Example
-------
>>> e = from_fits("../eFEDS/SrcCat_V2T.fits", maxN=100, verbose=0)
>>> e["ML00001"]["srcID"].strip()
'ML00001'
"""

from .astro_object import Astro_Object
from .astro_ensemble import Ensemble
from astropy.io import fits as pyfits
import numpy as np
import functools

gverbose=0

def multi_fits_support(n, verbose=gverbose):
    """
    Decorator to support that the `n`-first arguments can be filenames (strings) instead of ensembles.
    
    There are two modes:
      - The number of arguments is at least `n`: The `n`-th argument will be assumed to be the output fits-file. The number or "real" arguments is `n`-1 so that the `n`-th argument will be removed from the argument list.
      - The number of arguments is `n`-1: Only the first `n`-1 arguments will be inspected for filenames. No fits-file will be written
    
    Returns
    -------
    Original result of `func` : N/A
    """
    def multi_fits_support_func(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            args_tmp = list(args)
            if len(args)<n-1:
                raise TypeError("ensemle.tools::fits_support: Expecting at least n-1 arguments, but too few provided.")
            ens = []
            for i in range(n-1):
                if type(args[i]) == type("xxx"):
                    if verbose>0: print("ensemle.tools::multi_fits_support - Assuming ", args[0], " is a fits-filename.")
                    if "mapper" in kwargs:
                        e = from_fits(args[i], mapper=kwargs["mapper"])
                    else:
                        e = from_fits(args[i])
                else:
                    e = args[i]
                args_tmp[i] = e
            if len(args)>=n:
                args_tmp.pop(n-1)
            r = func(*tuple(args_tmp), **kwargs)
            #print("YYYYY", type(r), len(r))

            if len(args)>=n and type(args[n-1]) == type("xxx"):
                if verbose>0: print("ensemle.tools::multi_fits_support - Assuming ", args[n-1], " is a fits-filename.")
                to_fits(r, ofn=args[n-1], overwrite=True)
            return r    
        return wrapper
    return multi_fits_support_func

def fits_support(func, verbose=gverbose):
    """
    Decorator to support that the first arguments can be a filename (string) instead of an ensemble.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if len(args)>0:
            x = args[0]
        else:
            raise TypeError("ensemle.tools::fits_support: Expecting at least one argument, but none provided.")
        if type(x) == type("xxx"):
            if verbose>0: print("ensemle.tools::fits_support - Assuming ", args[0], " is a fits-filename.")
            if "mapper" in kwargs:
                e = from_fits(x, mapper=kwargs["mapper"])
                del kwargs["mapper"]
            else:
                e = from_fits(x)
        else:
            e = x
        r = func(e, **kwargs)
        if type(x) == type("xxx"):
            to_fits(r, ofn=x, overwrite=True)
        return r    
    return wrapper


def from_fits2(fn, mapper={}, verbose=10, extension=1, maxN=None):
    """
    OBSOLETE -- Generate Astro_Ensemble from fits-file
    
    Coordinates: It is important that columns srcID, RA, and Dec exist in input file. Otherwise, populate mapper accordingly.
    
    Parameters
    ----------
    fn : str
      Filename
    mapper : dictionary
      Mapping between column-names in fits-file and property-name in Ensemlbe
    extension : int
      Fits-extension to use
    maxN : int
      Maximum number of rows to read from fits-file
      
    Returns
    -------
    ensemble of objects : Ensemble
    """
    col_mapper = lambda x:mapper[x] if x in mapper else x

    if verbose>0: print("ensemble.tools::from_fits - Reading ",fn, " (mapper: ",mapper,")")
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
        #print("ensemble.tools::from_fits - working on ",col.name)
        if col.name in mapper.values(): continue
        if col_mapper(col.name) == "srcID":
            col_data.append(ff[extension].data[col.name][0:maxN].astype(str))
            
            if col.name in mapper.keys():
                col_data.append(ff[extension].data[col.name][0:maxN])
                names.append(col.name)
        else:
            if col.name in mapper.keys():
                col_data.append(ff[extension].data[col.name][0:maxN])
                names.append(col.name)
                
            col_data.append(ff[extension].data[col.name][0:maxN])

        if verbose>6: print("ensemble.tools::from_fits - Using col: ",col.name, " (-> ", col_mapper(col.name),") with format: ",col.format)
        names.append(col_mapper(col.name))
            
    names =     ",".join(names)    
    a = np.core.records.fromarrays(col_data,names=names)
    a = np.array(a)
    
    if verbose>1: 
        print("ensemble.tools::from_fits - Read ",np.shape(a), " entries with ",len(names.split(","))," properties from ",fn)
    
    e = Ensemble().from_array(a)
    if verbose>0:
        print("ensemble.tools::from_fits - Generated Ensemble with",np.shape(a), " entries and ",len(names.split(","))," properties")
    ff.close()    
    return e
    







def from_fits(fn, mapper={}, verbose=0, extension=1, maxN=None):
    """
    Generate Astro_Ensemble from fits-file
    
    Coordinates: It is important that columns srcID, RA, and Dec exist in input file. Otherwise, populate mapper accordingly.
    
    Parameters
    ----------
    fn : str
      Filename
    mapper : dictionary
      Mapping between column-names in fits-file and property-name in Ensemlbe
    extension : int
      Fits-extension to use
    maxN : int
      Maximum number of rows to read from fits-file
      
    Returns
    -------
    ensemble of objects : Ensemble
    """
    verbose=1
    
    col_mapper = lambda x:mapper[x] if x in mapper else x

    def fits_fmt_mapper(fmt):
        if "D" in fmt:
            return "f"
        elif "I" in fmt:
            return "i"
        elif "J" in fmt:
            return "i"
        elif "A" in fmt:
            tmp = '<U'+fmt[0:-1]
            return tmp
 
    if verbose>0: print("ensemble.tools::from_fits - Reading ",fn, " (mapper: ",mapper,")")
    ff = pyfits.open(fn)    
    cols = ff[extension].columns
    if verbose>5:
        print("ensemble.tools::from_fits - ", fn, " contains ",cols)
    
    if maxN is None:
        maxN = len(ff[extension].data[cols[0].name])
        
    # Generate np.recarray for array    
    #col_data = []
    #names=[]
    dtype=[]
    # asseble dtype
    for col in cols:
        #print("ensemble.tools::from_fits - working on ",col.name)
        if col.name in mapper.values(): continue
        if col_mapper(col.name) == "srcID":
            dt = ("srcID", '<U64')
            dtype.append(dt)
            if col.name in mapper.keys():
                dt = (col.name, fits_fmt_mapper(col.format))
                dtype.append(dt)    
        else:
            if col.name in mapper.keys():
                dtype.append( (col.name, fits_fmt_mapper(col.format)) )
            dtype.append( (col_mapper(col.name),    fits_fmt_mapper(col.format)) )
            
        if verbose>6: print("ensemble.tools::from_fits - Using col: ",col.name, " (-> ", col_mapper(col.name),") with format: ",col.format, " -> dtype=",dtype[-1])
        
  
    
    xxx = np.zeros(maxN, dtype=dtype) 
    
    for col in cols:
        #print("ensemble.tools::from_fits - populating ",col.name)
        if col.name in mapper.values(): continue
        if col_mapper(col.name) == "srcID":
            xxx["srcID"] = ff[extension].data[col.name][0:maxN].astype(str)
            
            if col.name in mapper.keys(): # Keep original column if mapped
                xxx[col.name] = ff[extension].data[col.name][0:maxN]
        else:
            if col.name in mapper.keys():
                xxx[col_mapper(col.name)] = ff[extension].data[col.name][0:maxN]                
            xxx[col.name] = ff[extension].data[col.name][0:maxN]

        if verbose>6: print("ensemble.tools::from_fits - Read data from col: ",col.name, " (-> ", col_mapper(col.name),") with format: ",col.format)
    
    
    names = [dt[0] for dt in list(dtype)]
    
    if verbose>1: 
        print("ensemble.tools::from_fits - Read ",np.shape(xxx), " entries with ",len(names)," properties from ",fn)
    
    ff.close()        
    e = Ensemble().from_array(xxx)
    if verbose>0:
        print("ensemble.tools::from_fits - Generated Ensemble with",np.shape(xxx), " entries and ",len(names)," properties")
    return e
    


    
def to_fits(ensemble, ofn=None, overwrite=False, verbose=0, mapper={}, maxN=None):
    """
    Write Ensemble-data to fits-file
    
    Parameters
    ----------
    ensemble : Ensemble
      The data
    ofn : str
      Filename
    overwrite : boolean
      Overwrite file if exists
    verbose : int
      Verbosity
    maxN : int
      Maximum number of rows to write
    mapper : dictionary
      Mapping between property-name in Ensemble and in fits-file; `x` if nothing in `mapper` for property `x` and `mapper` [`x`] otherwise.
    """
    if verbose>3: print("ensemble.tools::to_fits - len, shape, len(unique): ",len(ensemble), np.shape(ensemble.array), len(np.unique(ensemble.srcIDs())))
    
    fmt_mapper = {"i":"I","u":"I","U":"64A","f":"D"}
    col_mapper = lambda x:mapper[x] if x in mapper else x

    outcols = ensemble.known_cols
    try:
        outcols.remove("coord")
    except:
        pass
    #print(*outcols, *ensemble.known_cols)

    if verbose>1:
        print("ensemble.tools::to_fits - Using cols: ",outcols)
    array = ensemble.to_array(colnames=outcols, verbose=verbose)
    #print("ARRAY")
    if maxN is None:
        maxN = len(array)
    else:
        maxN = maxN if maxN<len(array) else len(array)
        
    cols = []
    for c in array.dtype.names:
        if c=="coord": continue
        #print(c)    
        fmt = fmt_mapper[array[c].dtype.kind]
        if verbose>5: print("ensemble.tools::to_fits - Using column ",c, "with dtype=",array[c].dtype, " and fits-format=",fmt)
        col = pyfits.Column(name=col_mapper(c), array=array[c][0:maxN], format=fmt)
        cols.append(col)
        
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    if ofn is not None:
        hdul.writeto(ofn, overwrite=overwrite)        
        if verbose>0:
            print("ensemble.tools::to_fits - Written ",maxN," objects with ",len(cols)," properties to ",ofn)
    return hdul



#def merge_catalogs(ero_fn, ero_index, gaia_fn, gaia_index, ofn, extra_columns={}, overwrite=True, maxN=None):
    #"""
    
      #Parameters
      #----------
      #ero_fn : str, filename for Gaia catalog
      #ero_index : indices in the Gaia catalog
    #"""
    #if len(ero_index) != len(gaia_index):
        #print("Indices must have same lengths (", len(ero_index), len(gaia_index),")")
        #return False
    
    #gaia_ff = pyfits.open(gaia_fn)
    #ero_ff = pyfits.open(ero_fn)
    
    #gaia_columns = {"source_id":"id_gaia", "ra":"ra_gaia", "dec":"dec_gaia",\
        #"parallax":"parallax", "phot_g_mean_mag":"Gmag", "bp_rp":"bp_rp",\
        #"phot_variable_flag": "phot_variable", "teff_val": "teff_gaia", "a_g_val":"a_g"}
    #ero_columns  = {"detUID": "id_ero", "RA":"ra_ero", "DEC":"dec_ero",\
        #"RADEC_ERR": "pos_err", "EXT":"ext",\
        #"ML_FLUX_0": "ML_FLUX_0", "ML_RATE_0":"ML_RATE_0", "ML_CTS_0":"ML_CTS_0",\
        #"ML_FLUX_1": "ML_FLUX_1", "ML_RATE_1":"ML_RATE_1", "ML_CTS_1":"ML_CTS_1",\
        #"ML_FLUX_2": "ML_FLUX_2", "ML_RATE_2":"ML_RATE_2", "ML_CTS_2":"ML_CTS_2",\
        #"ML_FLUX_3": "ML_FLUX_3", "ML_RATE_3":"ML_RATE_3", "ML_CTS_3":"ML_CTS_3"}
    
    #hdu = pyfits.PrimaryHDU()    
    #cols = []
    
    #if maxN:
        #gaia_index = gaia_index[0:maxN]
        #ero_index= ero_index[0:maxN]
    #else:
        #maxN = len(gaia_index)  
    #print("gaia_index length",len(gaia_index))
    
    #for c in gaia_columns:
        #fmt = "E"
        #if c == "source_id": fmt = "K"
        #elif c == "ra": fmt="D"
        #elif c == "dec": fmt="D"
        #elif c == "phot_variable_flag": fmt="16A"
        #print(c, fmt)
        #col = pyfits.Column(name=gaia_columns[c], array=gaia_ff[1].data[c][gaia_index], format=fmt)
        #cols.append(col)
    
    #for c in ero_columns:
        #fmt = "E"
        #if c == "detUID": fmt = "32A"
        #elif c == "RA": fmt="D"
        #elif c == "DEC": fmt="D"    
        #elif c == "source_type": fmt="2A" 
        #col = pyfits.Column(name=ero_columns[c], array=ero_ff[1].data[c][ero_index], format=fmt)
        #cols.append(col)
    
    #for ec in extra_columns:
        #fmt="E"
        #if ec == "source_type": fmt="2A" 
        #col = pyfits.Column(name=ec, array=extra_columns[ec][0:maxN], format=fmt)
        #cols.append(col)
        
    #cc = pyfits.ColDefs(cols)
    #xx = pyfits.BinTableHDU.from_columns(cc)
    ##xx.header["min_offset"] = min_offset
    #hdul = pyfits.HDUList([hdu, xx])
    #hdul.writeto(ofn, overwrite=overwrite)
