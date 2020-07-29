import functools
from eroML.ensemble import from_fits,to_fits
import numpy as np

def fits_support(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if len(args)>0:
            x = args[0]
        else:
            raise TypeError("fits_support: Expecting at least one argument, but none provided.")
        if type(x) == type("xxx"):
            print("enrich.py::fits_support - Assuming ", args[0], " is a fits-filename.")
            if "mapper" in kwargs:
                e = from_fits(x, mapper=kwargs["mapper"])
            else:
                e = from_fits(x)
        else:
            e = x
        func(e)
        if type(x) == type("xxx"):
            to_fits(e, ofn=x, overwrite=True)
    return wrapper

        
@fits_support
def enrich_Gaia(e):
    arr = e.to_array(colnames=["phot_g_mean_mag"])
    FG = 10**(-0.4* arr["phot_g_mean_mag"])*3.660e-08*720
    e.add_col("F_G", FG)
    
    ra = e.to_array(colnames=["ra"], array_type='array')
    dec = e.to_array(colnames=["dec"], array_type='array')
    #print(ra, type(ra))
    e.add_col("RA", ra)
    e.add_col("Dec", dec)
    
@fits_support
def enrich_eROSITA(e):
    err = e.to_array(colnames="RADEC_ERR", array_type="array")
    err[err<1.] = 1.
    e.set_col("RADEC_ERR", err)

@fits_support
def enrich_merged(e):
    offset_sig = e.to_array(colnames="RADEC_ERR", array_type="array")
    d2d = e.to_array(colnames="match_dist", array_type="array")
    e.add_col("offset_sig", d2d/offset_sig)
    #ec={"offset_sig": np.concatenate([d2d.arcsec/err, d2d2[gi2].arcsec/err[gi2], d2d3[gi3].arcsec/err[gi3]])}

