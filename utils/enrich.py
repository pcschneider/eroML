from eroML.ensemble import Ensemble, from_fits,to_fits, fits_support
import numpy as np
import astropy.units as u
from .iso_tools import add_iso_column
from .gaia_tools import add_quality_column
from .estimators import sky_dens4coordinates
import logging

logger = logging.getLogger('eroML')

activity_poly = [-3.22, 3.6/5.5] # in log scale

def activity_filter(color, FxFg, log_margin=0):
    """
      Parameters
      ----------
      color : array (Gaia Bp-Rp color)
      FxFg : array (ratio, linear)
      
      Returns 
      -------
        array : True for sources above threshold
    """
    ys = activity_poly[1]
    ye = color*ys + activity_poly[0]
    dy = np.log10(FxFg)-ye
    return (dy>log_margin).astype(int)

@fits_support
def eligible_Gaia(e, out_col="eligible_Gaia", verbose=5):
    """
    Add property reflecting if the source could be an eligible stellar counter part.
    """
    iso = e.to_array(colnames="iso_compatible", array_type="array")
    gaia = e.to_array(colnames="Gaia_quality", array_type="array")
    #iso = np.array([True if ii=="True" else False for ii in iso])
    #gaia = np.array([True if ii=="True" else False for ii in gaia])
    gi = np.where( (iso==1) & (gaia==1) )[0]
    logger.debug("Number of eligible sources: %i (%5.3f%%)" % (len(gi),len(gi)/len(e) * 100))
    el = np.zeros(len(e))
    el[gi] = 1

    if out_col not in e.known_cols: e.add_col(out_col, el.astype(int))
    else: e.set_col(out_col, el.astype(int))

    return e


@fits_support
def eligible_eROSITA(e, out_col="eligible_X"):
    det_likeli = e.to_array(colnames="DET_LIKE_0", array_type="array")
    
    ext_likeli = e.to_array(colnames="EXT_LIKE", array_type="array")
    gi = np.where((det_likeli > 6) & (ext_likeli < 6))[0]
    
    
    el = np.zeros(len(e))
    el[gi] = 1
    logger.debug("Number of eligible sources: %i (%5.3f%%)" % (len(gi),len(gi)/len(e) * 100))

    if out_col not in e.known_cols: e.add_col(out_col, el.astype(int))
    else: e.set_col(out_col, el.astype(int))

    return e
    


@fits_support
def eligible_ROSAT(e, out_col="eligible_X"):
    ext_ml = e.to_array(colnames="EXT_ML", array_type="array")
    det_likeli = e.to_array(colnames="EXI_ML", array_type="array")
    gi = np.where((det_likeli > 6) & (ext_ml == 0))[0]
    el = np.zeros(len(e))
    el[gi] = 1
    logger.debug("Number of eligible sources: %i (%5.3f%%)" % (len(gi),len(gi)/len(e) * 100))

    if out_col not in e.known_cols: e.add_col(out_col, el.astype(int))
    else: e.set_col(out_col, el.astype(int))

    return e

@fits_support
def sky_density(e, around=5, filter_prop="eligible_Gaia", filter_value=1, out_col="eligible_sky_density", verbose=1):
    """
    Calculate the local sky density for each source
    
    The sky density is in #stars/arcmin^-2
    
    Parameters
    ----------
    e : Ensemble
    around : float
        The on-sky radius in arcmin
    filter_prop : str
        The column to filter the input sample, i.e., use only entries which have the given `filter_value`.
    filter_value : int
    outcol : str
        Name for the new column containing the sky density
    """
    
    dens = np.zeros(len(e))
    dens[:] = np.nan
    
    if filter_prop is not None:
        gi = np.where(e.to_array(colnames=filter_prop, array_type="array") == filter_value)[0]
    else:
        gi = np.arange(len(e))
    
    
    coord = e.skyCoords()[gi]

    dens[gi] = sky_dens4coordinates(coord, around=around)
    
    if out_col not in e.known_cols: e.add_col(out_col, dens)
    else: e.set_col(out_col, dens)
    print("sky_dens::  outcol: ",out_col," nanmean: ",np.nanmean(dens))
    return e
    


@fits_support
def NN_Max(e):
    """
    Count number of matches for each source
    """
    sids = e.srcIDs()
    oids = e.to_array(colnames="original_srcID", array_type="array")
    #print(sids, oids)
    cnt = np.zeros(len(sids))
    unique_elements, counts_elements = np.unique(oids, return_counts=True)
    for a, b in zip(unique_elements, counts_elements):
        gi = np.where(oids == a)[0]
        #print(a, gi, b)
        cnt[gi] = b
    if "NN_max" in e.known_cols:
        e.set_col("NN_max", cnt.astype(int))
    else:
        e.add_col("NN_max", cnt.astype(int))
    return e

 
@fits_support
def calc_gaia_quality(e, colname="Gaia_quality", overwrite=False, verbose=10, filter_Nr=2):
    """
    """
     
    q = quality_filter(ff, filter_Nr=filter_Nr)
    c = pyfits.Column(name=colname, array=q, format="L")
    if verbose>1: print("gaia_tools::add_quality_column - #rows: %i, good quality: %i, fraction: %f" % (len(ff[1].data[cols[0].name]), np.sum(q),np.sum(q)/len(ff[1].data[cols[0].name]) ))
    cols.add_col(c)
    hdu = copy.copy(ff[0]) 
    #hdu = ff[0]
    #hdu = pyfits.PrimaryHDU()
    hdx = pyfits.BinTableHDU.from_columns(cols)
    hdul = pyfits.HDUList([hdu, hdx])
    hdul.writeto(ofile, overwrite=overwrite)
    ff.close()
 
 
@fits_support
def enrich_Gaia(e, filterNr=5, calc_sky_density=True):
    """
    Add G-band flux, compatibility with an isochrone, the `eligible` column, and the sky density of eligible sources
    """
    arr = e.to_array(colnames=["phot_g_mean_mag"])
    FG = 10**(-0.4* arr["phot_g_mean_mag"])*1.01324e-5
    e.add_col("Fg", FG, force=True)
    
    add_quality_column(e, filter_Nr=filterNr)    
    add_iso_column(e)        
    eligible_Gaia(e)
    
    if calc_sky_density:
        sky_density(e, around=3, filter_prop="eligible_Gaia", filter_value=1, out_col="eligible_sky_density")
    #sky_density(e, around=3, filter_prop=None, out_col="sky_density")
    return e

    
@fits_support
def enrich_eROSITA(e):
    """
    Enrich the eROSITA data.
    
    Add:
      - Fx
      - eligible_X
    as well as the dummy columns:
      - pm_RA, pm_Dec
      - ref_epoch
    """
    err = e.to_array(colnames="RADEC_ERR", array_type="array")
    err[err<1.] = 1.
    e.set_col("RADEC_ERR", err)

    Fx = e.to_array(colnames="ML_FLUX_0", array_type="array")
    if "Fx" in e.known_cols:
        e.set_col("Fx", Fx)
    else:
        e.add_col("Fx", Fx)
     
    t0 = e.to_array(colnames="TSTART", array_type="array")
    t1 = e.to_array(colnames="TSTOP", array_type="array") 
    
    e.add_col("pm_RA", np.zeros(len(err)))
    e.add_col("pm_Dec",np.zeros(len(err)))
    e.add_col("ref_epoch", np.ones(len(err))*2020.25)
    
    eligible_eROSITA(e)    
    return e



@fits_support
def enrich_ROSAT(e):
    """
    Enrich the ROSAT data.
    
    Add:
      - Fx
      - eligible_X
    as well as the dummy columns:
      - pm_RA, pm_Dec
      - ref_epoch
    """
    Xerr = e.to_array(colnames="XERR", array_type="array")
    Yerr = e.to_array(colnames="YERR", array_type="array")
    err = np.sqrt(Xerr**2+Yerr**2) * 45 / 2.
    err[err<3.] = 3.
    e.add_col("RADEC_ERR", err)

    HR = e.to_array(colnames="HR_1", array_type="array")
    CR = e.to_array(colnames="RATE", array_type="array")
    Fx = (5.3*HR + 8.31) * 1e-12
    ni = np.where(np.isfinite(Fx) == False)[0]
    Fx[ni] = 6e-12 * CR[ni]
    
    if "Fx" in e.known_cols:
        e.set_col("Fx", Fx)
    else:
        e.add_col("Fx", Fx)
     
    #t0 = e.to_array(colnames="TSTART", array_type="array")
    #t1 = e.to_array(colnames="TSTOP", array_type="array") 
    
    e.add_col("pm_RA", np.zeros(len(err)))
    e.add_col("pm_Dec",np.zeros(len(err)))
    e.add_col("ref_epoch", np.ones(len(err))*1991.5)
    
    eligible_ROSAT(e)    
    return e



#@fits_support
#def enrich_ROSAT(e):
    #"""
    #Enrich the ROSAT sources.
    
    #Add:
      #- Fx
      #- eligible_X
    #as well as the dummy columns:
      #- pm_RA, pm_Dec
      #- ref_epoch
    #"""
    #err = e.to_array(colnames="RADEC_ERR", array_type="array")
    #err[err<1.] = 1.
    #e.set_col("RADEC_ERR", err)

    #Fx = e.to_array(colnames="ML_FLUX_0", array_type="array")
    #if "Fx" in e.known_cols:
        #e.set_col("Fx", Fx)
    #else:
        #e.add_col("Fx", Fx)
     
    #t0 = e.to_array(colnames="TSTART", array_type="array")
    #t1 = e.to_array(colnames="TSTOP", array_type="array") 
    
    #e.add_col("pm_RA", np.zeros(len(err)))
    #e.add_col("pm_Dec",np.zeros(len(err)))
    #e.add_col("ref_epoch", np.ones(len(err))*2020.25)
    
    #eligible_eROSITA(e)    
    #return e


@fits_support
def enrich_merged(e):
    offset_sig = e.to_array(colnames="RADEC_ERR", array_type="array")
    d2d = e.to_array(colnames="match_dist", array_type="array")
    e.add_col("offset_sig", d2d/offset_sig)
    
    Fx = e.to_array(colnames="Fx", array_type="array")
    Fg = e.to_array(colnames="Fg", array_type="array")
    FxFg = Fx/Fg
    e.add_col("FxFg", FxFg)
    
    color = e.to_array(colnames="bp_rp", array_type="array")
    activity = activity_filter(color, FxFg)
    e.add_col("too_active", activity)
    return e
    
    
    
    #ec={"offset_sig": np.concatenate([d2d.arcsec/err, d2d2[gi2].arcsec/err[gi2], d2d3[gi3].arcsec/err[gi3]])}

