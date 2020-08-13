from eroML.ensemble import Ensemble, from_fits,to_fits, fits_support
import numpy as np
import astropy.units as u
from .iso_tools import add_iso_column
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
def eligible(e, out_col="eligible", verbose=5):
    """
    Add property reflecting if the source could be an eligible stellar counter part.
    """
    iso = e.to_array(colnames="iso_compatible", array_type="array")
    gaia = e.to_array(colnames="Gaia_quality", array_type="array")
    iso = np.array([True if ii=="True" else False for ii in iso])
    gaia = np.array([True if ii=="True" else False for ii in gaia])

    gi = np.where( (iso) & (gaia) )[0]
    print("enrich::eligible - Number of eligible sources: ",len(gi)," (",len(gi)/len(e) * 100,"%)")
    el = np.zeros(len(e))
    el[gi] = 1

    if out_col not in e.known_cols: e.add_col(out_col, el.astype(int))
    else: e.set_col(out_col, el.astype(int))

    return e


@fits_support
def sky_density(e, around=3, filter_prop="eligible", filter_value=1, out_col="eligible_sky_density", verbose=1):
    """
    
    Parameters
    ----------
    e : Ensemble
    around : float
        The on-sky radius in arcmin
    """
    #return None
    dens = np.zeros(len(e))
    dens[:] = np.nan
    
    if filter_prop is not None:
        gi = np.where(e.to_array(colnames=filter_prop, array_type="array") == filter_value)[0]
    else:
        gi = np.arange(len(e))
    
    
    coord = e.skyCoords()[gi]
    if verbose>0:
        print("sky_dens:: filter_prop: ",filter_prop, " out_col:",out_col)    
        print("sky_dens:: Searching around ",around, "arcmin.")
        print("sky_dens:: Coord range: (",max(coord.ra.degree), min(coord.ra.degree), " ; ", max(coord.dec.degree), min(coord.dec.degree),")")
    
    
    ra_range0 = abs(max(coord.ra.degree) - min(coord.ra.degree))
    ra_range1 = abs(max((coord.ra.degree  + 180) % 360)- min((coord.ra.degree + 180) % 360))
    #print("ra_range0, ra_range1",ra_range0, ra_range1)
    ra_range = ra_range0 if ra_range0<ra_range1 else ra_range1
    dec_range = abs(max(coord.dec.degree) - min(coord.dec.degree))
    sky_area = ra_range *dec_range * np.cos(np.nanmean(coord.dec.degree)/180*np.pi)
    #print("cos", np.cos(np.nanmedian(coord.dec.degree/180*np.pi)))
    sky_dens = len(coord)/sky_area/3600 # per armin^2
    if verbose>0: 
        print("sky_dens:: Sky area of Ensemble: ",sky_area, " (center: ",np.nanmedian(coord.ra.degree), np.nanmedian(coord.dec.degree),")")         
        print("sky_dens::    Mean sky density: ",sky_dens," #stars/arcmin^2")
    
    if sky_dens > 30:
        #return None
        print("sky_dens:: splitting Ensemble...")
        s = e.split(3)
        f = Ensemble()
        for ss in s:
            x = sky_density(ss, around=around, filter_prop=filter_prop, filter_value=filter_value, out_col=out_col)
            skd = x.to_array(out_col, array_type="array")
            x.set_col(out_col, skd/3)
            f.append(x)
        return f
    
    elif sky_dens > 10:
        ta = 0.5*u.arcmin
        idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, ta)
        uni, cnt = np.unique(idxc, return_counts=True)
        cnt=np.array(cnt)*(around/0.5)**2
        print("sky_dens:: Using 0.5 arcmin search radius and extrapolating.")        
    #elif sky_dens > 3:
        #ta = 0.7*u.arcmin
        #idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, ta)
        #uni, cnt = np.unique(idxc, return_counts=True)
        #cnt=np.array(cnt)*(around/0.7)**2
        #print("sky_dens:: Using 0.7 arcmin search radius and extrapolating.")
    else:
        idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, around*u.arcmin)
        uni, cnt = np.unique(idxc, return_counts=True)
        print("sky_dens:: Using ",around," arcmin search radius.")

    dens[gi] = 1/cnt
    if out_col not in e.known_cols: e.add_col(out_col, dens)
    else: e.set_col(out_col, dens)
    print("sky_dens::  outcol: ",out_col," nanmean: ",np.nanmean(dens), "returning ",e)
    return e
        
@fits_support
def enrich_Gaia(e):
    """
    Add G-band flux
    """
    arr = e.to_array(colnames=["phot_g_mean_mag"])
    FG = 10**(-0.4* arr["phot_g_mean_mag"])*3.660e-08*720
    e.add_col("Fg", FG, force=True)
    
    
    add_iso_column(e)        
    eligible(e)
    sky_density(e, around=3, filter_prop="eligible", filter_value=1, out_col="eligible_sky_density")
    sky_density(e, around=3, filter_prop=None, out_col="sky_density")
    return e
    
    #ra = e.to_array(colnames=["ra"], array_type='array')
    #dec = e.to_array(colnames=["dec"], array_type='array')
    ##print(ra, type(ra))
    #e.add_col("RA", ra)
    #e.add_col("Dec", dec)
    
@fits_support
def enrich_eROSITA(e):
    err = e.to_array(colnames="RADEC_ERR", array_type="array")
    err[err<1.] = 1.
    e.set_col("RADEC_ERR", err)

    Fx = e.to_array(colnames="ML_FLUX_0", array_type="array")
    if "Fx" in e.known_cols:
        e.set_col("Fx", Fx)
    else:
        e.add_col("Fx", Fx)
    return e

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

