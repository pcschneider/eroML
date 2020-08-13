from eroML.ensemble import Ensemble
from eroML.ensemble import from_fits,to_fits,multi_fits_support, fits_support
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import copy
from .gaia_tools import gaia4ero
from .estimators import NN_distribution
from .enrich import enrich_merged, activity_filter
import glob

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
    
@multi_fits_support(3)
def major_set(ero0, gaia, keep_ero_cols=None, keep_gaia_cols=None, NN=3, verbose=10, overwrite=True):
    """
    Merge eROSITA and Gaia sources into one fits-file. 
    
    Matching is based on the sky-distance. The `NN`-argument determines the n-th nearest match that enters 
    the merged table. Since sources must be unique, the second nearest neighbor will be added as `srcID`_NN2, i.e., 
    a new source is created. The original srcID is kept in the field `original_srcID`.
    
    Parameters
    ----------
    keep_[ero/gaia] cols : refers to respective ifn
    """


    #gaia = from_fits(ifn_gaia, mapper={"source_id":"srcID", "ra":"RA", "dec":"Dec"}, maxN=20000)
    #gaia = from_fits(ifn_gaia, maxN=20000)
    #print("gaia known cols: ",gaia.known_cols)
    #print("Gaia len",len(gaia), np.shape(gaia.array))
    #enrich_Gaia(gaia)
    
    #ero0 = from_fits(ifn_ero, mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=2000)
    #ero0 = from_fits(ifn_ero, maxN=2000)
    #enrich_eROSITA(ero0)
    #print("len 0: ", len(ero0))
    #print("XXXXXXXX", len(ero0), len(gaia))
    #print(type(ero0), type(gaia))
    #return copy.deepcopy(ero0)  
    #return
          
    eros = []
    for i in range(NN):
        if verbose>1: print("datasets::major_set - Merging NN=",i+1)
        ero1 = copy.deepcopy(ero0)
        ero1.merge_add(gaia, NN=i+1)
        ero1.add_col("NN", np.array(len(ero1)*[i+1]))
        ero1.add_col("original_srcID", np.array(ero1.srcIDs()))
        eros.append(ero1)
    ero1 = eros[0]
    
    for i, e in enumerate(eros[1::]):
        ero1.append(e, postfix="_NN"+str(i+2))
    enrich_merged(ero1)    
    
    #print("len 1: ", len(eros[0]))
    #gi2 = np.where(d2d2.arcsec < 5*err)[0]
    #gi3 = np.where(d2d3.arcsec < 5*err)[0]

    offset_sig = ero1.to_array(colnames="offset_sig", array_type="array")
    d2d = ero1.to_array(colnames="match_dist", array_type="array")
    gi = np.where((offset_sig<50) & (d2d<1000))[0]
    srcIDs = np.array(ero1.srcIDs())
    #print(len(gi), gi, srcIDs )
    good_ids = srcIDs[gi]
    #print("Keeping: ",good_ids)
    ero1.keep(good_ids)
    
    NN_Max(ero1)
    
    #sids = ero1.srcIDs()
    #oids = ero1.to_array(colnames="original_srcID", array_type="array")
    ##print(sids, oids)
    #cnt = np.zeros(len(sids))
    #unique_elements, counts_elements = np.unique(oids, return_counts=True)
    #for a, b in zip(unique_elements, counts_elements):
        #gi = np.where(oids == a)[0]
        ##print(a, gi, b)
        #cnt[gi] = b
    #ero1.add_col("NN_max", cnt.astype(int))
    
        
    #print(ero1.to_array(colnames="srcID_NN"))
    
    #to_fits(ero1, ofn, overwrite=True)
    return ero1


@multi_fits_support(2)
def training_set(e, abs_dist_cutoff=3, rel_dist_cutoff=2):
    """
    Generate training set based on good positional matches and a criterium on Fx/Fg
    
    Parameters
    ----------
    e : Ensemble
    abs/rel_dist_cutoff : float
        Use only sources, which are closer on sky than these cutoff values.
        
        Units: 
          - abs_dist_cutoff : arcsec
          - rel_dist_cutoff : sigma
    """
    
    #print("XXXXXXXXX")
    abs_dist = e.to_array("match_dist", array_type="array")
    rel_dist = e.to_array("offset_sig", array_type="array")
    
    gi = np.where( (abs_dist < abs_dist_cutoff) & (rel_dist < rel_dist_cutoff) )[0]
    
    
    color = e.to_array("bp_rp", array_type="array")
    FxFg = e.to_array("FxFg", array_type="array")
    below = (1-activity_filter(color, FxFg)).astype(bool)
    well_above = activity_filter(color, FxFg, log_margin=0.5).astype(bool)
    #print(len(below), len(well_above), len(e))
    
    cl = np.ones(len(e)).astype(int) * 2
    cl[below] = 1
    cl[well_above] = 0
    #print("unique: ",np.unique(cl), below)
    gi2 = np.where(cl < 2)[0]
    
    
    #print("#below acitvity threshold: ",np.sum(below), ", well above: ",np.sum(well_above), " all: ",len(e))
    #print(gi2, len(e))
    #print(type(gi), type(gi2), len(gi), len(gi2))
    keep = np.intersect1d(gi, gi2).astype(int).tolist()
        
    #print("keep",len(keep))
    
    e.add_col("category", cl)
    
    sids = np.array(e.srcIDs())
    #print("keep",keep, type(keep), " sids",sids, type(sids))
    #print(sids[keep])
    e.keep(sids[keep])
    
    
    #print("YYYYY")
    NN_Max(e)
    return e
    
    
@multi_fits_support(3)    
def training_random_set(training, random, abs_dist_cutoff=3, rel_dist_cutoff=2):
    """
    Merge the training set with the random set applying the same selection criteria.
    
    Random sources get `category` 2 assigned.
    
    Parameters
    ----------
    training, random : Ensemble
    abs/rel_dist_cutoff : float
        In arcsec
        
    Returns
    -------
    Merged : Ensemble
    
    """
    random.add_col("category", np.array(len(random)*[2]))
    rand_srcIDs = np.array(random.srcIDs())
    abs_dist = random.to_array("match_dist", array_type="array")
    rel_dist = random.to_array("offset_sig", array_type="array")
    gi = np.where( (abs_dist <=abs_dist_cutoff) & (rel_dist <= rel_dist_cutoff))[0]
    random.keep(rand_srcIDs[gi])
    training.append(random, postfix="rand")
    return training
    
    
@multi_fits_support(3)
def random_set(ero, gaia):
    """
    Generate a random set, i.e., shuffle the X-ray positions
   
    
      Parameters
      ----------
      ifn, ofn : string, filenames for in- and output files
      offsets : float, in arcsec
      overwrite : boolean
    """
    min_offset = 60 #arcsec
    max_offset = 180
    
    #gaia = from_fits(ifn_gaia, maxN=20000)
    #print("gaia known cols: ",gaia.known_cols)
    #print("Gaia len",len(gaia), np.shape(gaia.array))
    #enrich_Gaia(gaia)
    
    #ero0 = from_fits(ifn_ero, mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=2000)
    #ero = from_fits(ifn_ero, maxN=2000)
    
    #print(len(ra))
    coords = ero.skyCoords()
    #print(len(coords))
    N = len(coords)
    angle = np.random.rand(N)*360
    offset = np.random.rand(N)*(max_offset- min_offset)+min_offset
    #print(offset)

    doords = coords.directional_offset_by(angle*u.degree, offset*u.arcsec)

    nra, ndec = doords.ra.degree, doords.dec.degree
    ero.set_col("RA", nra)
    ero.set_col("Dec", ndec)
    
    return major_set(ero, gaia)
    
def prep_classify(ifn, extension=1, ofn=None, overwrite=False, verbose=1):
    """
    Keep only relevant columns
    """
    relevant_cols = ["srcID", "Fx","Fg","bp_rp","offset_sig","parallax","sky_density"]
    
    ff = pyfits.open(ifn)
    
    cols = []
    columns = {col.name:col.format for col in ff[extension].columns}
    
    col = pyfits.Column(name="srcID", array=ff[extension].data["srcID"], format=columns["srcID"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fx"])+13
    print("log Fx",np.mean(arr), np.median(arr), np.std(arr))
    col = pyfits.Column(name="logFx", array=arr, format=columns["Fx"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["Fg"])+12
    print("log Fg",np.mean(arr), np.median(arr), np.std(arr))
    col = pyfits.Column(name="logFg", array=arr, format=columns["Fg"])    
    cols.append(col)
    
    arr = ff[extension].data["bp_rp"]
    print("bp_rp",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="bp_rp", array=arr, format=columns["bp_rp"])    
    cols.append(col)
    
    dst = ff[extension].data["match_dist"]
    err = ff[extension].data["RADEC_ERR"]
    val = np.zeros(len(dst))
    for i, (x, s) in enumerate(zip(dst, err)):
        xxx = np.arange(0,x,0.001)
        y = len(xxx)*[0.1]
        y = xxx/s**2 * np.exp(-xxx**2/(2*s**2))
        val[i] = 1- np.trapz(y, xxx)
    #print(val)
    #import matplotlib.pyplot as plt
    #plt.scatter(ff[extension].data["offset_sig"], val)
    #plt.show()
    print("pos",np.nanmean(val), np.nanmedian(val), np.nanstd(val))
    col = pyfits.Column(name="pos", array=val, format=columns["match_dist"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["parallax"])
    gi = np.where(np.isnan(arr))[0]
    arr[gi] = -1
    print("plx",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_plx", array=arr, format=columns["parallax"])    
    cols.append(col)
    
    arr = np.log10(ff[extension].data["sky_density"])
    print("sky_density",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="log_sk", array=arr, format=columns["sky_density"])    
    cols.append(col)
        
    arr = ff[extension].data["category"]
    print("category",np.nanmean(arr), np.nanmedian(arr), np.nanstd(arr))
    col = pyfits.Column(name="category", array=arr, format=columns["category"])    
    cols.append(col)
         
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    hdul = pyfits.HDUList([hdu, xx])
    print("XXXXXX")
    if ofn is not None:
        hdul.writeto(ofn, overwrite=overwrite)        
        if verbose>0:
            print("datasets::prep_classify - Written ",len(dst)," objects with ",len(cols)," properties to ",ofn)
    ff.close()            
    return hdul
