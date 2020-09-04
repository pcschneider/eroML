from eroML.ensemble import Ensemble
from eroML.ensemble import from_fits,to_fits,multi_fits_support, fits_support
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import copy
from .gaia_tools import gaia4ero
from .estimators import NN_distribution
from .enrich import enrich_merged, activity_filter, NN_Max
import glob
import logging

logger = logging.getLogger('eroML')
    
@multi_fits_support(2)    
def shrink(e, cols=[]):
    """
    """
    standard_cols = ["srcID", "RA", "Dec"]
    merged_cols = np.unique(standard_cols + cols)
    #print(merged_cols)
    tmp = e.to_array(colnames=standard_cols)
    f = Ensemble()
    f.from_array(tmp)
    for c in e.known_cols:
        if c in merged_cols:
          f.add_col(c, e.array[c])
    return f      



def file_loop_1to1(idx, prefix="", postfix="", ofn_prefix="", ofn_postfix="", method=None, **kwargs):
    for j,i in enumerate(idx):
        fn = prefix+str(i)+postfix+".fits"
        ofn = ofn_prefix+str(i)+ofn_postfix+".fits"
        logger.debug("Creating data set for fn=%s (ofn=%s; file# %i/%i)." % (fn, ofn, j+1, len(idx)))
        method(fn, ofn, **kwargs)


def file_loop_2to1(idx, prefix1="", postfix1="", prefix2="", postfix2="", ofn_prefix="", ofn_postfix="", method=None, **kwargs):
    for j, i in enumerate(idx):
        fn1 = prefix1+str(i)+postfix1+".fits"
        fn2 = prefix2+str(i)+postfix2+".fits"
        ofn = ofn_prefix+str(i)+ofn_postfix+".fits"
        logger.debug("Creating data set for fn1=%s and fn2=%s (ofn=%s; file# %i/%i)." % (fn1, fn2, ofn, j+1, len(idx)))
        method(fn1, fn2, ofn, **kwargs)

def major_loop(idx, ero_prefix=None, ero_postfix=None, gaia_prefix=None, gaia_postfix=None, major_prefix=None, major_postfix=None):
    """
    """
    for i in idx:
        ero_fn = ero_prefix+str(i)+ero_postfix+".fits"
        gaia_fn = gaia_prefix+str(i)+gaia_postfix+".fits"
        ofn = major_prefix+str(i)+major_postfix+".fits"
        logger.debug("Creating major set for ero=%s and Gaia=%s (ofn=%s)." % (ero_fn, gaia_fn, ofn))
        major_set(ero_fn, gaia_fn, ofn)



def random_loop(idx, ero_prefix=None, ero_postfix=None, gaia_prefix=None, gaia_postfix=None, random_prefix=None, random_postfix=None, multi=1, min_offset=60, max_offset=180):
    """
    """
    for i in idx:
        ero_fn = ero_prefix+str(i)+ero_postfix+".fits"
        gaia_fn = gaia_prefix+str(i)+gaia_postfix+".fits"
        ofn = random_prefix+str(i)+random_postfix+".fits"
        logger.debug("Creating random data set for ero=%s and Gaia=%s (ofn=%s)." % (ero_fn, gaia_fn, ofn))
        random_set(ero_fn, gaia_fn, ofn, multi=multi, min_offset=min_offset, max_offset=max_offset)

def training_loop(idx, major_prefix=None, major_postfix=None, random_prefix=None, random_postfix=None, training_prefix=None, training_postfix=None, abs_dist=2, rel_dist=1):
    """
    """
    for j, i in enumerate(idx):
        major_fn = major_prefix+str(i)+major_postfix+".fits"
        random_fn = random_prefix+str(i)+random_postfix+".fits"
        ofn = training_prefix+str(i)+training_postfix+".fits"
        logger.debug("Creating training data set for major=%s and random=%s (ofn=%s; file# %i/%i)." % (major_fn, random_fn, ofn, j+1, len(idx)))
        training_set(major_fn, random_fn, ofn, abs_dist_cutoff=abs_dist, rel_dist_cutoff=rel_dist)

        
@multi_fits_support(3)
def major_set(ero, gaia, eligible_ero="eligible_eROSITA", eligible_gaia="eligible_Gaia", NN=3, verbose=10, overwrite=True):
    """
    Merge eROSITA and Gaia sources into one fits-file. 
    
    Matching is based on the sky-distance. The `NN`-argument determines the n-th nearest match that enters 
    the merged table. Since sources must be unique, the second nearest neighbor will be added as `srcID`_NN2, i.e., 
    a new source is created. The original srcID is kept in the field `original_srcID`.
    
    Parameters
    ----------
    keep_[ero/gaia] cols : refers to respective ifn
    """         
    
    ero0 = copy.deepcopy(ero)
    ero_e = ero0.to_array(colnames=eligible_ero, array_type="array")
    ero_ids = np.array(ero0.srcIDs())
    gi = np.where(ero_e == 1)[0]
    ero0.keep(ero_ids[gi])
    logger.debug("Keeping %i of % i as eligble eROSITA sources." % (len(gi), len(ero)))
    
    gaia0 = copy.deepcopy(gaia)
    gaia_e = gaia0.to_array(colnames=eligible_gaia, array_type="array")
    gaia_ids = np.array(gaia0.srcIDs())
    gi = np.where(gaia_e == 1)[0]
    gaia0.keep(gaia_ids[gi])
    logger.debug("Keeping %i of % i as eligble Gaia sources." % (len(gi), len(gaia)))
    
    eros = []
    for i in range(NN):
        if verbose>1: print("datasets::major_set - Merging NN=",i+1)
        ero1 = copy.deepcopy(ero0)
        
        ero1.merge_add(gaia0, NN=i+1, epoch=2019)
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
    return ero1


@multi_fits_support(3)
def training_set(major0, random0, abs_dist_cutoff=3, rel_dist_cutoff=2):
    major = copy.deepcopy(major0)
    random = copy.deepcopy(random0)
    
    #print("len", len(major), len(random))
    
    for e in [major, random]:
        abs_dist = e.to_array("match_dist", array_type="array")
        rel_dist = e.to_array("offset_sig", array_type="array")
        #print(abs_dist)
        gi = np.where( (abs_dist < abs_dist_cutoff) & (rel_dist < rel_dist_cutoff) )[0]
        sids = np.array(e.srcIDs())[gi]
        e.keep(sids)
    
    random.add_col("category", np.ones(len(random))*2)
    
    color = major.to_array("bp_rp", array_type="array")
    FxFg = major.to_array("FxFg", array_type="array")
    below = (1-activity_filter(color, FxFg)).astype(bool)
    well_above = activity_filter(color, FxFg, log_margin=0.5).astype(bool)
    #print(len(below), len(well_above))
    #print(below, well_above)
    
    cl = np.ones(len(major)).astype(int) *2
    cl[below] = 0
    cl[well_above] = 1
    gi2 = np.where(cl < 2)[0]
    major.add_col("category", cl)
    sids = np.array(major.srcIDs())[gi2]
    #print(type(major), len(major), len(sids))
    major.keep(sids)
    
    major.append(random, postfix="_rnd")
    #print(len(major))
    return major
    
    #print("unique: ",np.unique(cl), below)
    
    
    
    

@multi_fits_support(2)
def training_set2(e, abs_dist_cutoff=3, rel_dist_cutoff=2):
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
def random_set(ero0, gaia, multi=1, min_offset=60, max_offset=180):
    """
    Generate a random set, i.e., shuffle the X-ray positions
   
    
    Parameters
    ----------
    ifn, ofn : string, filenames for in- and output files
    offsets : float, in arcsec
    multi : int
          Increase size of sample by N
    overwrite : boolean
    """
    
    #min_offset = 60 #arcsec
    #max_offset = 180
    
    #gaia = from_fits(ifn_gaia, maxN=20000)
    #print("gaia known cols: ",gaia.known_cols)
    #print("Gaia len",len(gaia), np.shape(gaia.array))
    #enrich_Gaia(gaia)
    
    #ero0 = from_fits(ifn_ero, mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=2000)
    #ero = from_fits(ifn_ero, maxN=2000)
    logger.debug("Random data set for min_offset=%f, max_offset=%f, multi=%i" % (min_offset, max_offset, multi))
    ero = copy.deepcopy(ero0)
    
    #NN_random = np.array(len(tmp)*multi*[1])
    NN_random = np.ones(len(ero)*multi)
    N = len(ero)
    
    osids = ero.srcIDs()
    
    #print(N)
    for i in range(multi-1):
        tmp = copy.deepcopy(ero0)
        ero.append(tmp, postfix="_rnd"+str(i+2))
        i0, i1 = (i+1)*N, (i+2)*N-1
        #print(i0, i1)
        NN_random[i0:i1] = i+2

    for si in osids:
        ero.rename(si, si+"_rnd1")
        
        
    ero.add_col("N_random", NN_random)
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
    
