from eroML.ensemble import Ensemble
from eroML.ensemble import from_fits,to_fits
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import copy
from .gaia_tools import gaia4ero
from .estimators import NN_distribution
from .enrich import enrich_merged

def major_set(ifn_ero, ifn_gaia, ofn, keep_ero_cols=None, keep_gaia_cols=None, NN=3, verbose=1, overwrite=True):
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
    gaia = from_fits(ifn_gaia, maxN=20000)
    #print("gaia known cols: ",gaia.known_cols)
    #print("Gaia len",len(gaia), np.shape(gaia.array))
    #enrich_Gaia(gaia)
    
    #ero0 = from_fits(ifn_ero, mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=2000)
    ero0 = from_fits(ifn_ero, maxN=2000)
    #enrich_eROSITA(ero0)
    #print("len 0: ", len(ero0))
          
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
    
    sids = ero1.srcIDs()
    oids = ero1.to_array(colnames="original_srcID", array_type="array")
    #print(sids, oids)
    cnt = np.zeros(len(sids))
    unique_elements, counts_elements = np.unique(oids, return_counts=True)
    for a, b in zip(unique_elements, counts_elements):
        gi = np.where(oids == a)[0]
        #print(a, gi, b)
        cnt[gi] = b
    ero1.add_col("NN_max", cnt.astype(int))
    
        
    #print(ero1.to_array(colnames="srcID_NN"))
    
    to_fits(ero1, ofn, overwrite=True)
    return ero1
    
    print(len(err),np.shape(ide), np.shape(idg), len(ec["offset_sig"]))
    
    if keep_ero_cols is not None:
        for c in np.atleast_1d(keep_ero_cols):
            print(c)
            ec[c] = np.concatenate([ero_ff[1].data[c], ero_ff[1].data[c][gi2], ero_ff[1].data[c][gi3]])
    if keep_gaia_cols is not None:
        for c in np.atleast_1d(keep_ero_cols):
            print(c)
            ec[c] = np.concatenate([gaia_ff[1].data[c], gaia_ff[1].data[c][gi2], gaia_ff[1].data[c][gi3]])                                   
                                   
    
    merge_catalogs(ifn_ero, ide, ifn_gaia, idg, ofn, extra_columns=ec, overwrite=overwrite)
    
    plt.hist(ec["offset_sig"], bins=50, range=(0, 10))
    plt.xlabel("Offset in sigma")
    plt.ylabel("N")
    plt.show()

def training_set(ifn):
    """
    Generate training set based on good positional matches and a criterium on Fx/Fg
    """
    pass

def random_set(ifn):
    """
    Generate a random set, i.e., shuffle the X-ray positions
    """
    pass
