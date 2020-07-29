from eroML.ensemble import Ensemble
from eroML.ensemble import from_fits,to_fits,multi_fits_support
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import copy
from .gaia_tools import gaia4ero
from .estimators import NN_distribution
from .enrich import enrich_merged


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
    return dy>log_margin
        



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
    
    #to_fits(ero1, ofn, overwrite=True)
    return ero1
    
    #print(len(err),np.shape(ide), np.shape(idg), len(ec["offset_sig"]))
    
    #if keep_ero_cols is not None:
        #for c in np.atleast_1d(keep_ero_cols):
            #print(c)
            #ec[c] = np.concatenate([ero_ff[1].data[c], ero_ff[1].data[c][gi2], ero_ff[1].data[c][gi3]])
    #if keep_gaia_cols is not None:
        #for c in np.atleast_1d(keep_ero_cols):
            #print(c)
            #ec[c] = np.concatenate([gaia_ff[1].data[c], gaia_ff[1].data[c][gi2], gaia_ff[1].data[c][gi3]])                                   
                                   
    
    #merge_catalogs(ifn_ero, ide, ifn_gaia, idg, ofn, extra_columns=ec, overwrite=overwrite)
    
    #plt.hist(ec["offset_sig"], bins=50, range=(0, 10))
    #plt.xlabel("Offset in sigma")
    #plt.ylabel("N")
    #plt.show()

def training_set(ifn):
    """
    Generate training set based on good positional matches and a criterium on Fx/Fg
    """
    ff_stellar = pyfits.open(stellar_fn)
    ff_non_stellar = pyfits.open(non_stellar_fn)
    ff_random = pyfits.open(random_fn)
    N_stellar = len(ff_stellar[1].data["id_gaia"])
    N_non_stellar = len(ff_non_stellar[1].data["id_gaia"])
    N_random = len(ff_random[1].data["id_gaia"])
    print("%s:%i, %s:%i, %s:%i" % (stellar_fn,N_stellar, non_stellar_fn, N_non_stellar, random_fn, N_random))
    hdu = pyfits.PrimaryHDU()    
    cols = []
    for c in ff_stellar[1].columns:
        col = pyfits.Column(name=c.name, array=np.concatenate((ff_stellar[1].data[c.name], ff_non_stellar[1].data[c.name], ff_random[1].data[c.name])), format=c.format)
        cols.append(col)
    
    category = np.zeros(N_stellar+N_non_stellar+N_random)
    category[0:N_stellar] = 1
    ccol = pyfits.Column(name="class", array=category, format="J")
    cols.append(ccol)
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    xx.header["stellar_file"] = stellar_fn
    xx.header["non_stellar_file"] = non_stellar_fn
    xx.header["random_file"] = random_fn
    hdul = pyfits.HDUList([hdu, xx])
    hdul.writeto(ofn, overwrite=overwrite)     
    print("Written: ",ofn, " with ", len(category), " sources.")

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
    print(len(coords))
    N = len(coords)
    angle = np.random.rand(N)*360
    offset = np.random.rand(N)*(max_offset- min_offset)+min_offset
    print(offset)

    doords = coords.directional_offset_by(angle*u.degree, offset*u.arcsec)

    nra, ndec = doords.ra.degree, doords.dec.degree
    ero.set_col("RA", nra)
    ero.set_col("Dec", ndec)
    
    return major_set(ero, gaia)
    

    idx, d2d,d3d = coords.match_to_catalog_sky(doords)
    print(d2d.arcsec)

    hdu = pyfits.PrimaryHDU()    
    cols = []
    for c in ff[1].columns:
        print(c.name)
        if c.name == "RA":
            col = pyfits.Column(name=c.name, array=doords.ra.degree, format=c.format)
        elif c.name == "RA":
            col = pyfits.Column(name=c.name, array=doords.dec.degree, format=c.format)
        else:
            col = pyfits.Column(name=c.name, array=ff[1].data[c.name], format=c.format)
        cols.append(col)
    cc = pyfits.ColDefs(cols)
    xx = pyfits.BinTableHDU.from_columns(cc)
    xx.header["min_offset"] = min_offset
    xx.header["max_offset"] = max_offset
    hdul = pyfits.HDUList([hdu, xx])
    hdul.writeto(ofn, overwrite=overwrite)
