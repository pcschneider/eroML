from eroML.ensemble import Ensemble
from eroML.ensemble import from_fits,to_fits
import numpy as np

def enrich_Gaia(e):
    arr = e.to_array(colnames=["phot_g_mean_mag"])
    FG = 10**(-0.4* arr["phot_g_mean_mag"])*3.660e-08*720
    e.add_col("F_G", FG)

def enrich_eROSITA(e):
    err = e.to_array(colnames="RADEC_ERR", array_type="array")
    err[err<1.] = 1.
    e.set_col("RADEC_ERR", err)


def enrich_merged(e):
    offset_sig = e.to_array(colnames="RADEC_ERR", array_type="array")
    d2d = e.to_array(colnames="match_dist", array_type="array")
    e.add_col("offset_sig", d2d/offset_sig)
    #ec={"offset_sig": np.concatenate([d2d.arcsec/err, d2d2[gi2].arcsec/err[gi2], d2d3[gi3].arcsec/err[gi3]])}


def major_catalog(ifn_gaia, ifn_ero, ofn, keep_ero_cols=None, keep_gaia_cols=None, overwrite=True):
    """
      Parameters
      ----------
      keep_[ero/gaia] cols : refers to ifn
    """


    gaia = from_fits(ifn_gaia, mapper={"source_id":"srcID", "ra":"RA", "dec":"Dec"})#, maxN=2000)
    enrich_Gaia(gaia)
    
    ero = from_fits(ifn_ero, mapper={"detUID":"srcID", "DEC":"Dec"})#, maxN=2000)
    enrich_eROSITA(ero)
    
    ero.merge_add(gaia, NN=1)
    #gi2 = np.where(d2d2.arcsec < 5*err)[0]
    #gi3 = np.where(d2d3.arcsec < 5*err)[0]

    enrich_merged(ero)
    
    offset_sig = ero.to_array(colnames="offset_sig", array_type="array")
    d2d = ero.to_array(colnames="match_dist", array_type="array")
    gi = np.where((offset_sig<50) & (d2d<1000))[0]
    srcIDs = np.array(ero.srcIDs())
    #print(len(gi), gi, srcIDs )
    good_ids = srcIDs[gi]
    print("Keeping: ",good_ids)
    ero.keep(good_ids)
    
    to_fits(ero, "test.fits", overwrite=True)
    return
    
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

class Tile():
    """
    Orchestrates the analyzes of a sky region.
    """
    pass
