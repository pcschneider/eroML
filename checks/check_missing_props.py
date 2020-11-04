import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u

missing_ids = np.genfromtxt("missing_srcIDs.txt", dtype=str)
fn0 = "../../ero_data/SrcCat_V2T_hp.fits"
fn1 = "../../ero_data/major_eFEDS.fits"

ff0 = pyfits.open(fn0)
ff1 = pyfits.open(fn1)

srcIDs0 = ff0[1].data["srcID"]
srcIDs1 = ff1[1].data["srcID"]

columns = ["ID_SRC","RA","Dec", "QUALITY", "RADEC_ERR", "EXT", "EXT_LIKE", "ML_FLUX_0", "ML_FLUX_ERR_0", "DET_LIKE_0", "Fx", "healpix"]
for i in missing_ids:#[0:2]:
    gi = np.where(srcIDs0 == i)[0]
    if len(gi) == 1: gi=gi[0]
    else:
        print("Argh")
        continue
    print(i, gi)
    oo = []
    for c in columns:
        oo.append("%s: %s " % (c, str(ff0[1].data[c][gi])))
    print("; ".join(oo))
    gaia_glob = str("../../ero_data/Gaia_ID2_nside32_%i.fits" % ff0[1].data["healpix"][gi])
    ra, dec = ff0[1].data["RA"][gi], ff0[1].data["Dec"][gi]
    co = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
    print(gaia_glob)
    gfn = glob.glob(gaia_glob)    
    if len(gfn)==1:
        ffg = pyfits.open(gfn[0])
        rag, decg = ffg[1].data["RA"], ffg[1].data["Dec"]
        cog = SkyCoord(ra=rag, dec=decg, unit=(u.degree, u.degree))
        idx, d2d, d3d = co.match_to_catalog_sky(cog)
        print("Nearest : ",idx, d2d.arcsec)
        for cc in ["eligible_Gaia","Gaia_quality", "iso_compatible", "phot_g_mean_mag","parallax", "parallax_error"]:
            print("Nearest : ",cc,"-",ffg[1].data[cc][idx])
        ei = np.where(ffg[1].data["eligible_Gaia"] == 1)[0]
        idx, d2d, d3d = co.match_to_catalog_sky(cog[ei])
        print("Eligible: ",idx, d2d.arcsec)
        for cc in ["eligible_Gaia","Gaia_quality", "iso_compatible", "phot_g_mean_mag","parallax", "parallax_error"]:
            print("Eligible: ",cc,"-",ffg[1].data[cc][ei[idx]])
    else:
        print("Argh2")
        continue
    print()
