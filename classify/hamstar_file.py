from astropy.io import fits as pyfits
import numpy as np

mfn = "../classify/major_proba.fits"
ifn = "../../ero_data/merged_eFEDS.fits"

mff = pyfits.open(mfn)
ff = pyfits.open(ifn)

props = {"ero_ID":"original_srcID","Gaia_ID":"srcID_NN",\
    "match_dist":"match_dist","p_stellar":None,"p_ij":None,\
    "ero_RA":"RA","ero_Dec":"Dec",\
    "Gaia_RA":"RA_NN","Gaia_Dec":"Dec_NN", "det_likeli":"DET_LIKE_0",\
    "rate":None, "rate_error":None,
    "Fx":"Fx","Fx_error":None,"HR":None, "Gmag":"phot_g_mean_mag","BP_RP":"bp_rp","plx":"parallax","plx_err":"parallax_error",\
    "pm_RA":None,"pm_Dec":None}    

cols = {}
for p in props.keys():
    cols[p] = []



gi = np.where(mff[1].data["predicted"] == 0)[0]
for i in gi:
    srcID = mff[1].data["srcID"][i]
    NN = mff[1].data["NN"][i]    
    if len(srcID) > 9:
        osrcID = srcID[0:7]
    else:
        osrcID = srcID
    #print(srcID, osrcID, NN)
    ii = np.where((ff[1].data["original_srcID"] == osrcID) & (ff[1].data["NN"] == NN))[0]
    print(ii)
    for p in props.keys():
        mm = props[p]
        if mm is not None:
            val = ff[1].data[mm][ii]
        else:
            val = np.nan
        if p == "p_stellar":
            val = mff[1].data["proba"][i]
        cols[p].append(val)
    #print()
    
fcols = []    
for p in props.keys():
    fmt = "D"
    if p in ["ero_ID","Gaia_ID"]:
        fmt = "32A"
    if p == "p_stellar2":    
        col = pyfits.Column(name=p, array=len(gi)*[1], format=fmt)        
    elif p == "p_ij":    
        col = pyfits.Column(name=p, array=len(gi)*[1], format=fmt)        
    else:
        col = pyfits.Column(name=p, array=cols[p], format=fmt)        
    fcols.append(col)

cc = pyfits.ColDefs(fcols)    
xx = pyfits.BinTableHDU.from_columns(cc)
hdu = pyfits.PrimaryHDU()    
hdul = pyfits.HDUList([hdu, xx])
hdul.writeto("SVM_new2.fits", overwrite=True)

print(len(np.unique(xx.data["ero_ID"])))
exit()


#[eROSITA source ID]
#{1} ero_ID, string, 32 characters

#[Gaia source ID]
#{2} Gaia_ID, string, 32 characters

#[Match properties]
#{3} match_dist, float (match distance in arcsec)
#{4} p_stellar, float (between 0 and 1)
#{5} p_ij, float (likelihood for this particular association; between 0 
#and 1)

#[eROSITA source position (in degree)]
#{6} ero_RA, float
#{7} ero_Dec, float

#[Gaia source position at epoch of eROISTA observation (in degree)]
#{8} Gaia_RA, float
#{9} Gaia_Dec, float

#[eROSITA X-ray properties]
#{10} det_likeli, float (detection likelihood)
#{11} rate, float ((linear) net count rate in 1/s)
#{12} rate_error, float
#{13} Fx, float, float (linear in erg/s/cm2)
#{14} Fx_error, float
#{15} HR, float

#[Gaia source properties]
#{16} Gmag, float (mag); no erros, because these are small
#{17} BP_RP, float (mag); no erros, because these are small
#{18} plx, float
#{19} plx_err, float
#{20} pm_RA, float (in mas/year)
#{21} pm_Dec, float (in mas/year)

