from astropy.io import fits as pyfits

fn = "../eFEDS/SrcCat_V2T.fits"
ff = pyfits.open(fn)
    
cols = [] 
cols.append( pyfits.Column(name='RA_CORR', format='D', array=ff[1].data["RA"]) )
cols.append( pyfits.Column(name='DEC_CORR', format='D', array=ff[1].data["DEC"]) )

orig_cols = ff[1].columns
new_cols = pyfits.ColDefs(cols)
hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)

ff[1] = hdu
ff[1].header["TELESCOP"] = "eROSITA"

ff.writeto("../eFEDS/SrcCat_V2T_q.fits", overwrite=True)
