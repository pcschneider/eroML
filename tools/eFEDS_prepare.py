from astropy.io import fits as pyfits
import numpy as np

ifn = "../../ero_data/eFEDS_c001_clean_V18C.fits"
ofn = "../../ero_data/eFEDS_c001_hs.fits"

ff = pyfits.open(ifn)

gi = np.where((ff[1].data["EXT"]<12) & (ff[1].data["RADEC_ERR"]>0))[0]
print(len(gi), " of ",len(ff[1].data["EXT"]))
cols = []
for c in ff[1].columns:
    if c.name == "detUID":
        arr = ff[1].data[c.name]
        cols.append(
            pyfits.Column(name="DETUID", format=c.format, unit=c.unit, array=arr[gi])
            )
        continue
        
    if c.name == "RADEC_ERR":
        arr = ff[1].data["RADEC_ERR_CORR"][gi]
    else:
        arr = ff[1].data[c.name][gi]
    cols.append(
        pyfits.Column(name=c.name, format=c.format, unit=c.unit, array=arr)
        )

            
new_cols = pyfits.ColDefs(cols)
hdu = pyfits.BinTableHDU.from_columns(new_cols)
hdu.header["TELESCOP"] = "eROSITA"
hdu.writeto(ofn, overwrite=True)
