from astropy.io import fits as pyfits
import glob

gs = "../../ero_data/ero_ID2_nside32_*.fits"

fnames = glob.glob(gs)
#print(len(fnames))

for i,fn in enumerate(fnames):
    print(fn, "(",i+1,"/",len(fnames),")")
    ff = pyfits.open(fn)
    #print(ff[1].columns)
    ff[1].columns.change_name("eligible_eROSITA", "eligible_X")
    ff.writeto(fn, overwrite=True)
    ff.close()
