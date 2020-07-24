import healpy as hp
from astropy.io import fits as pyfits

NSIDE = 16
print(
    "Approximate resolution at NSIDE {} is {:.2} deg".format(
        NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
    )
)
    
NPIX = hp.nside2npix(NSIDE)
print(NPIX)

    
fn = "../../eFEDS/SrcCat_V2T.fits"
ff = pyfits.open(fn)
