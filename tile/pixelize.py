import healpy as hp
from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import copy

def extract_healpix_range(ifn, ofn=None, min_index=0, max_index=None, colname="healpix", extension=1, overwrite=True, verbose=4):
    """
    Generate a new fits-file containing only sources within specified index range
    
    Parameters
    -----------
    ifn, ofn : str
        Input and output fits-filenames (ofn=ifn if ofn==None)
    extension : int
        Extension in `ifn` to be used
    min_index, max_index : int
        Indices with min_index<=idx<=max_index will be extracted
        If max_index is None, max_index=N+1 with N being the largest index in `ifn`; this is useful to extract, e.g., all indices
    colname : str
        Name for column containing the healpix index
    """
    if ofn is None:
        raise IOError("pixelize::extract_healpix_range - No `ofn` provided.")
    if ofn==ifn and overwrite==False:
        raise IOError("pixelize::extract_healpix_range - `ofn`==`ifn` and `overwrite`==False.")
    
    if verbose>1:
        print("pixelize::extract_healpix_range - Reading ",ifn)
    
    ff = pyfits.open(ifn)

    index = ff[extension].data[colname]
    if max_index is None:
        max_index = max(index)+1
    gi = np.where(np.logical_and(index >= min_index, index <= max_index))[0]
    if verbose>2:
        print("pixelize::extract_healpix_range - ",len(gi), " objects in pixel range.")
    if verbose>3:
        for i in range(min_index, max_index+1):
            ti = np.where(index==i)[0]
            print("pixelize::extract_healpix_range - #entries for healpix %i:%i" % (i, len(ti)))
    
    hdu = copy.copy(ff[0]) 
    cols = ff[extension].columns
    
    ll = []
    for col in ff[extension].columns:
        if verbose>5: print(col, " --- ",col.name, col.format)
        cl = pyfits.Column(name=col.name, array=ff[extension].data[col.name][gi], format=col.format)
        ll.append(cl)
    hdx = pyfits.BinTableHDU.from_columns(ll)
    hdul = pyfits.HDUList([hdu, hdx])
        
    hdul.writeto(ofn, overwrite=overwrite)
    if verbose>0: print('"pixelize::add_healpix_col - Written: ',ofn)
    

def add_healpix_col(ifn, ofn=None, extension=1, overwrite=False, nside=16, verbose=3):
    """
    Add the healpix index based on RA, Dec
    
    Parameters
    ----------
    ifn, ofn : str
        Input and output fits-filenames (ofn=ifn if ofn==None)
    extension : int
        Extension in `ifn` to be used
    nside : int, must be 2**n
        Resolution of healpix map used to calculate the indices.
    """
    if ofn is None:
        ofn = ifn
    if ofn==ifn and overwrite==False:
        raise IOError("pixelize::add_healpix_col - No `ofn` provided and `overwrite`==False.")
    
    if verbose>1:
        print("pixelize::add_healpix_col - Approx. resolution at NSIDE {} is {:.2} deg".format(nside, hp.nside2resol(nside, arcmin=True) / 60))
    
    ff = pyfits.open(ifn)
    theta, phi = copy.deepcopy(ff[extension].data["RA"]), copy.deepcopy(ff[1].data["Dec"])

    theta[theta<0]+=360
    phi+=90

    if verbose>1:
        print("pixelize::add_healpix_col: theta: ",min(theta), max(theta))
        print("pixelize::add_healpix_col: phi: ",min(phi), max(phi))

    x = hp.pixelfunc.ang2pix(nside, (180-phi)/180*np.pi, theta/180*np.pi, nest=False, lonlat=False)
    
    if verbose>1: 
        uni, cnt = np.unique(x, return_counts=True)
        print("pixelize::add_healpix_col: Number of populated healpix pixel: ",len(uni))
    if verbose>2: 
        print("pixelize::add_healpix_col: Populated healpix pixel: ",np.unique(uni))

    
    colnames = [c.name for c in ff[extension].columns]     
    
    outcol = "healpix"
    if outcol in colnames:
        ff[extension].data[outcol] = x
    else:
        hdu = pyfits.PrimaryHDU()    
        col = pyfits.Column(name=outcol, array=x, format="I")
        ff[extension].columns.add_col(col)
        
    hdu = copy.copy(ff[0]) 
    cols = ff[extension].columns
    hdx = pyfits.BinTableHDU.from_columns(cols)
    hdul = pyfits.HDUList([hdu, hdx])
    
    hdul.writeto(ofn, overwrite=overwrite)
    if verbose>0: print('"pixelize::add_healpix_col - Written: ',ofn)


    #m = np.ones(NPIX)
    #m[uni] = cnt

    #hp.mollview(m, norm="log")
    #plt.show()
