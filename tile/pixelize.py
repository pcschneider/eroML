import healpy as hp
from astropy.io import fits as pyfits
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass    
import copy
import logging
import os

logger = logging.getLogger('eroML')

def populated_hpix(fn, colname="healpix", extension=1):
    """
    The populated healpix 
    """
    logger.debug("Reading: %s" % fn)
    ff = pyfits.open(fn)
    index = ff[extension].data[colname]    
    #hpix_range=min(index), max(index)
    ff.close()
    return np.unique(index)
 
 
def hpix2process(ifn, index0=0, index1=None, pix_file=None, colname="healpix", extension=1):
    """
    Returns list of healpix indices based on `index0`, `index1`, and `pix_file` for the healpix indices that are in `ifn`
    """
    hpix_all = populated_hpix(ifn, colname=colname, extension=extension)
    hpix_range=min(hpix_all), max(hpix_all)    

    logger.debug("Available HPIX range: %i - %i" % (hpix_range[0], hpix_range[1]))
    
    if pix_file is not None:
        logger.info("Reading healpix from '%s'" % pix_file)
        dd = np.genfromtxt(pix_file, dtype=int)
        indices = np.atleast_1d(dd)
    else:
        indices = np.arange(hpix_range[0], hpix_range[1])
        
    if index1==None: index1=max(indices)
    
    if index1<hpix_range[0]:
        logger.error("Cannot extract HPIX with index1=%i, because hpix range is  %i - %i" % (index1, hpix_range[0], hpix_range[1]))
        raise IndexError(str("Cannot extract HPIX with index1=%i, because hpix range is  %i - %i" % (index1, hpix_range[0], hpix_range[1])))
    
    if index0<hpix_range[0]: 
        index0=hpix_range[0]
        logger.warning("Adjusting lower bound of requested index range: %i - %i" % (index0, index1))

    if index1>hpix_range[1]: 
        index1=hpix_range[1]
        logger.warning("Adjusting upper bound of requested index range: %i - %i" % (index0, index1))

    gi = np.where((indices>=index0) & (indices <= index1) )[0]
    return indices[gi]
    
    
def generate_healpix_files(ifn, index0=0, index1=None, pix_file=None, prefix="", postfix="", colname="healpix", extension=1, overwrite=True, skip=True, verbose=4):
    """
    Split the full source list into sub-lists containing sources only for individual healpix
    
    Parameters
    -----------
    ifn : str
        The eROSITA source file
    prefix, postfix : str
        The resulting filenames will be `prefix`index`postfix`.fits
    index0, index1 : int 
        min and max indices to extract
    pix_file : str
        Extract only healpix that are in this file
    colname : str
        The column name containing the healpix index for each source
    extension : int
        The extension in the fits-file (in `ifn`)
    skip : boolean
        Skip file creation if file already exists.
    """
    idx = hpix2process(ifn, index0=index0, index1=index1, pix_file=pix_file, colname=colname, extension=extension)

    for j, i in enumerate(idx):
        ofn=prefix+str(i)+postfix+'.fits'
        if os.path.exists(ofn):
            logger.warning("Skipping %s as it already exists (and `skip`==True)..." % ofn)
            continue
        logger.info("Generating %s (file %i/%i)" % (ofn, j+1, len(idx)))
        extract_healpix_range(ifn, ofn, min_index=i, max_index=i, colname=colname, extension=extension, overwrite=overwrite)
        

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
    
    logger.info("Reading %s" % ifn)
    
    ff = pyfits.open(ifn)

    index = ff[extension].data[colname]
    if max_index is None:
        max_index = max(index)+1
    gi = np.where(np.logical_and(index >= min_index, index <= max_index))[0]
    
    logger.debug("     %i objects in pixel range." % len(gi))
    #if verbose>3:
        #for i in range(min_index, max_index+1):
            #ti = np.where(index==i)[0]
            #print("pixelize::extract_healpix_range - #entries for healpix %i:%i" % (i, len(ti)))
    
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
    logger.info("Written %s with %i entries." % (ofn, len(gi)))
    

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
    else:
        logger.debug("Using ofn=%s" % ofn)
        
    if verbose>1:
        logger.debug("pixelize::add_healpix_col - Approx. resolution at NSIDE {} is {:.2} deg".format(nside, hp.nside2resol(nside, arcmin=True) / 60))
    
    ff = pyfits.open(ifn)
    ra, dec = copy.deepcopy(ff[extension].data["RA_CORR"]), copy.deepcopy(ff[1].data["DEC_CORR"])

    ra[ra<0]+=360
    phi = 90-dec
    theta = ra

    #resol = hp.nside2resol(nside, arcmin=True)
    #pix_center = hp.pix2ang(nside, hpix, nest=True)
    #ra, dec = pix_center[1]/np.pi*180, 90. - pix_center[0]/np.pi*180

    if verbose>1:
        logger.debug("pixelize::add_healpix_col: theta: %f %f " % (min(theta), max(theta)))
        logger.debug("pixelize::add_healpix_col: phi: %f %f " % (min(phi), max(phi)))

    x = hp.pixelfunc.ang2pix(nside, phi/180*np.pi, theta/180*np.pi, nest=True, lonlat=False)
    
    uni, cnt = np.unique(x, return_counts=True)    
    if verbose>1: 
        logger.debug("pixelize::add_healpix_col: Number of populated healpix pixel: %i " % len(uni))
    if verbose>2: 
        logger.log(5, "pixelize::add_healpix_col: Populated healpix pixel: %s " % str(np.unique(uni)))

    
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
    if verbose>0: logger.info('"pixelize::add_healpix_col - Written: %s ' % ofn)

    return uni

    #m = np.ones(NPIX)
    #m[uni] = cnt

    #hp.mollview(m, norm="log")
    #plt.show()
