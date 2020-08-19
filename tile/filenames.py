from eroML.config import *
from .pixelize import hpix2process
import logging

logger = logging.getLogger('eroML')

def file4(which, cconfig=None):# which=""):
    """
    Options for `which`:
      - ero_filename_hp : File containing ALL eROSITA sources with healpix indices
    
    Parameters
    ----------
    
    Returns  
    -------
    filename(s) : str or array of str
    """    
    def ff4idx(prex, posx):
        rr = [] 
        for i in idx:
            fn = prex+str(i)+posx+'.fits'
            rr.append(fn)
        return rr    
    
    if cconfig is not None:
        logger.info("Using custom config-file: `%s` " % cconfig)
        read_config(cconfig)
    
    idx = None
    
    if "tiles" in which:    
        healpix_file = config["Healpix"].get("pix_file", None)
        index0 = config["Healpix"].getint("index0", 0)
        index1 = config["Healpix"].getint("index1", None)
        idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
        logger.debug("Creating filenames for %i tiles." % len(idx))

    if which=="ero_filename":
        return config["Sources"]["ero_filename"]
    
    elif which=="ero_filename_hp":
        return config["Sources"]["ero_filename_hp"]
    
    elif which=="gaia_tiles":
        prex = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="ero_tiles":
        prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
             
    elif which=="major_tiles":
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="random_tiles":
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="training_tiles":
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="major_tiles_small":
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = config["Merging"]["shrink_postfix"]
        return ff4idx(prex, posx)
    
    elif which=="random_tiles_small":
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = config["Merging"]["shrink_postfix"]
        return ff4idx(prex, posx)
    
    elif which=="training_tiles_small":
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = config["Merging"]["shrink_postfix"]
        return ff4idx(prex, posx)
    
    elif which=="major":
        return config["Sources"]["major_filename"]
    
    elif which=="training":
        return config["Sources"]["training_filename"]
    
    elif which=="random":
        return config["Sources"]["random_filename"]
    
    else:
        logger.error("I don't know anything about %s!" % str(which))
        raise ValueError(str(" 'which'==%s not known..." % which))
        
    
