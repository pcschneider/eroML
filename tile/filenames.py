#from eroML.config import *
from .pixelize import hpix2process
import logging

logger = logging.getLogger('eroML')

def discover_filenames(which, config=None):
    """
    Discover the EXISITING filenames for a certain data set
    
    Options for `which` are: 
      - gaia_tiles
      - ero_tiles
      
    Returns
    -------
    filenames : array of str
    """
    if config is None:
        from eroML.config import config
    
    if which.lower() == "gaia_tiles":
        glob_str = config["Gaia_Download"]["directory"]+"/Gaia_"+config["Gaia_Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_*.fits"
        fnames = glob.glob(glob_str)

    elif which.lower() == "ero_tiles":
        glob_str = config["Xdata_preparation"]["directory"]+"/"+config["Xdata_preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_*.fits"
        fnames = glob.glob(glob_str) 
        
    return fnames



def file4(which, cconfig=None):# which=""):
    """
    Options for `which`:
      - X_filename
      - X_filename_hp : File containing ALL eROSITA sources with healpix indices
      - gaia_tiles
      - ero_files
      - major_tiles (plus _small)
      - random_tiles  (plus _small)
      - training_tiles  (plus _small)
      - major 
      - training
      - random
    
    Example::
    
        from eroML.tile import file4
        fn = file4("major", cconfig="eFEDS_EDR3.ini")
        print(fn)
    
    Parameters
    ----------
    which : str (see above)
    cconfig : str
        Filename for config file
    
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
        from eroML.config import config, read_config
        logger.info("Using custom config-file: `%s` " % cconfig)
        tmp = read_config(cconfig)
        print(tmp)
    else:
        from eroML.config import config
        print("XXX")        
    
    idx = None
    
    if "tiles" in which:    
        healpix_file = config["Healpix"].get("pix_file", None)
        index0 = config["Healpix"].getint("index0", 0)
        index1 = config["Healpix"].getint("index1", None)
        idx = hpix2process(config["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
        logger.debug("Creating filenames for %i tiles." % len(idx))

    if which=="X_filename":
        return config["Sources"]["X_filename"]
    
    elif which=="X_filename_hp":
        return config["Sources"]["X_filename_hp"]
    
    elif which=="gaia_tiles":
        prex = config["Gaia_Download"]["directory"]+"/"+config["Gaia_Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="ero_tiles":
        prex = config["Xdata_preparation"]["directory"]+"/"+config["Xdata_preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
             
    elif which=="major_tiles":
        prex = config["Datasets"]["directory"]+"/"+config["Datasets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="random_tiles":
        prex = config["Datasets"]["directory"]+"/"+config["Datasets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="training_tiles":
        prex = config["Datasets"]["directory"]+"/"+config["Datasets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        return ff4idx(prex, posx)
    
    elif which=="major_tiles_small":
        prex = config["Datasets"]["directory"]+"/"+config["Datasets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = config["Merging"]["shrink_postfix"]
        return ff4idx(prex, posx)
    
    elif which=="random_tiles_small":
        prex = config["Datasets"]["directory"]+"/"+config["Datasets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = config["Merging"]["shrink_postfix"]
        return ff4idx(prex, posx)
    
    elif which=="training_tiles_small":
        prex = config["Datasets"]["directory"]+"/"+config["Datasets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
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
        
    
