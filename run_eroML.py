import argparse
import configparser
from eroML.config import *
from eroML.tile import loop
from eroML.tile import merger
from eroML.utils import download_Gaia_tiles
from eroML.utils import enrich_Gaia
import logging.handlers
import logging
import glob

def discover_filenames(which):
    """
    """
    if which.lower() == "gaia":
        glob_str = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_*.fits"
        fnames = glob.glob(glob_str)
    return fnames


def setup_logger():
    """Setup a Rotating File Handler"""
    formatter = logging.Formatter('%(asctime)s - %(funcName)s - %(levelname)s - %(message)s')

    logger = logging.getLogger("eroML")
    logger.setLevel(4)#logging.DEBUG)
    
    logging.VERBOSE = 5
    logging.addLevelName(logging.VERBOSE, "VERBOSE")
    
    #logging.Logger.verbose = lambda inst, msg, *args, **kwargs: inst.log(logging.VERBOSE, msg, *args, **kwargs)
    #logging.verbose = lambda msg, *args, **kwargs: logging.log(logging.VERBOSE, msg, *args, **kwargs)
    
    handler = logging.handlers.RotatingFileHandler("main.log", backupCount=5)
    handler.setLevel(level=logging.INFO)
    handler.setFormatter(formatter)
    handler.doRollover()
    logger.addHandler(handler)
    
    #x = hanlder.findCaller()
    
    handler = logging.handlers.RotatingFileHandler("debug.log", backupCount=5)
    handler.setLevel(logging.DEBUG)#level=logging.DEBUG)
    handler.setFormatter(formatter)
    handler.doRollover()
    logger.addHandler(handler)
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.VERBOSE)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    
    return logger



def read_config(fn):
    cc = configparser.ConfigParser()
    cc.read(fn)
    logger.info("Reading config-file: ",fn)
    for sec in cc.sections():
        #print("  [",sec,"]")
        for k in cc[sec].keys():
            #print("     ",k, ": ", cc[sec][k])
            config[sec][k] = cc[sec][k]
    #print()

# file logger
logger = setup_logger()
logger.info('Starting eroML')

parser = argparse.ArgumentParser()
parser.add_argument("conf_fn", nargs='?', default=None, help="Config-file")

args = parser.parse_args()
if args.conf_fn is not None:
    logger.info("Using config-file: `%s` " % args.conf_fn)
    read_config(args.conf_fn)


if config["Gaia Download"]["perform"].lower()=="true":
    logger.info("Downloading Gaia data")
    download_Gaia_tiles(outdir=config["Gaia Download"]["directory"],\
               prefix=config["Gaia Download"]["prefix"], idx=None,\
               nside=int(config["Healpix"]["nside"]),\
               overwrite=config["Gaia Download"]["overwrite"],\
               edge=float(config["Gaia Download"]["edge"]),\
               verbose=int(config["Gaia Download"]["verbose"]))
 
if config["Enrich Gaia"]["perform"].lower()=="true":
    fnames = discover_filenames("gaia")
    logger.info("Enriching Gaia data files")
    for fn in fnames:
        logger.debug("Working on %s" % fn)
        enrich_Gaia(fn)
    
exit()    

#loop(config["Sources"]["ero_fn"],ofn=config["Healpix"]["ero_fn_hp"], NSIDE=int(config["Healpix"]["hp_nside"]), rID=config["Sources"]["rID"])

# Merge
# Training:
glob_str = config["Sources"]["data_dir"] +"/*rID"+config["Sources"]["rID"]+"_training_random.fits"
print("Merging ",glob_str," -> ",config["Sources"]["training_fn"])
merger.merge_matching(glob_str,ofn=config["Sources"]["training_fn"], overwrite=True)

## Random
#glob_str = config["Sources"]["data_dir"] +"/*rID"+config["Sources"]["rID"]+"_random.fits"
#print("Merging ",glob_str," -> ",config["Sources"]["random_fn"])
#merger.merge_matching(glob_str,ofn=config["Sources"]["random_fn"], overwrite=True)



# Major
#glob_str = config["Sources"]["data_dir"] +"/*rID"+config["Sources"]["rID"]+"_major.fits"
#print("Merging ",glob_str)
#merger.merge_matching(glob_str,ofn=config["Sources"]["merged_fn"], overwrite=True)


