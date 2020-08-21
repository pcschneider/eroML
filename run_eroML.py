import argparse
#import configparser
from configparser import ConfigParser, ExtendedInterpolation
from eroML.config import *
from eroML.tile import loop, file4, merge_fits
from eroML.tile import merger, add_healpix_col, hpix2process, generate_healpix_files
from eroML.utils import download_Gaia_tiles, Gaia_tile_loop
from eroML.utils import enrich_Gaia, ero_tile_loop
from eroML.utils import major_loop, random_loop, training_loop
from eroML.utils import file_loop_1to1, shrink
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


def setup_logger(main_fn, debug_fn):
    """Setup a Rotating File Handler"""
    formatter = logging.Formatter('%(asctime)s - %(funcName)s - %(levelname)s - %(message)s')

    logger = logging.getLogger("eroML")
    logger.setLevel(4)#logging.DEBUG)
    
    logging.VERBOSE = 5
    logging.addLevelName(logging.VERBOSE, "VERBOSE")
    
    #logging.Logger.verbose = lambda inst, msg, *args, **kwargs: inst.log(logging.VERBOSE, msg, *args, **kwargs)
    #logging.verbose = lambda msg, *args, **kwargs: logging.log(logging.VERBOSE, msg, *args, **kwargs)
    
    handler = logging.handlers.RotatingFileHandler(main_fn, backupCount=5)
    handler.setLevel(level=logging.INFO)
    handler.setFormatter(formatter)
    handler.doRollover()
    logger.addHandler(handler)
    
    #x = hanlder.findCaller()
    
    handler = logging.handlers.RotatingFileHandler(debug_fn, backupCount=5)
    handler.setLevel(logging.DEBUG)#level=logging.DEBUG)
    handler.setFormatter(formatter)
    handler.doRollover()
    logger.addHandler(handler)
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.VERBOSE)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    
    return logger


# file logger
#logger = setup_logger()
#logger.info('Starting eroML')

parser = argparse.ArgumentParser()
parser.add_argument("conf_fn", nargs='?', default=None, help="Config-file")

unlogged=[]

args = parser.parse_args()
if args.conf_fn is not None:
    unlogged.append(("info", "Using custom config-file: `%s` " % args.conf_fn))
    tmp = read_config(args.conf_fn)
    unlogged.append(tmp)

# file logger
logger = setup_logger(config["General"]["main_log_file"], config["General"]["debug_log_file"])
logger.info('Starting eroML')
for m in unlogged:
    out = getattr(logger, m[0])
    out(m[1])

    
if config["Healpix"]["calculate"].lower()=="true":
    logger.info("Adding healpix index to eROSITA source list")
    add_healpix_col(config["Sources"]["ero_filename"], ofn=config["Sources"]["ero_filename_hp"], nside=config["Healpix"].getint("nside"), overwrite=True)


if config["eROSITA preparation"]["perform"].lower()=="true":
    logger.info("Preparing eROSITA data")
    prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    posx = ""
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)    
    generate_healpix_files(config["Sources"]["ero_filename_hp"], prefix=prex, postfix=posx, index0=index0, index1=index1, pix_file=healpix_file)

if config["eROSITA preparation"]["enrich"].lower()=="true":
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    logger.debug("Enriching %i ero-tiles." % len(idx))
    prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    posx = ""
    ero_tile_loop(idx, prefix=prex, postfix=posx)

if config["Gaia Download"]["perform"].lower()=="true":
    logger.info("Downloading Gaia data")
    
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)
    
    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    logger.debug("Downloading %i Gaia-tiles." % len(idx))
    
    download_Gaia_tiles(outdir=config["Gaia Download"]["directory"],\
               prefix=config["Gaia Download"]["prefix"], idx=idx,\
               nside=int(config["Healpix"]["nside"]),\
               overwrite=config["Gaia Download"]["overwrite"],\
               edge=float(config["Gaia Download"]["edge"]),\
               verbose=int(config["Gaia Download"]["verbose"]))


if config["Enrich Gaia"]["perform"].lower()=="true":
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    logger.debug("Enriching %i ero-tiles." % len(idx))
    prex = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    posx = ""
    
    filt = int(config["Enrich Gaia"]["filter"])
    #print("filt: ",filt)
    logger.info("Enriching Gaia data files (filt=%i)" % filt)
    Gaia_tile_loop(idx, prefix=prex, postfix=posx, filterNr=filt)
    

if config["Data sets"]["major"].lower()=="true":
    
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    g_prex = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    g_posx = ""
    
    e_prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    e_posx = ""
    
    m_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    m_posx = ""
    
    logger.debug("Major sets for %i tiles." % len(idx))
    major_loop(idx, ero_prefix=e_prex, ero_postfix=e_posx, gaia_prefix=g_prex, gaia_postfix=g_posx, major_prefix=m_prex, major_postfix=m_posx)
    
        


if config["Data sets"]["random"].lower()=="true":
    
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    g_prex = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    g_posx = ""
    
    e_prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    e_posx = ""
    
    r_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    r_posx = ""
    
    mino = float(config["Data sets"]["min_random_offset"])
    maxo = float(config["Data sets"]["max_random_offset"])
    multi= int(config["Data sets"]["random_multi"])
    logger.debug("Random data sets for %i tiles." % len(idx))
    random_loop(idx, ero_prefix=e_prex, ero_postfix=e_posx, gaia_prefix=g_prex, gaia_postfix=g_posx, random_prefix=r_prex, random_postfix=r_posx, min_offset=mino, max_offset=maxo, multi=multi)
    




if config["Data sets"]["training"].lower()=="true":
    
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    
    # random
    r_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    r_posx = ""
    
    #major
    m_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    m_posx = ""
    
    #training
    t_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    t_posx = ""
    
    
    ad = float(config["Data sets"]["training_abs_dist"])
    rd = float(config["Data sets"]["training_rel_dist"])
    
    logger.debug("Random data sets for %i tiles." % len(idx))
    training_loop(idx, major_prefix=m_prex, major_postfix=m_posx, random_prefix=r_prex, random_postfix=r_posx, training_prefix=t_prex, training_postfix=t_posx, abs_dist=ad, rel_dist=rd)

if config["Merging"]["shrink"].lower()=="true":
    healpix_file = config["Healpix"].get("pix_file", None)
    index0 = config["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    
    r_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    r_posx = ""
    
    #major
    m_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    m_posx = ""
    
    #training
    t_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    t_posx = ""
 
    sh = config["Merging"]["shrink_postfix"]
    cols = config["Columns"]["keep"].split(",")
    #print("cols: ",cols)
    
    for tt in ["random","training","major"]:
        prex = config["Data sets"]["directory"]+"/"+config["Data sets"][tt+"_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
        posx = ""
        file_loop_1to1(idx, prefix=prex, postfix=posx, method=shrink, ofn_prefix=prex, ofn_postfix=posx+sh,cols=cols)
    
    
if config["Merging"]["major"].lower() == "true":
    if config["Merging"]["shrink"].lower()=="true":
        which = "major_tiles_small"
    else:
        which = "major_tiles"
    ofn = file4("major")    
    fnames = file4(which)
    merge_fits(fnames, ofn=ofn)

if config["Merging"]["random"].lower() == "true":
    if config["Merging"]["shrink"].lower()=="true":
        which = "random_tiles_small"
    else:
        which = "random_tiles"
    ofn = file4("random")    
    fnames = file4(which)
    merge_fits(fnames, ofn=ofn)
    
if config["Merging"]["training"].lower() == "true":
    if config["Merging"]["shrink"].lower()=="true":
        which = "training_tiles_small"
    else:
        which = "training_tiles"
    ofn = file4("training")    
    fnames = file4(which)
    merge_fits(fnames, ofn=ofn)    

    
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


