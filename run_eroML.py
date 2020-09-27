import argparse
#import configparser
from configparser import ConfigParser, ExtendedInterpolation
from eroML.config import *
from eroML.tile import file4
from eroML.tile import merge_fits#, add_healpix_col, hpix2process, generate_healpix_files
#from eroML.utils import download_Gaia_tiles, Gaia_tile_loop
#from eroML.utils import enrich_Gaia, ero_tile_loop
#from eroML.utils import major_loop, random_loop, training_loop
#from eroML.utils import file_loop_1to1, shrink
from eroML.utils import setup_logger
from eroML.tools import calculate_healpix,prepare_Gaia_data, perform_Gaia_download,generate_ero_tiles, perform_ero_data_preparation
from eroML.tools import generate_major_sets, generate_random_sets, generate_training_sets, shrinking
import glob


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
    calculate_healpix(cconfig=config)
    
if config["X data preparation"]["perform"].lower()=="true":
    logger.info("Preparing eROSITA data")
    generate_ero_tiles(cconfig=config)

if config["X data preparation"]["enrich"].lower()=="true":
    logger.debug("Enriching ero-tiles.")
    perform_ero_data_preparation(cconfig=config)

if config["Gaia Download"]["perform"].lower()=="true":
    logger.info("Downloading Gaia data")
    perform_Gaia_download(cconfig=config)    

if config["Enrich Gaia"]["perform"].lower()=="true":
    logger.info("Enriching Gaia data")
    prepare_Gaia_data(cconfig=config)
    

if config["Data sets"]["major"].lower()=="true":
    logger.debug("Generating major sets")
    generate_major_sets(cconfig=config)
        

if config["Data sets"]["random"].lower()=="true":
    logger.debug("Generating random sets")
    generate_random_sets(cconfig=config)


if config["Data sets"]["training"].lower()=="true":
    logger.debug("Generating training sets.")
    generate_training_sets(cconfig=config)
    
    
if config["Merging"]["shrink"].lower()=="true":
    logger.debug("Shrinking files.")    
    merging(cconfig=config)
    

if config["Merging"]["training"].lower() == "true":
    if config["Merging"]["shrink"].lower()=="true":
        which = "training_tiles_small"
    else:
        which = "training_tiles"
    ofn = file4("training")    
    fnames = file4(which, cconfig=config)
    merge_fits(fnames, ofn=ofn)    

    
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


