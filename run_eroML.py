import argparse
#import configparser
from configparser import ConfigParser, ExtendedInterpolation
from eroML.config import *
#from eroML.tile import loop, file4, merge_fits
#from eroML.tile import merger, add_healpix_col, hpix2process, generate_healpix_files
#from eroML.utils import download_Gaia_tiles, Gaia_tile_loop
#from eroML.utils import enrich_Gaia, ero_tile_loop
#from eroML.utils import major_loop, random_loop, training_loop
#from eroML.utils import file_loop_1to1, shrink
from eroML.utils import setup_logger
from eroML.tools import calculate_healpix,prepare_Gaia_data, perform_Gaia_download
from eroML.tools import generate_ero_tiles, perform_ero_data_preparation
from eroML.tools import create_major_sets, create_random_sets, create_training_sets
from eroML.tools import fake_positions
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
    
if config["eROSITA preparation"]["perform"].lower()=="true":
    logger.info("Preparing eROSITA data")
    generate_ero_tiles(cconfig=config)

if config["eROSITA preparation"]["enrich"].lower()=="true":
    logger.debug("Enriching ero-tiles.")
    perform_ero_data_preparation(cconfig=config)

if config["Gaia Download"]["perform"].lower()=="true":
    logger.info("Downloading Gaia data")
    perform_Gaia_download(cconfig=config)    

if config["Enrich Gaia"]["perform"].lower()=="true":
    logger.info("Enriching Gaia data")
    prepare_Gaia_data(cconfig=config)

if config["Data sets"]["major"].lower()=="true":    
    logger.debug("Creating major sets.")
    create_major_sets(cconfig=config)
    
if config["Data sets"]["random"].lower()=="true":
    logger.debug("Creating random sets.")
    create_random_sets(cconfig=config)

if config["Data sets"]["training"].lower()=="true":
    logger.debug("Creating training sets.")
    create_training_sets(cconfig=config)
    
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

 
if config["Data preparation"]["fake_positions"].lower() == "true":
    fake_positions(cconfig=config)
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


