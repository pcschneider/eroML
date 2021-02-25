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
from eroML.tools import calculate_healpix,prepare_Gaia_data, perform_Gaia_download,calculate_Gaia_sky_density
from eroML.tools import generate_ero_tiles, perform_ero_data_preparation
from eroML.tools import create_major_sets, create_random_sets, create_training_sets
from eroML.tools import fake_positions, generate_major_sets, shrinking
import glob
from classify import prepare_classify

parser = argparse.ArgumentParser()
parser.add_argument("conf_fn", nargs='?', default=None, help="Config-file")
parser.add_argument("--steps", dest='steps', nargs='*', action='append', default=None, help="Steps to process")

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

if args.steps is not None and "Datasets" in args.steps[0]:
    args.steps[0].remove("Datasets")
    args.steps[0].append("Datasets/major")
    args.steps[0].append("Datasets/random")
    args.steps[0].append("Datasets/training")

if args.steps is not None and "Merging" in args.steps[0]:
    args.steps[0].remove("Merging")
    args.steps[0].append("Merging/major")
    args.steps[0].append("Merging/random")
    args.steps[0].append("Merging/training")
    

if args.steps is not None:
    logger.info("Performing only requested processing steps: ")
    for i, s in enumerate(args.steps[0]):
        logger.info("      (%i)  %s" % (i+1, s))
    logger.info("Other processing steps selected in the ini-file are ignored.")    
        
if (config["Healpix"]["calculate"].lower()=="true" and args.steps==None) or (args.steps and "Healpix/calculate" in args.steps[0]):
    logger.info("Adding healpix index to X-ray source list")
    calculate_healpix(cconfig=config)
    
if (config["Xdata_preparation"]["perform"].lower()=="true" and args.steps==None) or (args.steps and "Xdata_preparation/perform" in args.steps[0]):
    logger.info("Preparing X-ray data")
    generate_ero_tiles(cconfig=config)

if (config["Xdata_preparation"]["enrich"].lower()=="true" and args.steps==None) or (args.steps and "Xdata_preparation/enrich" in args.steps[0]):
    logger.debug("Enriching ero-tiles.")
    perform_ero_data_preparation(cconfig=config)

if (config["Gaia_Download"]["perform"].lower()=="true" and args.steps==None) or (args.steps and "Gaia_Download/perform" in args.steps[0]):
    logger.info("Downloading Gaia data")
    perform_Gaia_download(cconfig=config)    

if (config["Enrich_Gaia"]["perform"].lower()=="true" and args.steps==None) or (args.steps and "Enrich_Gaia/perform" in args.steps[0]):
    logger.info("Enriching Gaia data")
    prepare_Gaia_data(cconfig=config)

if (config["Enrich_Gaia"]["calculate_density"].lower()=="true" and args.steps==None) or (args.steps and "Enrich_Gaia/calculate_density" in args.steps[0]):
    logger.info("Calculating sky density for (eligible) Gaia sources")
    calculate_Gaia_sky_density(cconfig=config)


if (config["Datasets"]["major"].lower()=="true" and args.steps==None) or (args.steps and "Datasets/major" in args.steps[0]):
    logger.debug("Generating major sets")
    generate_major_sets(cconfig=config)

if (config["Datasets"]["random"].lower()=="true" and args.steps==None) or (args.steps and "Datasets/random" in args.steps[0]):
    logger.debug("Generating random sets")
    generate_random_sets(cconfig=config)

if (config["Datasets"]["training"].lower()=="true" and args.steps==None) or (args.steps and "Datasets/training" in args.steps[0]):
    logger.debug("Generating training sets.")
    generate_training_sets(cconfig=config)

    
if (config["Merging"]["shrink"].lower()=="true" and args.steps==None) or (args.steps and "Merging/shrink" in args.steps[0]):
    logger.debug("Shrinking files.")    
    shrinking(cconfig=config)
    
if (config["Merging"]["training"].lower() == "true" and args.steps==None) or (args.steps and "Merging/training" in args.steps[0]):
    if config["Merging"]["shrink"].lower()=="true" or "Merging/shrink" in args.steps[0]:
        which = "training_tiles_small"
    else:
        which = "training_tiles"
    ofn = file4("training")    
    fnames = file4(which, cconfig=config)
    merge_fits(fnames, ofn=ofn)    
    
if (config["Merging"]["major"].lower() == "true" and args.steps==None) or (args.steps and "Merging/major" in args.steps[0]):
    if config["Merging"]["shrink"].lower()=="true" or "Merging/shrink" in args.steps[0]:
        which = "major_tiles_small"
    else:
        which = "major_tiles"
    ofn = file4("major")    
    fnames = file4(which)
    merge_fits(fnames, ofn=ofn)

if (config["Merging"]["random"].lower() == "true" and args.steps==None) or (args.steps and "Merging/random" in args.steps[0]):
    if config["Merging"]["shrink"].lower()=="true" or "Merging/shrink" in args.steps[0]:
        which = "random_tiles_small"
    else:
        which = "random_tiles"
    ofn = file4("random")    
    fnames = file4(which)
    merge_fits(fnames, ofn=ofn)
    
    
if (config["Classification"]["prepare"].lower() == "true" and args.steps==None) or (args.steps and "Classification/prepare" in args.steps[0]):
    
    logger.info("Preparing major file for classification...")
    ifn = file4("major")
    ofn = config["Classification"]["major_filename"]
    ovwr = config["Classification"]["overwrite"].lower() == "true"
    prepare_classify(ifn, ofn=ofn, overwrite=ovwr)
    
    logger.info("Preparing training file for classification...")
    ifn = file4("training")
    ofn = config["Classification"]["training_filename"]
    ovwr = config["Classification"]["overwrite"].lower() == "true"
    prepare_classify(ifn, ofn=ofn, overwrite=ovwr)
    
    logger.info("Preparing random file for classification...")
    ifn = file4("random")
    ofn = config["Classification"]["random_filename"]
    ovwr = config["Classification"]["overwrite"].lower() == "true"
    prepare_classify(ifn, ofn=ofn, overwrite=ovwr)
    
 
if (config["Data preparation"]["fake_positions"].lower() == "true" and args.steps==None) or (args.steps and "Data_preparation/fake_positions" in args.steps[0]):
    fake_positions(cconfig=config)


