from eroML.tile import loop, file4, merge_fits
from eroML.tile import merger, add_healpix_col, hpix2process, generate_healpix_files
from eroML.utils import download_Gaia_tiles, Gaia_tile_loop
from eroML.utils import enrich_Gaia, X_tile_loop
from eroML.utils import major_loop, random_loop, training_loop
from eroML.utils import file_loop_1to1, shrink
from eroML.utils import setup_logger
import logging
from configparser import ConfigParser
from eroML.config import config, read_config

logger = logging.getLogger('eroML')

def custom_config(fn=None):
    

    if fn is not None:
        if type(fn) == type("XXX"):
            read_config(fn)
        elif type(fn) == ConfigParser:
            return fn
    return config    

def calculate_healpix(cconfig=None):
    cconfig = custom_config(cconfig)
    logger.info("Calculating healpix indices for eROSITA sources with ")
    logger.info("    ero_filename: %s" % cconfig["Sources"]["X_filename"])
    logger.info("    ofn: %s" % cconfig["Sources"]["X_filename_hp"])
    logger.info("    nside: %i" % cconfig["Healpix"].getint("nside"))
    add_healpix_col(cconfig["Sources"]["X_filename"], ofn=cconfig["Sources"]["X_filename_hp"], nside=cconfig["Healpix"].getint("nside"), overwrite=True)
    

def generate_ero_tiles(cconfig=None):
    cconfig = custom_config(cconfig)
    logger.info("Preparing Xray data")
    prex = cconfig["X data preparation"]["directory"]+"/"+cconfig["X data preparation"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    posx = ""
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)    
    generate_healpix_files(cconfig["Sources"]["X_filename_hp"], prefix=prex, postfix=posx, index0=index0, index1=index1, pix_file=healpix_file)

def perform_ero_data_preparation(cconfig=None):
    cconfig = custom_config(cconfig)
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    logger.debug("Enriching %i ero-tiles." % len(idx))
    prex = cconfig["X data preparation"]["directory"]+"/"+cconfig["X data preparation"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    posx = ""
    X_tile_loop(idx, prefix=prex, postfix=posx)

def perform_Gaia_download(cconfig=None):
    cconfig = custom_config(cconfig)
        
    logger.info("Downloading Gaia data")
    
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)
    
    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    logger.debug("Downloading %i Gaia-tiles." % len(idx))
    
    download_Gaia_tiles(outdir=cconfig["Gaia Download"]["directory"],\
               prefix=cconfig["Gaia Download"]["prefix"], idx=idx,\
               nside=int(cconfig["Healpix"]["nside"]),\
               overwrite=cconfig["Gaia Download"]["overwrite"],\
               edge=float(cconfig["Gaia Download"]["edge"]),\
               verbose=int(cconfig["Gaia Download"]["verbose"]),\
               Glim=float(cconfig["Gaia Download"]["Glim"]) )


def prepare_Gaia_data(cconfig=None):
    cconfig = custom_config(cconfig)
        
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    logger.debug("Enriching %i ero-tiles." % len(idx))
    prex = cconfig["Gaia Download"]["directory"]+"/"+cconfig["Gaia Download"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    posx = ""
    
    filt = int(config["Enrich Gaia"]["filter"])
    
    logger.info("Enriching Gaia data files (filt=%i)" % filt)
    Gaia_tile_loop(idx, prefix=prex, postfix=posx, filterNr=filt)
    


def generate_major_sets(cconfig=None):
    cconfig = custom_config(cconfig)
    
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    g_prex = cconfig["Gaia Download"]["directory"]+"/"+cconfig["Gaia Download"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    g_posx = ""
    
    e_prex = cconfig["X data preparation"]["directory"]+"/"+cconfig["X data preparation"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    e_posx = ""
    
    m_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["major_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    m_posx = ""
    
    logger.debug("Major sets for %i tiles." % len(idx))
    major_loop(idx, ero_prefix=e_prex, ero_postfix=e_posx, gaia_prefix=g_prex, gaia_postfix=g_posx, major_prefix=m_prex, major_postfix=m_posx)
    
   
def generate_random_sets(cconfig=None):
    cconfig=custom_config(cconfig)
    
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    g_prex = cconfig["Gaia Download"]["directory"]+"/"+cconfig["Gaia Download"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    g_posx = ""
    
    e_prex = cconfig["X data preparation"]["directory"]+"/"+cconfig["X data preparation"]["prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    e_posx = ""
    
    r_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["random_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    r_posx = ""
    
    mino = float(cconfig["Data sets"]["min_random_offset"])
    maxo = float(cconfig["Data sets"]["max_random_offset"])
    multi= int(cconfig["Data sets"]["random_multi"])
    logger.debug("Random data sets for %i tiles." % len(idx))
    random_loop(idx, ero_prefix=e_prex, ero_postfix=e_posx, gaia_prefix=g_prex, gaia_postfix=g_posx, random_prefix=r_prex, random_postfix=r_posx, min_offset=mino, max_offset=maxo, multi=multi)
       
   
def generate_training_sets(cconfig=None):
    cconfig = custom_config(cconfig)
    
    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    
    # random
    r_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["random_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    r_posx = ""
    
    #major
    m_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["major_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    m_posx = ""
    
    #training
    t_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["training_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    t_posx = ""
    
    
    ad = float(cconfig["Data sets"]["training_abs_dist"])
    rd = float(cconfig["Data sets"]["training_rel_dist"])
    
    logger.debug("Random data sets for %i tiles." % len(idx))
    training_loop(idx, major_prefix=m_prex, major_postfix=m_posx, random_prefix=r_prex, random_postfix=r_posx, training_prefix=t_prex, training_postfix=t_posx, abs_dist=ad, rel_dist=rd)
   

def shrinking(cconfig=None):
    cconfig = custom_config(cconfig)

    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = cconfig["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    
    r_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["random_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    r_posx = ""
    
    #major
    m_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["major_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    m_posx = ""
    
    #training
    t_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["training_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    t_posx = ""
 
    sh = cconfig["Merging"]["shrink_postfix"]
    cols = cconfig["Columns"]["keep"].split(",")
    #print("cols: ",cols)
    
    for tt in ["training","major","random"]:
        prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"][tt+"_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
        posx = ""
        file_loop_1to1(idx, prefix=prex, postfix=posx, method=shrink, ofn_prefix=prex, ofn_postfix=posx+sh,cols=cols)
    
   
#if config["Data sets"]["major"].lower()=="true":
    
    #healpix_file = config["Healpix"].get("pix_file", None)
    #index0 = config["Healpix"].getint("index0", 0)
    #index1 = config["Healpix"].getint("index1", None)

    #idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    #g_prex = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #g_posx = ""
    
    #e_prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #e_posx = ""
    
    #m_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #m_posx = ""
    
    #logger.debug("Major sets for %i tiles." % len(idx))
    #major_loop(idx, ero_prefix=e_prex, ero_postfix=e_posx, gaia_prefix=g_prex, gaia_postfix=g_posx, major_prefix=m_prex, major_postfix=m_posx)
    
        


#if config["Data sets"]["random"].lower()=="true":
    
    #healpix_file = config["Healpix"].get("pix_file", None)
    #index0 = config["Healpix"].getint("index0", 0)
    #index1 = config["Healpix"].getint("index1", None)

    #idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    #g_prex = config["Gaia Download"]["directory"]+"/"+config["Gaia Download"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #g_posx = ""
    
    #e_prex = config["eROSITA preparation"]["directory"]+"/"+config["eROSITA preparation"]["prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #e_posx = ""
    
    #r_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #r_posx = ""
    
    #mino = float(config["Data sets"]["min_random_offset"])
    #maxo = float(config["Data sets"]["max_random_offset"])
    #multi= int(config["Data sets"]["random_multi"])
    #logger.debug("Random data sets for %i tiles." % len(idx))
    #random_loop(idx, ero_prefix=e_prex, ero_postfix=e_posx, gaia_prefix=g_prex, gaia_postfix=g_posx, random_prefix=r_prex, random_postfix=r_posx, min_offset=mino, max_offset=maxo, multi=multi)
    




#if config["Data sets"]["training"].lower()=="true":
    
    #healpix_file = config["Healpix"].get("pix_file", None)
    #index0 = config["Healpix"].getint("index0", 0)
    #index1 = config["Healpix"].getint("index1", None)

    #idx = hpix2process(config["Sources"]["ero_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    
    ## random
    #r_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["random_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #r_posx = ""
    
    ##major
    #m_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["major_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #m_posx = ""
    
    ##training
    #t_prex = config["Data sets"]["directory"]+"/"+config["Data sets"]["training_prefix"]+"_nside"+config["Healpix"]["nside"]+"_"
    #t_posx = ""
    
    
    #ad = float(config["Data sets"]["training_abs_dist"])
    #rd = float(config["Data sets"]["training_rel_dist"])
    
    #logger.debug("Random data sets for %i tiles." % len(idx))
    #training_loop(idx, major_prefix=m_prex, major_postfix=m_posx, random_prefix=r_prex, random_postfix=r_posx, training_prefix=t_prex, training_postfix=t_posx, abs_dist=ad, rel_dist=rd)

def shrinking(cconfig=None):
    cconfig = custom_config(cconfig)

    healpix_file = cconfig["Healpix"].get("pix_file", None)
    index0 = cconfig["Healpix"].getint("index0", 0)
    index1 = config["Healpix"].getint("index1", None)

    idx = hpix2process(cconfig["Sources"]["X_filename_hp"], index0=index0, index1=index1, pix_file=healpix_file)
    
    r_prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"]["random_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
    r_posx = ""
 
    sh = cconfig["Merging"]["shrink_postfix"]
    cols = cconfig["Columns"]["keep"].split(",")
    logger.debug("Shrinking files and keeping only: %s " % str(cols))
    for tt in ["random","training","major"]:
        prex = cconfig["Data sets"]["directory"]+"/"+cconfig["Data sets"][tt+"_prefix"]+"_nside"+cconfig["Healpix"]["nside"]+"_"
        posx = ""
        file_loop_1to1(idx, prefix=prex, postfix=posx, method=shrink, ofn_prefix=prex, ofn_postfix=posx+sh,cols=cols)
    
    
#if config["Merging"]["major"].lower() == "true":
    #if config["Merging"]["shrink"].lower()=="true":
        #which = "major_tiles_small"
    #else:
        #which = "major_tiles"
    #ofn = file4("major")    
    #fnames = file4(which)
    #merge_fits(fnames, ofn=ofn)

#if config["Merging"]["random"].lower() == "true":
    #if config["Merging"]["shrink"].lower()=="true":
        #which = "random_tiles_small"
    #else:
        #which = "random_tiles"
    #ofn = file4("random")    
    #fnames = file4(which)
    #merge_fits(fnames, ofn=ofn)
    
#if config["Merging"]["training"].lower() == "true":
    #if config["Merging"]["shrink"].lower()=="true":
        #which = "training_tiles_small"
    #else:
        #which = "training_tiles"
    #ofn = file4("training")    
    #fnames = file4(which)
    #merge_fits(fnames, ofn=ofn)    

    
