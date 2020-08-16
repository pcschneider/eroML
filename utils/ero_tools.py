import logging
from .enrich import enrich_eROSITA

logger = logging.getLogger('eroML')


def ero_tile_loop(idx, prefix=None, postfix=None):
    """
    Loop through ero tiles and enrich them
    """
    logger.info("Enrichting %i eROISTA source tiles." % len(idx))
    
    for i in idx:
        fn = prefix+str(i)+postfix+'.fits'
        logger.debug("Enriching eROISTA tile: %s." % fn) 
        enrich_eROSITA(fn, mapper={"DETUID":"srcID", "DEC":"Dec"})
