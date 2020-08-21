import logging
from .enrich import enrich_eROSITA

logger = logging.getLogger('eroML')


def ero_tile_loop(idx, prefix=None, postfix=None):
    """
    Loop through ero tiles and enrich them
    """
    logger.info("Enrichting %i eROISTA source tiles." % len(idx))
    
    for j, i in enumerate(idx):
        fn = prefix+str(i)+postfix+'.fits'
        logger.debug("Enriching eROISTA tile: %s (file %i/%i." % (fn, j+1, len(idx))) 
        try:
            enrich_eROSITA(fn, mapper={"DETUID":"srcID", "RA_CORR":"RA", "DEC_CORR":"Dec"})
        except:
            enrich_eROSITA(fn)
            
