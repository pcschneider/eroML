import logging
from .enrich import enrich_eROSITA, enrich_ROSAT
#from astropy.io import fits as pyfits
from ..ensemble import from_fits, to_fits

logger = logging.getLogger('eroML')


def X_tile_loop(idx, prefix=None, postfix=None):
    """
    Loop through ero tiles and enrich them
    """
    logger.info("Enrichting %i X-ray source tiles." % len(idx))
    
    for j, i in enumerate(idx):
        fn = prefix+str(i)+postfix+'.fits'
        logger.info("Enriching X-ray tile: %s (file %i/%i." % (fn, j+1, len(idx))) 
        
        
        
        e = None
        
        try: # Genuie eROSITA source file
            logger.debug("Trying eROSITA source file...")
            e = from_fits(fn, mapper={"DETUID":"srcID", "RA_CORR":"RA", "DEC_CORR":"Dec"})
            e.instrument = "eROSITA"
            logger.debug("Successfully read regular eROSITA file: \'%s'\." % fn)
        except:
            logger.debug("Not an original eROSITA source file (%s)" % fn)
        
        
        try:    
            logger.debug("Trying ROSAT source file...")
            e = from_fits(fn, mapper={"IAU_NAME":"srcID", "RA_DEG":"RA", "DEC_DEG":"Dec"})
            e.instrument = "ROSAT"
            logger.debug("Successfully read regular ROSAT file: \'%s\'." % fn)
        except:
            logger.debug("Not an original ROSAT source file (%s)" % fn)
            
        try: # Already enriched source files
            logger.debug("Trying regular source file...")
            e = from_fits(fn)
            logger.debug("Successfully read regular source file: \'%s\'." %fn)
        except:
            logger.debug("Not a regular source file (%s)" % fn)
            raise Exception("Cannot read source file (%s)" % fn)
        
        if "RATE_1rxs" in e.known_cols:
            logger.debug("Enriching a ROSAT source file...")
            enrich_ROSAT(e)
        else:            
            logger.debug("Enriching a ROSAT source file...")
            enrich_eROSITA(e)
            
        to_fits(e, ofn=fn, overwrite=True)
