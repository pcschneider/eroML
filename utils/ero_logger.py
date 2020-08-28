import logging.handlers
import logging

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
