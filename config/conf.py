import logging
from configparser import ConfigParser, ExtendedInterpolation


def read_config(fn):
    cc = ConfigParser(interpolation=ExtendedInterpolation())
    cc.read(fn)
    logger.debug("Reading config-file: %s" %fn)
    for sec in cc.sections():
        for k in cc[sec].keys():
            config[sec][k] = cc.get(sec, k)

logger = logging.getLogger('eroML')

config = ConfigParser(interpolation=ExtendedInterpolation())
config.read("default.ini")

if __name__=="__main__":
    pass
