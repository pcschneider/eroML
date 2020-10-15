import logging
from configparser import ConfigParser, ExtendedInterpolation
import eroML

def read_config(fn):
    cc = ConfigParser(interpolation=ExtendedInterpolation())
    rr = cc.read(fn)
    if len(rr) == 0:
        return ("error","Config-file \'%s\' does not contain info (it may even not exist)." % fn)
    #logger.debug("Reading config-file: %s" %fn)
    for sec in cc.sections():
        for k in cc[sec].keys():
            config[sec][k] = cc.get(sec, k)
    return ("debug","Reading config-file: %s" %fn)


dd = eroML.__path__[0]
config = ConfigParser(interpolation=ExtendedInterpolation())
config.read(dd+"/"+"default.ini")

if __name__=="__main__":
    pass
