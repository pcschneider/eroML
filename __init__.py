from .ensemble import *

def test():
    modules = ["ensemble"]

    for m in modules:
        mod = __import__(m, globals(), locals(), ['test'])
        tstfunc = getattr(mod, 'test')
        tstfunc()
