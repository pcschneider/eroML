from .ensemble import *

def test():
    modules = ["ensemble"]
    x = []
    for m in modules:
        mod = __import__(m, globals(), locals(), ['test'])
        print(mod)
        tstfunc = getattr(mod, 'test')
        tstfunc()
        x.append(m)
    print(x)    
