from .ensemble import *

def test():
    modules = ["ensemble"]
    x = []
    for m in modules:
        mod = __import__(m, globals(), locals(), ['test'])
        tstfunc = getattr(mod, 'test')
        print(mod)
        x.append(m)
    print(x)    
