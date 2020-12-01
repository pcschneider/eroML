from .ensemble import *
from .utils import *
from .tile import *
from .classify import *
from. positions import *

def test():
    modules = ["ensemble", "utils","tile"]
    x = []
    for m in modules:
        mod = __import__(m, globals(), locals(), ['test'])
        print("running test() from ",mod)
        
        tstfunc = getattr(mod, 'test')
        tstfunc()
        x.append(m)
        print()
    print(x)    


