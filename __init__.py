from .ensemble import *
from .utils import *

def test():
    modules = ["ensemble", "utils"]
    x = []
    for m in modules:
        mod = __import__(m, globals(), locals(), ['test'])
        print("running test() from ",mod)
        
        tstfunc = getattr(mod, 'test')
        tstfunc()
        x.append(m)
        print()
    print(x)    
