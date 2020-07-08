from .astro_object import Astro_Object
from .astro_object import from_Simbad
from .astro_ensemble import Ensemble

def test():
    import unittest
    tl = unittest.TestLoader()
    test_suits = tl.discover(".")
    print("tests:", test_suits)
    for i, t in enumerate(test_suits):
        print(i,t)
        #print(len(t))
        print()
    testRunner = unittest.runner.TextTestRunner()
    r = testRunner.run(test_suits)
    print(r)
    return tests
