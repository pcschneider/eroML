from .astro_object import Astro_Object
from .astro_object import from_Simbad
from .astro_ensemble import Ensemble

def test():
    import unittest
    tl = unittest.TestLoader()
    test_suits = tl.discover(".")
    print("tests:")
    for t in test_suits:
        print(t)
        for tst in t._tests: 
            print(tst, tst._tests)
        print()
    testRunner = unittest.runner.TextTestRunner()
    r = testRunner.run(test_suits)
    print(r)
    return tests

#def print_suite(suite):
    #if hasattr(suite, '__iter__'):
        #for x in suite:
            #print_suite(x)
    #else:
        #print(suite)

#print_suite(unittest.defaultTestLoader.discover('.'))
