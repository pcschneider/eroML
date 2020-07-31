from eroML.utils.datasets import major_set, training_set, random_set
from eroML.utils.enrich import enrich_Gaia, enrich_eROSITA, enrich_merged
from eroML.utils.estimators import NN_distribution
from eroML.utils.gaia_tools import prepare_gaia
from eroML.utils.iso_tools import *

def test():
    import doctest
    
    import eroML.utils.estimators
    r = doctest.testmod(eroML.utils.estimators)
    print("Result of doctest for `NN_distribution`:", r)
    
    #import eroML.ensemble.astro_ensemble
    #r = doctest.testmod(astro_ensemble)
    #print("Result of doctest for `astro_ensemble`:", r)

    #import eroML.ensemble.tools
    #r = doctest.testmod(tools)
    #print("Result of doctest for `tools`:", r)
    #print()
    
    import unittest
    tl = unittest.TestLoader()
    test_suits = tl.discover(".")
    for i, t in enumerate(test_suits):
        for j, tt in enumerate(t):
            print(t)
            for ttt in tt._tests:
                print("Discovered test functions:", ttt._testMethodName)
                print("   ",ttt.id())
            print()
    testRunner = unittest.runner.TextTestRunner()
    r = testRunner.run(test_suits)
    return r

