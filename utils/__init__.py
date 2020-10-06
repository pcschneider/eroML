from eroML.utils.datasets import major_loop, random_loop, training_loop, major_set, training_set, random_set,file_loop_1to1,file_loop_2to1, shrink
from eroML.utils.enrich import enrich_Gaia, enrich_eROSITA, enrich_merged, sky_density
from eroML.utils.estimators import NN_distribution, Dist_model
from eroML.utils.gaia_tools import get_gaia, download_Gaia_tiles, Gaia_tile_loop
from eroML.utils.iso_tools import *
from eroML.utils.ero_tools import X_tile_loop
from .gaia_tools import vo2fits
from .ero_logger import setup_logger

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
            #for ttt in tt._tests:
                #print("Discovered test functions:", ttt._testMethodName)
                #print("   ",ttt.id())
            print()
    testRunner = unittest.runner.TextTestRunner()
    r = testRunner.run(test_suits)
    return r

