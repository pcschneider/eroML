"""
This module handles Ensembles, i.e., collections of astronomical objects, which are mainly characterized by a srcID and coordinates.    
"""

from .astro_object import Astro_Object
from .astro_object import from_Simbad
from .astro_ensemble import Ensemble
from .tools import from_fits, to_fits
from astropy.coordinates import SkyCoord
from astropy import units as u

def test():
    import doctest
    
    import eroML.ensemble.astro_object
    r = doctest.testmod(astro_object,   extraglobs={'a': Astro_Object({"srcID":"a", "coord":SkyCoord(ra=150.0*u.degree, dec=20.0*u.degree), "pm":(10,-50)}, pm_name="pm")})
    print("Result of doctest for `astro_object`:", r)
    
    import eroML.ensemble.astro_ensemble
    r = doctest.testmod(astro_ensemble)
    print("Result of doctest for `astro_ensemble`:", r)

    import eroML.ensemble.tools
    r = doctest.testmod(tools)
    print("Result of doctest for `tools`:", r)
    print()
    
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
    return tests

