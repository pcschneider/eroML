from .astro_object import Astro_Object
from .astro_ensemble import Ensemble

def test():
    import unittest
    tl = unittest.TestLoader()
    tests = tl.discover(".")
    testRunner = unittest.runner.TextTestRunner()
    r = testRunner.run(tests)
    print(r)
