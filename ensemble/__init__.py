from .astro_object import Astro_object
from .ensemble import Ensemble

def test(level=1, verbosity=1):
    from numpy.testing import Tester
    return Tester().test(level, verbosity)
