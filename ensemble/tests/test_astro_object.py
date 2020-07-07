import unittest
from eroML.ensemble import Astro_object
from astropy.coordinates import SkyCoord
import astropy.units as u

class TestAstro_object(unittest.TestCase):

    def test_init(self):
        a = Astro_object({"srcID":"a", "coord":SkyCoord(ra=150.0*u.degree, dec=20.0*u.degree), "pm":(10,-50)}, pm_name="pm")        
        self.assertEqual(a.srcID, "a", "Should be \"a\"")

    def test_sum_tuple(self):
        pass
        #self.assertEqual(sum((1, 2, 2)), 6, "Should be 6")

if __name__ == '__main__':
    unittest.main()
