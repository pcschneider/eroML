import unittest
from eroML.ensemble import Astro_Object, from_Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

class TestAstro_Object(unittest.TestCase):

    def test_init(self):
        a = Astro_Object({"srcID":"a", "coord":SkyCoord(ra=150.0*u.degree, dec=20.0*u.degree), "pm":(10,-50)}, pm_name="pm")        
        self.assertEqual(a.srcID, "a", "Should be \"a\"")

    def test_pm(self):
        a = Astro_Object({"srcID":"AU Mic", "coord":SkyCoord(ra=311.2897182064525 *u.degree, dec=-31.3409004752072*u.degree), "pm_Gaia":(281.424,-359.895)}, pm_name="pm_Gaia")  
        tmp = a.coord_tuple(epoch=2021)
        self.assertTrue(np.allclose(tmp, (311.29164037, -31.34299987), rtol=1e-9), "Should be True, because the difference should be small.")

    def test_Simbad(self):
        a = from_Simbad("AU Mic")
        
if __name__ == '__main__':
    unittest.main()
