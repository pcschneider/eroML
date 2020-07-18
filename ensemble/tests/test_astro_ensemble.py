import unittest
from eroML.ensemble import Ensemble, Astro_Object
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

class TestEnsemble(unittest.TestCase):

    def setUp(self):
        e = Ensemble()

        a = Astro_Object({"srcID":"a", "coord":SkyCoord(ra=1.0*u.degree, dec=-1.0*u.degree), "pm":(10,-10)}, pm_name="pm")
        b = Astro_Object({"srcID":"b", "coord":SkyCoord(ra=2.0*u.degree, dec=-2.0*u.degree), "pm":(20,-20)}, pm_name="pm")
        c = Astro_Object({"srcID":"c", "coord":SkyCoord(ra=3.0*u.degree, dec=-3.0*u.degree), "pm":(30,-30)}, pm_name="pm")
        d = Astro_Object({"srcID":"d", "coord":SkyCoord(ra=4.0*u.degree, dec=-4.0*u.degree), "pm":(40,-40)}, pm_name="pm")
        e = Ensemble()
        #print("Setup: ",e.known_cols)
        for x in [a,b,c,d]:
            e.add_one_object(x)
        
        self.e = e
        
    def test_getitem(self):
        self.assertEqual(self.e["c"]["srcID"],"c")
    
    def test_add_object(self):
        with self.assertRaises(TypeError):
            self.e.add_one_object("1")

    def test_remove_non_existing_object(self):
        with self.assertRaises(IndexError):
            self.e.del_object("1")
            
    def test_add_existing_object_no_auto_resolve(self):
        a = Astro_Object({"srcID":"a", "coord":SkyCoord(ra=1.0*u.degree, dec=-1.0*u.degree), "pm":(10,-10)}, pm_name="pm")
        with self.assertRaises(KeyError):
            self.e.add_one_object(a, auto_resolve=False)

    def test_add_existing_object(self):
        a = Astro_Object({"srcID":"a", "coord":SkyCoord(ra=1.0*u.degree, dec=-1.0*u.degree), "pm":(10,-10)}, pm_name="pm")
        self.assertEqual(self.e.add_one_object(a), "a_2")


    def test_remove_existing_object(self):
        self.assertEqual(self.e["c"]["srcID"],"c")
        self.e.del_object("b")
        self.assertEqual(len(self.e), 3)
        self.assertEqual(self.e["c"]["srcID"],"c")

    def test_init(self):
        self.assertEqual(len(self.e),4,"Should be 5.")
    
    def test_columns(self):
        xx = self.e.to_array(colnames=("RA","Dec"), array_type="array")
        #print("xxx",xx, np.shape(xx), type(xx))
        #print(np.array(xx), xx.shape)
        self.assertTrue(np.allclose(np.shape(xx), (2,4)))
        
        
        
if __name__ == '__main__':
    unittest.main()
