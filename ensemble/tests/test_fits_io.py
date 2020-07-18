from eroML.ensemble.tools import from_fits, to_fits
import unittest

class Test_Fits_IO(unittest.TestCase):
    def test_fits_read(self):
        x = from_fits("../eFEDS/SrcCat_V2T.fits", mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=100)
        self.assertEqual(x["ML00001"]["srcID"].strip(), "ML00001")
        
    def test_fits_io(self):
        x = from_fits("../eFEDS/SrcCat_V2T.fits", mapper={"detUID":"srcID", "DEC":"Dec"}, maxN=100)
        to_fits(x, "test.fits", overwrite=True, maxN = 10)
        y = from_fits("test.fits")
        #print("yyy")
        xx = y["ML00001"]
        #print("xx",xx)
        self.assertEqual(y["ML00001"]["srcID"].strip(), "ML00001")
        
if __name__ == '__main__':
    #pass
    unittest.main()
        
