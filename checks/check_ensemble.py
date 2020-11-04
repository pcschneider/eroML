from eroML.ensemble import Ensemble, Astro_Object
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from eroML.ensemble.tools import from_fits, to_fits

#a = Astro_Object({"srcID":"a", "coord":SkyCoord(ra=1.0*u.degree, dec=-1.0*u.degree), "pm":(10,-10)}, pm_name="pm")
#b = Astro_Object({"srcID":"b", "coord":SkyCoord(ra=2.0*u.degree, dec=-2.0*u.degree), "pm":(20,-20)}, pm_name="pm")
#c = Astro_Object({"srcID":"c", "coord":SkyCoord(ra=3.0*u.degree, dec=-3.0*u.degree), "pm":(30,-30)}, pm_name="pm")
#e = Ensemble()

#for x in [a,b,c]:
    #e.add_object(x)
    
#print(len(e), e._N)


#xx = e.array(colnames=("ra","dec"))
#print(np.shape(xx))

#print(e.skyCoords())



x = from_fits("../eFEDS/SrcCat_V2T.fits", mapper={"detUID":"srcID"}, maxN=100)
to_fits(x, "test.fits", overwrite=True, maxN = 10)
y = from_fits("test.fits")
