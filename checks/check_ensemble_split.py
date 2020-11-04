from eroML.ensemble import Ensemble, Astro_Object, fake_ensemble
from eroML.utils import sky_density
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

e = fake_ensemble(300, center=(0,0), width=(0.5, 0.2), random_pos=False)
#print(e)
g = sky_density(e, around=3, filter_prop=None, out_col="sky_density")

#print()
#print("aaaaaaaa",g.known_cols, len(g), g)
skd = g.to_array("sky_density", array_type="array")
#print(skd)

c = g.skyCoords()
plt.scatter(c.ra.degree, c.dec.degree, c=skd)
plt.show()

#a,b,c,d = e.split(4)
#print(a)

#print(a.srcIDs())
