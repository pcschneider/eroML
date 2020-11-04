from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt


c = SkyCoord(ra=219.87383306, dec=-60.83222194, unit=(u.degree, u.degree))
pmra, pmdec = -3608, 686
dtime=21


angle = np.arctan2(pmra, pmdec)# + np.pi*1
print(angle, angle*180/np.pi)

d_ra = pmra * dtime / 1000
d_dec = pmdec * dtime /1000
dd = (d_ra**2 + d_dec**2)**0.5

d = c.directional_offset_by(angle/np.pi*180*u.degree, dd*u.arcsec )
print("pm: ",dd)
print(c.ra.degree, c.dec.degree)
print(d.ra.degree, d.dec.degree)

# 2000
#219.87383306 -60.83222194
plt.scatter([219.87383306], [-60.83222194], label="2000")

## 2021
             #219.83065418598989 -60.82821334851694
             #219.83064879, -60.83220918  
plt.scatter([219.83065419], [-60.82821335], label="Simbad 2021") # <- Simbad
plt.scatter(d.ra.degree, d.dec.degree, label="this")
plt.legend()
plt.show()
