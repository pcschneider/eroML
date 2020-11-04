import numpy as np
from astropy.time import Time,TimeDelta

t = Time( 5.154387500000000E+04, format='mjd')
t = Time(54101, format="mjd")

t.format = 'isot'
print(t)
dd = TimeDelta(6.3e8  - 5.154388e4*365*86400, format='sec')
print(dd)
t+=dd
print(t)
