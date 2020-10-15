from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import glob


dd = np.genfromtxt("dens_log.dat", usecols=(1,2))
print(np.shape(dd))
