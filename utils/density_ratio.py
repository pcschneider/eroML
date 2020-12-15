import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

dd = np.genfromtxt("offs.dat", unpack=True)
dd[1]*=100
X = np.transpose(dd[0:3])

gi = np.where(dd[4]==0)[0]

real = KernelDensity(kernel='gaussian', bandwidth=1).fit(X[gi])

