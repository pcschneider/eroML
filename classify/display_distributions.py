from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
from eroML.tile import file4

fn = ""
fn = file4("major", cconfig="eFEDS_EDR3.ini")
print("Reading: ",fn)
ff = pyfits.open(fn)
fd = ff[1].data
print("Which contains ",fd.columns)

props = ["match_dist", "parallax", "Fg", "Fx", "bp_rp"]
for p in props:
    q0, q1 = np.nanquantile(fd[p], 0.1), np.nanquantile(fd[p], 0.9)
    print(p, q0, q1)
    plt.hist(fd[p],range=(0.8*q0, 1.2*q1), bins=20)
    plt.title(p)

    plt.show()


