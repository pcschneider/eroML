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

props = ["RADEC_ERR", "match_dist", "parallax", "Fg", "Fx", "bp_rp"]
for p in props:
    
    #fig = plt.figure()
    #ax0, ax
    fig, (ax0, ax1) = plt.subplots(2)
    q0, q1 = np.nanquantile(fd[p], 0.02), np.nanquantile(fd[p], 0.95)
    print(p, q0, q1)
    print("    mean, median: ",np.nanmean(fd[p]), np.nanmedian(fd[p]))
    ax0.hist(fd[p],range=(0.5*q0, 5*q1), bins=20, log=True)
    ax1.hist(np.log10(fd[p]),range=(np.log10(0.5*q0), np.log10(5*q1)), bins=20, log=True)
    ax0.set_title(p)

    plt.show()


