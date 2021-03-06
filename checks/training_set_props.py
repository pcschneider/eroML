import matplotlib.pyplot as plt
import numpy as np
from eroML.tile import file4
from eroML.ensemble import from_fits, to_fits
from eroML.positions import analytic_probability, calc_sigma_from_RADEC_ERR
from eroML.utils import activity_filter

plt.rcParams.update({'font.size': 16})
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.96, top=0.96, wspace=0, hspace=0)

activity_poly = [-3.22, 3.6/5.5]

mfn = file4("major", cconfig="eFEDS_EDR3_HamStar.ini")
major = from_fits(mfn)
print("Using merged file ",mfn," with ",len(major), " entries")

tfn0 = "svm_training_IDs.txt"
tfn0 = "positions/training_IDs4.txt"
tfn1 = "Bayes_eFEDS_good_pos.txt"
tfn1 = "Bayes_eFEDS_final_training_set.txt"
t0 = np.genfromtxt(tfn0, dtype=str)
t1tmp = np.genfromtxt(tfn1)
t1 = np.array([str("ML%05i" % d) for d in t1tmp])
print("and training IDs ",tfn0, tfn1, " with ", len(t0), len(t1), " entries.")

color = major.to_array("bp_rp", array_type="array")
FxFg = major.to_array("FxFg", array_type="array")
below = (1-activity_filter(color, FxFg)).astype(bool)
well_above = activity_filter(color, FxFg, log_margin=0.5).astype(bool)
cl = np.ones(len(major)).astype(int) *2
cl[below] = 0
cl[well_above] = 1
gi2 = np.where(cl < 1)[0]

#plt.scatter(color, FxFg, label="all ("+str(len(color))+")")
#plt.scatter(color[gi2], FxFg[gi2], label="bona fide stars ("+str(len(gi2))+")")
#plt.xlabel("BP - RP")
#plt.ylabel("FxFg")
#plt.yscale("log")
#plt.ylim(1e-7,1)
#plt.legend()
#plt.show()

dist = 1000/major.to_array("parallax", array_type="array")
Fx = major.to_array("Fx", array_type="array")
Lx = 4*np.pi*(3.18e-2*dist)**2 * Fx * 1e40
gi = np.where(((Lx > 2e31) & (cl==0)) | (dist>1500))[0]
print("Below saturation limit, but too high Lx", len(gi))
cl[gi] = 2
print("Leaving ",len(np.where(cl==0)[0]), " bona fide stars.")


major.add_col("category", cl)
to_fits(major,"x.fits", overwrite=True)



c = major.to_array("category", array_type="array", srcIDs=t0)
#print(np.shape(c))
gi = np.where(c[0]==0)[0]
#print(c)
srcIDs0 = major.to_array("srcID", array_type="array", srcIDs=t0[gi])
srcIDs1 = major.to_array("srcID", array_type="array", srcIDs=t1)
color0 = major.to_array("bp_rp", array_type="array", srcIDs=t0[gi])[0]
FxFg0 = major.to_array("FxFg", array_type="array", srcIDs=t0[gi])[0]
color1 = major.to_array("bp_rp", array_type="array", srcIDs=t1)[0]
FxFg1 = major.to_array("FxFg", array_type="array", srcIDs=t1)[0]

print("stars in t0: ",len(srcIDs0[0]))
print("stars in t1: ",len(srcIDs1[0]))


srcIDs = np.unique(np.concatenate((srcIDs0[0], srcIDs1[0])))
color = major.to_array("bp_rp", array_type="array", srcIDs=srcIDs)[0]
FxFg = major.to_array("FxFg", array_type="array", srcIDs=srcIDs)[0]

print("Number of sources within t0 and t1: ",len(srcIDs))
e = major
pp = e.to_array(["match_dist", "eligible_sky_density", "RADEC_ERR"], array_type="dict", srcIDs=srcIDs)
pp0 = e.to_array(["match_dist", "eligible_sky_density", "RADEC_ERR"], array_type="dict", srcIDs=t0[gi])
pp1 = e.to_array(["match_dist", "eligible_sky_density", "RADEC_ERR"], array_type="dict", srcIDs=t1)
#print(pp0)

ap = analytic_probability(pp["match_dist"], sigma=calc_sigma_from_RADEC_ERR(pp["RADEC_ERR"]), sky_density=pp["eligible_sky_density"], ps=0.07)
ap0 = analytic_probability(pp0["match_dist"], sigma=calc_sigma_from_RADEC_ERR(pp0["RADEC_ERR"]), sky_density=pp0["eligible_sky_density"], ps=0.07)
ap1 = analytic_probability(pp1["match_dist"], sigma=calc_sigma_from_RADEC_ERR(pp1["RADEC_ERR"]), sky_density=pp1["eligible_sky_density"], ps=0.07)

#print(ap0[si[0]], ap0[si[-1]])
#print(max(ap1), np.quantile(ap1, 0.01))
x = np.linspace(0,4, 100)
ys = activity_poly[1]
y = x*ys + activity_poly[0]
plt.plot(x, 10**y, color='k')

#srcIDs = major.to_array("srcID", array_type="array")
print("Merged training set size:",len(srcIDs))
svm = np.zeros(len(srcIDs))
svm[np.in1d(srcIDs, srcIDs0)] = 1
#major.add_col("svm_training", svm)
seb = np.zeros(len(srcIDs))
seb[np.in1d(srcIDs, srcIDs1)] = 1
#major.add_col("bayes_training", seb)
#to_fits(major,"x.fits", overwrite=True)

si = np.where((svm == 1) & (seb == 0))[0]
bi = np.where((svm == 0) & (seb == 1))[0]
plt.scatter(color, FxFg, label="All", c=ap, vmin=0.6, vmax=1.)
plt.scatter(color[si], FxFg[si], label="SVM", fc="None", ec="k", s=28)
plt.scatter(color[bi], FxFg[bi], label="Bayesian", fc="None", ec="r", s=24)
cb = plt.colorbar()
cb.set_label("Analytic probability")

from astropy.io import fits as pyfits
fn = "train_HamStar_preprocessed.fits"
ff = pyfits.open(fn)
fd = ff[1].data
gi = np.where(fd["category"]==0)[0]
print(" #Stellar training objects in ",fn,":",len(gi))
a, b = fd["bp_rp"][gi], 10**fd["FxFg"][gi]
plt.scatter(a, b, s=15, marker='x', color='g', label=fn+" stars")
gi = np.where(fd["category"]>0)[0]
a, b = fd["bp_rp"][gi], 10**fd["FxFg"][gi]
plt.scatter(a, b, s=15, marker='x', color='r', label=fn +" random")


plt.yscale("log")
plt.ylim(2e-7, 1e-1)
plt.legend()
plt.xlabel("BP-RP (mag)")
plt.ylabel("Fx/Fbol")


plt.show()
exit()


#md0 = e.to_array("match_dist", array_type="array", srcIDs=t0)[0]
#md1 = e.to_array("match_dist", array_type="array", srcIDs=t1)[0]
##print(len(md0), len(md1))

#plt.hist(ap, bins=30, range=(0,1), alpha=0.3, label="All")
#plt.hist(ap0, bins=30, range=(0,1), alpha=0.5, label="SVM")
#plt.hist(ap1, bins=30, range=(0,1), alpha=0.5, label="Bayes")
#plt.legend()
#plt.xlim(0.6,1.)
#plt.ylim(0, 230)
#plt.xlabel("Geometric match probability")
#plt.ylabel("N")
#plt.show()


#plt.hist(pp0["match_dist"]/calc_sigma_from_RADEC_ERR(pp0["RADEC_ERR"]), label="md0", range=(0,2), bins=30, density=True, alpha=0.5)
#plt.hist(pp1["match_dist"]/calc_sigma_from_RADEC_ERR(pp1["RADEC_ERR"]), label="md1", range=(0,2), bins=30, density=True, alpha=0.5)
#plt.xlabel(r"r($\sigma$)")
#plt.ylabel("N")
#plt.legend()
#plt.show()
