import matplotlib.pyplot as plt
import numpy as np
from eroML.tile import file4
from eroML.ensemble import from_fits
from eroML.positions import analytic_probability, calc_sigma_from_RADEC_ERR
#import matplotlib 
plt.rcParams.update({'font.size': 16})
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.96, top=0.96, wspace=0, hspace=0)

mfn = file4("major", cconfig="eFEDS_EDR3.ini")
e = from_fits(mfn)
print("Using merged file ",mfn," with ",len(e), " entries")
tfn0 = "svm_training_IDs.txt"
tfn1 = "eFEDS_final_training_set.txt"
tfn1 = "eFEDS_good_pos.txt"

t0 = np.genfromtxt(tfn0, dtype=str)
t1tmp = np.genfromtxt(tfn1)
t1 = np.array([str("ML%05i" % d) for d in t1tmp])

print("and training IDs ",tfn0, tfn1, " with ", len(t0), len(t1), " entries.")
md0 = e.to_array("match_dist", array_type="array", srcIDs=t0)[0]
md1 = e.to_array("match_dist", array_type="array", srcIDs=t1)[0]
#print(len(md0), len(md1))

pp = e.to_array(["match_dist", "eligible_sky_density", "RADEC_ERR"], array_type="dict")
pp0 = e.to_array(["match_dist", "eligible_sky_density", "RADEC_ERR"], array_type="dict", srcIDs=t0)
pp1 = e.to_array(["match_dist", "eligible_sky_density", "RADEC_ERR"], array_type="dict", srcIDs=t1)
#print(pp0)

ap = analytic_probability(pp["match_dist"], sigma=calc_sigma_from_RADEC_ERR(pp["RADEC_ERR"]), sky_density=pp["eligible_sky_density"], ps=0.07)
ap0 = analytic_probability(pp0["match_dist"], sigma=calc_sigma_from_RADEC_ERR(pp0["RADEC_ERR"]), sky_density=pp0["eligible_sky_density"], ps=0.07)
ap1 = analytic_probability(pp1["match_dist"], sigma=calc_sigma_from_RADEC_ERR(pp1["RADEC_ERR"]), sky_density=pp1["eligible_sky_density"], ps=0.07)

plt.hist(ap, bins=30, range=(0,1), alpha=0.3, label="All")
plt.hist(ap0, bins=30, range=(0,1), alpha=0.5, label="SVM")
plt.hist(ap1, bins=30, range=(0,1), alpha=0.5, label="Bayes")
plt.legend()
plt.xlim(0.6,1.)
plt.ylim(0, 230)
plt.xlabel("Geometric match probability")
plt.ylabel("N")
plt.show()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.96, wspace=0, hspace=0)

plt.hist(pp0["match_dist"]/calc_sigma_from_RADEC_ERR(pp0["RADEC_ERR"]), label="md0", range=(0,2), bins=30, density=True, alpha=0.5)
plt.hist(pp1["match_dist"]/calc_sigma_from_RADEC_ERR(pp1["RADEC_ERR"]), label="md1", range=(0,2), bins=30, density=True, alpha=0.5)
plt.xlabel(r"r($\sigma$)")
plt.ylabel("N")
plt.legend()
plt.show()
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.98, top=0.96, wspace=0, hspace=0)

plt.scatter(calc_sigma_from_RADEC_ERR(pp0["RADEC_ERR"]), pp0["match_dist"], label="SVM", alpha=0.5)
plt.scatter(calc_sigma_from_RADEC_ERR(pp1["RADEC_ERR"]), pp1["match_dist"], label="Bayes", alpha=0.5)
plt.legend()
plt.xlabel("$\sigma$ (arcsec)")
plt.ylabel("Match distance (arcsec)")
plt.ylim(0,3.5)
plt.xlim(0.9,6.2)
plt.show()
