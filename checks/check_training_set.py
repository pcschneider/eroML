from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt

fn = "Tile_training.fits"
fn = "tile/test_merged.fits"

fn = "../../ero_data/training_eFEDS.fits"
#fn = "../ero_data/merged_random_test120.fits"

ff = pyfits.open(fn)[1].data

gi0 = np.where(ff["category"] >0)[0]
plt.scatter(ff["bp_rp"][gi0], np.log10(ff["logFxFg"][gi0]), color='r', marker='.', label="others")

gi1 = np.where(ff["category"] == 0)[0]

plt.scatter(ff["bp_rp"][gi1], np.log10(ff["logFxFg"][gi1]), color='b', marker='.', label="Stars")

#print(ff["Gaia_Quality"], ff["Gaia_Quality"].dtype, np.shape( ff["Gaia_Quality"]))
#print(ff["category"], ff["category"].dtype, np.shape( ff["category"]))

#gi3 = np.where( (ff["category"] == 1) & (ff["Gaia_Quality"] == 1))[0]
#plt.scatter(ff["bp_rp"][gi3], np.log10(ff["FxFg"][gi3]), label="stars")

#gi1 = np.where( (ff["category"] == 0) & (ff["Gaia_Quality"] == 1))[0]
#plt.scatter(ff["bp_rp"][gi1], np.log10(ff["FxFg"][gi1]), label="others")


#gi2 = np.where( (ff["category"] == 1) & (ff["Gaia_Quality"] == 0) & (ff["iso_compatible"]=="True") )[0]
#print(len(gi2))
#plt.scatter(ff["bp_rp"][gi2], np.log10(ff["logFxFg"][gi2]), label="stars/iso", s=15, facecolor="none", edgecolor='b')


plt.xlabel("bp_rp")
plt.ylabel("log Fx/Fg")

#print("#stars: ",len(gi3), " #others ",len(gi1), " #all ",len(ff["category"]))

plt.legend()
plt.show()

gi = np.where(ff["category"] > 0)
plt.hist(ff["match_dist"][gi], label="others", bins=20, range=(0,100))

gi = np.where(ff["category"] == 0)
plt.hist(ff["match_dist"][gi],label="stars", bins=20, range=(0,100))

plt.legend()
plt.show()

#plt.figure()
#gi1 = np.where(ff["category"] == 0)[0]
#plt.hist(ff["parallax"][gi1]/ ff["parallax_error"][gi1], range=(-5,10), bins=31, alpha=0.5)
#print("a",len(gi1))

#gi1 = np.where((ff["category"] == 1) & (ff["iso_compatible"] == "True")  )[0]
#print("b",len(gi1))
#plt.hist(ff["parallax"][gi1]/ ff["parallax_error"][gi1], range=(-5,10), bins=31, alpha=0.5)

#gi1 = np.where((ff["category"] == 1) & (ff["Gaia_Quality"] == 1)  )[0]
#print("c",len(gi1))
#plt.hist(ff["parallax"][gi1]/ ff["parallax_error"][gi1], range=(-5,10), bins=31, alpha=0.5)

#gi1 = np.where(ff["category"] == 1)[0]
#print("d",len(gi1))
#plt.hist(ff["parallax"][gi1]/ ff["parallax_error"][gi1], range=(-5,10), bins=31, alpha=0.5)



#plt.figure()
#gi1 = np.where((ff["category"] == 1) & (ff["Gaia_Quality"] == 0) & (ff["iso_compatible"] == "False")  )[0]
#print("x",len(gi1))
#plt.hist(ff["phot_g_mean_mag"][gi1], bins=31, alpha=0.5)
#ll = np.where(ff["phot_g_mean_mag"][gi1] < 17)[0]
#print("y",len(ll))
#print(ff["srcID_NN"][gi1[ll]])

#plt.figure()
#gi1 = np.where((ff["category"] == 1) & (ff["Gaia_Quality"] == 0) & (ff["iso_compatible"] == "False")  )[0]
#print("z",len(gi1))
#plt.hist(ff["Fx"][gi1], bins=31, alpha=0.5)
#print("y",len(np.where(ff["Fx"][gi1] < 17)[0]))
plt.show()

