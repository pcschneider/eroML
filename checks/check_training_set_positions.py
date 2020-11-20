from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt


fn0 = "../../ero_data/training_eFEDS.fits"
#fn = "../ero_data/merged_random_test120.fits"
ff0 = pyfits.open(fn0)[1].data

fn1 = "../../ero_data/random_eFEDS.fits"
#fn = "../ero_data/merged_random_test120.fits"
ff1 = pyfits.open(fn1)[1].data


gi = np.where(ff1["NN"] ==1)[0]
plt.hist(ff1["match_dist"][gi], label='random', bins=20, range=(0,100), density=True)



gi = np.where((ff0["category"] > 0) & (ff0["NN"] == 1))
plt.hist(ff0["match_dist"][gi], label="training random", bins=20, range=(0,100), density=True, alpha=0.5)
plt.xlabel("match_dist")

#gi = np.where(ff["category"] == 0)
#plt.hist(ff["match_dist"][gi],label="stars", bins=20, range=(0,100))

plt.legend()
plt.show()




gi = np.where(ff1["NN"] >0)[0]
plt.hist(ff1["sigma_r"][gi], label='random', bins=20, range=(0,10), density=True)


gi = np.where((ff0["category"] > 0) & (ff0["NN"] > 0))
plt.hist(ff0["sigma_r"][gi], label="training random", bins=20, range=(0,10), density=True, alpha=0.5)

gi = np.where((ff0["category"] == 0) & (ff0["NN"] >0))
plt.hist(ff0["sigma_r"][gi], label="training real", bins=20, range=(0,10), density=True, alpha=0.5)

plt.xlabel("sigma_r")
plt.ylabel("N")
#gi = np.where(ff["category"] == 0)
#plt.hist(ff["match_dist"][gi],label="stars", bins=20, range=(0,100))

plt.legend()
plt.show()




gi = np.where(ff1["NN"] >0)[0]
plt.hist(ff1["offset_sig"][gi], label='random', bins=20, range=(0,10), density=True)


gi = np.where((ff0["category"] > 0) & (ff0["NN"] > 0))
plt.hist(ff0["offset_sig"][gi], label="training random", bins=20, range=(0,10), density=True, alpha=0.5)

gi = np.where((ff0["category"] == 0) & (ff0["NN"] == 1))
plt.hist(ff0["offset_sig"][gi], label="training real", bins=20, range=(0,10), density=True, alpha=0.5)

plt.xlabel("offset_sig")
plt.ylabel("N")
#gi = np.where(ff["category"] == 0)
#plt.hist(ff["match_dist"][gi],label="stars", bins=20, range=(0,100))

plt.legend()
plt.show()





gi = np.where(ff1["NN"] >0)[0]
plt.hist(ff1["sigma_r"][gi], label='random', bins=20, range=(0,10), density=True)


gi = np.where((ff0["category"] > 0) & (ff0["NN"] > 0))
plt.hist(ff0["sigma_r"][gi], label="training random", bins=20, range=(0,10), density=True, alpha=0.5)

gi = np.where((ff0["category"] == 0) & (ff0["NN"] > 0))
plt.hist(ff0["sigma_r"][gi], label="training real", bins=20, range=(0,10), density=True, alpha=0.5)

plt.xlabel("sigma_r")

#gi = np.where(ff["category"] == 0)
#plt.hist(ff["match_dist"][gi],label="stars", bins=20, range=(0,100))

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

