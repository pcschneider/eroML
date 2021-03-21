from joblib import dump, load
from eroML.classify import get_props
from astropy.io import fits as pyfits
import numpy as np
from eroML.positions import calc_sigma_from_RADEC_ERR

fn = "svc.joblib"
#fn = "svc_uniform.joblib"
#fn = "svc_linear.joblib"
#fn = "svc_LR.joblib"

print("Reading ",fn)
clf = load(fn)


#props = ["offset_sig",]
#props =  ["bp_rp", "logFxFg", "expected_rnd","offset_sig", "log_plx"]

#props = ["bp_rp", "logFg", "logFxFg", "pos","log_plx","skd"]

#X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")
#X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")

props = ["match_dist","RADEC_ERR","eligible_sky_density"]


fn = "../../ero_data/merged_major_eFEDS_EDR3.fits"
fn = "../../ero_data/merged_major_eFEDS_EDR3_HamStar.fits"
#fn = "../../ero_data/merged_random_eFEDS_EDR3.fits"
ofn = "training_IDs4.txt"
ff = pyfits.open(fn)
fd = ff[1].data
X = np.transpose([calc_sigma_from_RADEC_ERR(fd["RADEC_ERR"]), fd["match_dist"], fd["eligible_sky_density"]])
print(np.shape(X))

y = clf.predict(X)

gi = np.where((X[::,0] > 7) | (X[::,1] > 10))[0] # Only nearest neighbour
print("Filter: ",len(gi))
y[gi] = 1
#p = clf.predict_proba(X)
gi = np.where(y==0)[0]
#print(fd["srcID"][gi])
print("SVM: ",len(gi))

np.savetxt(ofn, fd["srcID"][gi], fmt="%s")

dd = np.genfromtxt("../Bayes_eFEDS_good_pos.txt")
#dd = np.genfromtxt("../eFEDS_final_training_set.txt")
print("Sebastian",len(dd))

a = fd["srcID"][gi]
b = np.array([str("ML%05i" % d) for d in dd])
print(a.dtype)
print(b.dtype)

shared = np.intersect1d(a, b)
print("both: ",len(shared))
#print(shared)

