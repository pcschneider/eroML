import matplotlib.pyplot as plt
from sklearn import svm, tree
import numpy as np
from joblib import dump, load
from eroML.classify import recovery

fn = "svc.joblib"
#fn = "svc_uniform.joblib"
#fn = "svc_linear.joblib"
#fn = "svc_LR.joblib"

print("Reading ",fn)

clf = load(fn)
print("    trained with ",clf.filename_)
b = clf.predict(clf.X_)
recovery(clf.y_, b)
print()

fn = "../offs2.dat"
print("Evaluating: ",fn)
dd = np.genfromtxt(fn, unpack=True)
gi = np.where(dd[3] == 1)[0] # Only nearest neighbour
X = np.transpose(dd[0:3,gi])
y = dd[4][gi]
y[y>0] = 1
b = clf.predict(X)
recovery(y, b)
print()


fn = "simu.dat"
print("Evaluating: ",fn)
dd = np.genfromtxt(fn, unpack=True)
gi = np.where(dd[3] == 1)[0] # Only nearest neighbour
X = np.transpose(dd[0:3,gi])
y = dd[4][gi]
y[y>0] = 1
b = clf.predict(X)
recovery(y, b)

