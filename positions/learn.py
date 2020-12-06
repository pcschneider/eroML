import matplotlib.pyplot as plt
from sklearn import svm, tree
import numpy as np
import sklearn.utils as sku
from sklearn.linear_model import SGDClassifier
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.text import Text
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from eroML.classify import multidim_visualization
from sklearn.model_selection import GridSearchCV


def my_custom_loss_func(y_true, y_pred):
     
     return np.log1p(diff)

#oo = np.transpose([sigout, pos_off, skdens*3600, nth, cls])
dd = np.genfromtxt("../offs.dat", unpack=True)
#dd[1]*=3
gi = np.where(dd[3] == 1)[0] # Only nearest neighbour
X = np.transpose(dd[0:2,gi])
#for i in range(X.shape[1]):
    #std = np.std(X[::,i])
    #print(i, std)
    #X[::,i]/=std
##X[::,1]*=30

#X[::,0]*=2
#X[::,1]*=10

#for i in range(3):
    #plt.hist(X[::,i])
    #plt.title(str(i))
    #plt.show()
print(np.shape(X))
y = dd[4][gi]
y[y>0] = 1
#X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")

#clf = svm.SVC(class_weight={1: 3}, probability=True)
clf = svm.SVC(C=100, probability=True, kernel='poly', degree=2,class_weight={1: 10, 0:0.5}, gamma=50)


parameters = {'C':np.linspace(1,100,3), 'class_weight':[{1:2}, {1:2.5},{1:3}]}

sv = svm.SVC(probability=True, kernel='poly', degree=2)
grid = GridSearchCV(sv, parameters)
grid.fit(X, y)
print("optimized")

clf = svm.SVC(C=25, probability=True, kernel='poly', degree=3, class_weight={1: 2.7}, gamma=2.77e-4)


#clf = PCA(n_components=2)
#clf = tree.DecisionTreeClassifier()
#clf = svm.SVC(kernel='linear', probability=True,class_weight={1: 3})
#clf = SGDClassifier(loss='hinge')
gi = np.where(X[::,1] < 0.7*X[::,0])[0]
sw = np.ones(len(y))
print("#Weighted: ",len(gi))
sw[gi] = 2.5
print(dir(clf))
print("gamma:", clf.gamma, "coeff0: ",clf.coef0)
print("class_weight",clf.class_weight, "degree", clf.degree)
#clf.fit(X, y)

clf.fit(X, y, sample_weight=sw)
b = clf.predict(X)
pp = clf.predict_proba(X)
#print(pp.shape)
#plt.scatter(pp[::,0], b, color='r')
#plt.scatter(pp[::,1], b, color='b')
#plt.show()
#recovery(y, b)
print()
#kernel = 1.0 * RBF(1.0)
#gpc = GaussianProcessClassifier(kernel=kernel, random_state=0)
#gpc.fit(X, y)
#print(gpc.score(X, y))
#c = gpc.predict(X)
#recovery(y, c)
gi = np.where(b == y)[0]
print("Correct: ",len(gi), " of ",len(b))

gi = np.where((b == 0) & (y==0))[0]
print("star as star: ",len(gi), " of ",len(np.where(y==0)[0]))

gi = np.where((b == 1) & (y==1))[0]
print("random as random: ",len(gi), " of ",len(np.where(y==1)[0]))

gi = np.where((b == 0) & (y==1))[0]
print("random as star: ",len(gi), " of ",len(b))
gi = np.where((b == 1) & (y==0))[0]
print("star as random: ",len(gi), " of ",len(b))

gi = np.where(b==0)[0]
print("Objects as star classified: ",len(gi))
#exit()
multidim_visualization(clf, X, y, names={1:"offset", 2:"sky_dens", 0:"sig"})
#multidim_visualization(clf, X, y, names={1:"offset", 0:"sig"})
#exit()    
