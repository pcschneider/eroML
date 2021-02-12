from joblib import dump, load
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
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import make_scorer 
from scipy.stats import uniform
from sklearn.pipeline import make_pipeline, Pipeline
from sklearn.preprocessing import StandardScaler, FunctionTransformer, PolynomialFeatures
from sklearn.model_selection import KFold
from sklearn.kernel_approximation import PolynomialCountSketch
import time
from eroML.classify import my_custom_loss_func

scoring = make_scorer(my_custom_loss_func, greater_is_better=False)

#oo = np.transpose([sigout, pos_off, skdens*3600, nth, cls])
dd = np.genfromtxt("../offs.dat", unpack=True)
dd = np.genfromtxt("simu.dat", unpack=True)

#dd[1]*=3
gi = np.where(dd[3] == 1)[0] # Only nearest neighbour
X = np.transpose(dd[0:3,gi])
gg = 1 / X.shape[1] / X.var()
print("Var: ",X.var(), " gamma (for 'scale') = ", gg)

# gamma = 1 / (n_features * X.var() if gamma = 'scale (default)'
      
      
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

gi = np.where(X[::,1] < 0.7*X[::,0])[0]
sw_train = np.ones(len(y))
print("#Weighted: ",len(gi))
#sw_train[gi] = 2.1

print("#true stars: ",len(y) - np.sum(y))
#X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")

#clf = svm.SVC(class_weight={1: 3}, probability=True)
#clf = svm.SVC(C=100, probability=True, kernel='poly', degree=2,class_weight={1: 10, 0:0.5}, gamma=50)
clf = SGDClassifier(loss="hinge", alpha=0.01, max_iter=200)

#https://scikit-learn.org/stable/tutorial/statistical_inference/putting_together.html



parameters = {'clf__C':np.linspace(30,80,10), 'clf__class_weight':[{1:x} for x in np.linspace(2,3,10)], 'clf__gamma':np.logspace(np.log10(gg)-1, np.log10(gg)+1,20)}
#print("parameters: ",parameters)
sv = svm.SVC(probability=True, kernel='poly', degree=2)
#grid = GridSearchCV(sv, parameters, scoring=scoring)#, fit_params={'sample_weight': sw_train})

##distributions = dict(C=uniform(loc=10, scale=90), class_weight=[{1:x} for x in np.linspace(1,5,15)])#, gamma=np.logspace(-7,-5,20))
##distributions = dict(C=uniform(loc=10, scale=90), class_weight={0:1, 1:uniform(loc=1, scale=3)})
#print("distributions: ",distributions)

#grid = RandomizedSearchC
#grid = GridSearchCV(ppl, parameters, scoring=scoring)

#grid.fit(X, y, clf__sample_weight= sw_train)
#print(grid.cv_results_.keys())
#for k in grid.cv_results_.keys():
    #print(k, grid.cv_results_[k])
#print()
#print(grid.best_params_)
#print()
#print(grid.best_index_)
#print("optimized")


clf = svm.SVC(C=65, probability=True, kernel='poly', degree=2,class_weight={1: 2.3})#, gamma=50)

#clf = svm.SVC(C=65, probability=True, kernel='poly', degree=3,class_weight={1: 3.7}, gamma=0.0008*25)

#exit()
#clf = svm.SVC(C=25, probability=True, kernel='poly', degree=3, class_weight={1: 2.7}, gamma=2.77e-4)
#clf = SGDClassifier(loss="hinge", alpha=0.01, max_iter=200)


#clf = PCA(n_components=2)
#clf = tree.DecisionTreeClassifier()
clf = svm.SVC(kernel='linear', probability=True,class_weight={1: 3})

print(dir(clf))
#print("gamma:", clf.gamma, "coeff0: ",clf.coef0)
#print("class_weight",clf.class_weight, "degree", clf.degree)
#clf.fit(X, y)


#ps = PolynomialCountSketch(degree=3, random_state=1)


#X_features = poly.fit_transform(X)
#print("shape:",np.shape(X_features))
#X_old = X
#X = X_features
N = len(X)
train_index = np.random.choice(np.arange(N), N)
print(len(train_index))


sgd = SGDClassifier(loss='hinge')


print("start")       
t0 = time.perf_counter()

#svc = svm.SVC(C=65, probability=False, kernel='poly', degree=3,class_weight={1: 2.3})
#svc.fit(X[train_index], y[train_index], sample_weight=sw_train[train_index])
#print("x - x", time.perf_counter()-t0)

poly = PolynomialFeatures(3)
X_features = poly.fit_transform(X)
print("featured", np.shape(X_features[train_index]))
clf = Pipeline(steps=[('poly', poly), ('clf', svm.LinearSVC(C=65,class_weight={0: 2}, max_iter=10000, dual=False))])

parameters = {'clf__C':np.linspace(5,50,4), 'clf__class_weight':[{0:x} for x in np.linspace(0.2,5,3)]}
grid = GridSearchCV(clf, parameters, scoring=scoring, cv=3, n_jobs=-1)
grid.fit(X[train_index], y[train_index], clf__sample_weight=sw_train[train_index])

#print(grid.clf_results_.keys())
for k in grid.cv_results_.keys():
    print(k, grid.cv_results_[k])
#print()
print(grid.best_params_)
print()
print(grid.best_index_)
print("optimized")

#clf = svm.LinearSVC(C=65,class_weight={1: 2.3})

#clf = MLPClassifier(hidden_layer_sizes=(30,10))


clf = grid.best_estimator_

grid.X_ = X
grid.y_ = y

dump(grid, 'grid.joblib') 

#clf.fit(X[train_index], y[train_index], clf__sample_weight=sw_train[train_index])
print("xxx", time.perf_counter()-t0)    
b = clf.predict(X)
#pp = clf.predict_proba(X)
#print(pp.shape)
#plt.scatter(pp[::,0], b, color='r')
#plt.scatter(pp[::,1], b, color='b')
#plt.show()
#recovery(y, b)
print("yyy steps:", clf["clf"].n_iter_)
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
#multidim_visualization(clf, X_features, y, names={i:"param"+str(i) for i in range(np.shape(X_features)[1])})
#multidim_visualization(clf, X, y, names={1:"offset", 0:"sig"})
#exit()    
