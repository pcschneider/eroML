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
from eroML.classify import multidim_visualization, scaler
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
from sklearn.linear_model import LogisticRegression
from eroML.classify import recovery

scoring = make_scorer(my_custom_loss_func, greater_is_better=False)

#oo = np.transpose([sigout, pos_off, skdens*3600, nth, cls])
fn = "../offs4.dat"
#fn = "../offs2.dat"
#fn = "simu.dat"
print("Reading property data ",fn)


#mfn = 'svc_linear.joblib'
mfn = 'svc.joblib'

    

dd = np.genfromtxt(fn, unpack=True)

#dd[1]*=3
gi = np.where(dd[3] == 1)[0] # Only nearest neighbour
#gi = np.where((dd[3] == 1) & (dd[1] < 3*dd[0]))[0] # Only nearest neighbour
gi = np.where((dd[3] == 1) & (dd[1] < 3*dd[0]))[0] # Only nearest neighbour
X = np.transpose(dd[0:3,gi])
y = dd[4][gi]
y[y>0] = 1
print("#true stars: ",len(y) - np.sum(y))

N = len(gi)
print("X.shape: ",np.shape(X), "N: ",N)
gg = 1 / X.shape[1] / X.var()
print("Var: ",X.var(), " gamma (for 'scale') = ", gg)
# gamma = 1 / (n_features * X.var() if gamma = 'scale (default)'
      

gi = np.where(X[::,1] < 0.7*X[::,0])[0]
sw_train = np.ones(len(y))
sw_train[gi] = 2.3
print("#Weighted: ",len(gi), " of ",N)
train_index = np.random.choice(np.arange(N), N, replace=False)
#print(train_index[0:10])
##clf = svm.SVC(C=65, kernel='linear', probability=True, degree=2, class_weight={0: 0.33}, cache_size=2000)
#clf = svm.SVC(C=65, probability=False, kernel='poly', degree=3 ,class_weight={0: 1.0}, gamma=3e-4)
#clf = svm.SVC(C=50, kernel='poly', degree=3)
clf = Pipeline(steps=[('scale',  FunctionTransformer(func=scaler, kw_args={"factor":1, "axis":1})), ('clf', svm.SVC(C=0.5, class_weight={0: 0.17}, cache_size=2000))])
#weight: 1.15

#clf = Pipeline(steps=[('poly', PolynomialFeatures(4)), ('clf', svm.LinearSVC(C=1,class_weight={0: 3.2}, max_iter=10000, dual=False))])
#clf = LogisticRegression(C=1, penalty='l2', solver='saga', multi_class='multinomial', max_iter=10000,class_weight={0: 1.5})

print("start, Ntrain=",len(train_index))       
t0 = time.perf_counter()
clf.fit(X[train_index],  y[train_index])#, clf__sample_weight=sw_train[train_index])
#clf.fit(X,  y)#, sample_weight=sw_train)
print("fitted (in %5.2f s)" % float(time.perf_counter()-t0))
clf.X_ = X
clf.y_ = y
clf.filename_ = fn
dump(clf, mfn)

b = clf.predict(X)
recovery(y, b)

