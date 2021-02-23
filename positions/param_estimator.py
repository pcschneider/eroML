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


def recovery(y, b):
    print("Stars in training set: ",len(np.where(y==0)[0]))
    print("Others in training set: ",len(np.where(y==1)[0]))
    print("Random in training set: ",len(np.where(y==2)[0]))

    print("Stars predicted: ",len(np.where(b==0)[0]))
    print("Others predicted: ",len(np.where(b>0)[0]))

    
    i0 = len(np.where(np.logical_and(y==0, b==0))[0])
    print("Stars as stellar recovered: ",i0)
    
    i1 = len(np.where(np.logical_and(y>0, b>0))[0])
    print("Others as others recovered: ",i1)
    
    i3 = len(np.where(np.logical_and(y==0, b>0))[0])
    print("Stars as others recovered: ",i3)
    
    i2 = len(np.where(np.logical_and(y>0, b==0))[0])
    print("Others as stars recovered: ",i2)
    return i0, i1, i3, i2



if __name__ == "__main__":
    dd = np.genfromtxt("../offs3.dat", unpack=True)
    #dd = np.genfromtxt("simu.dat", unpack=True)

    #dd[1]*=3
    gi = np.where(dd[3] == 1)[0] # Only nearest neighbour
    X = np.transpose(dd[0:3,gi])
    y = dd[4][gi]
    y[y>0] = 1

    print(np.shape(X), len(y))
    Nstars = len(y) - np.sum(y)
    print("#true stars: ",Nstars)
    
    gg = 1 / X.shape[1] / X.var()
    print("Var: ",X.var(), " gamma (for 'scale') = ", gg)
    
    # Sample Weights
    gi = np.where(X[::,1] < 0.7*X[::,0])[0]
    sw_train = np.ones(len(y))
    #sw_train[gi] = 2.1
    #print("#Weighted: ",len(gi))
    
    #X_train, X_test, y_train, y_test = train_test_split(
    #X, y, test_size=0.25, random_state=42)
    N = len(X)
    train_index = np.random.choice(np.arange(N), N, replace=False)
    print("Using ",len(train_index), " of ",N," samples for training.")
    Ntotal = len(train_index)
    Nstars = len(train_index) - np.sum(y[train_index])
    print(" so that there are ",Nstars," stars")
    
    poly = PolynomialFeatures(3)
    X_features = poly.fit_transform(X)
    print("featured", np.shape(X_features[train_index]))
    
    oo = open("params2.txt","w")
    oo.write("# C w i0 i1 i2 i3 Nstars Ntotal\n")
    oo.write("#  0 - C\n")
    oo.write("#  1 - class weight for class 0 \n")
    oo.write("#  2 - i0: Stars as stellar recovered\n")
    oo.write("#  3 - i1: Others as others recovered\n")
    oo.write("#  4 - i2: Stars as others recovered\n")
    oo.write("#  5 - i3: Others as stars recovered\n")
    oo.write("#  6 - Nstars: Stars in sample\n")
    oo.write("#  7 - Ntotal: Sample size\n")
  
    
    print("start")       
    t0 = time.perf_counter()

    for c in np.logspace(-2,2,15):
    #for c in np.linspace(0.01,0.12,3):
        for w in np.linspace(0.2,2.5,50):
            print("C=",c, "weight for class=0: ",w)
            #clf = Pipeline(steps=[('poly', poly), ('clf', svm.LinearSVC(C=c,class_weight={0: w}, max_iter=10000, dual=False))])
            clf = svm.SVC(C=c, kernel='linear', probability=True, degree=3,class_weight={0: w})
            clf = Pipeline(steps=[('poly', PolynomialFeatures(3)), ('clf', svm.LinearSVC(C=c,class_weight={0: w}, max_iter=10000, dual=False))])             
            #clf = svm.NuSVC(nu=c, kernel='linear', probability=True, degree=3,class_weight={0: w})

            #clf = Pipeline(steps=[('poly', poly), ('clf', svm.LinearSVC(C=c,class_weight={0: w}, max_iter=10000, dual=False))])
            #clf = svm.SVC(C=c, kernel='rbf', probability=True, degree=3,class_weight={0: w})
            clf.fit(X,  y)
            b = clf.predict(X)
            i0,i1,i2,i3 = recovery(y, b)
            oo.write("%f %f %i %i %i %i %i %i \n" % (c, w, i0, i1, i2, i3, Nstars, N))
            oo.flush()
            print()
    oo.close()            
            
            
