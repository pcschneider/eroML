
from astropy.io import fits as pyfits
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pyplot as plt
from sklearn import svm, tree
import numpy as np
import sklearn.utils as sku
from sklearn.linear_model import SGDClassifier
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer
from sklearn.metrics import accuracy_score, hinge_loss, balanced_accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.pipeline import Pipeline
from eroML.classify import get_props, multidim_visualization, recovery, rescale
from joblib import dump, load
from sklearn.preprocessing import FunctionTransformer

def fp(y_true, y_pred, pp = False): 
    a = confusion_matrix(y_true, y_pred)[0, 0]
    b = confusion_matrix(y_true, y_pred)[0, 1]
    c = confusion_matrix(y_true, y_pred)[1, 0]
    d = confusion_matrix(y_true, y_pred)[1, 1]
    if pp: print(a,b,c,d)
    return a + d - (b + c) * (b-c)**2


    
  

if __name__ == "__main__":
    
    props = ["bp_rp", "logFg","logFx", "pos","log_plx","skd"]
    props = ["sigma_r","match_dist", "bp_rp", "logFg","logFx", "log_plx","skd"]
    props = ["sigma_r","match_dist", "skd"]
    props = ["sigma_r","match_dist", "skd", "bp_rp", "logFg","logFx", "log_plx"]
    
    props = ["offset_sig", "expected_rnd","bp_rp", "logFg","logFx", "log_plx"]
    #props =  ["bp_rp", "logFxFg", "expected_rnd","offset_sig", "log_plx"]
    
    #props = ["bp_rp", "logFg", "logFxFg", "pos","log_plx","skd"]
    #props = ["pos","skd"]
    #X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")
    #X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")
    X, y = get_props("../../ero_data/training_eFEDS_clean.fits", prop_cols=props,category_column="category", pandas=True)
    y[y>0] = 1
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    
    #N = 10000
    #print(np.shape(X), np.shape(y))
    #idx = np.random.choice(range(len(y)), size=N)
    #X = X[idx]
    #y = y[idx]
    
    #clf = MLPClassifier(solver='adam', alpha=1e-1, hidden_layer_sizes=(6, 4), random_state=1, max_iter=10000)
    
    
<<<<<<< HEAD
=======
    i2 = len(np.where(np.logical_and(y>0, b==0))[0])
    print("Others as stars recovered: ",i2)


if __name__ == "__main__":

    #props = ["bp_rp", "log_FxFg", "offset_sig", "log_distance"]
    #props = ["log_FG", "log_FX", "bp_rp", "offset_sig", "log_distance"]
    #X, y = get_props("../merged_training.fits", prop_cols=props)

    props = ["logFx","logFg","pos","log_plx","bp_rp"]
    props = ["pos", "logFxFg","bp_rp","log_plx"]
    
    #props = ["bp_rp", "logFg","logFx", "pos","log_plx"]
>>>>>>> pcs
    
    #clf = svm.SVC(class_weight={1: 3}, probability=True)
<<<<<<< HEAD
    #clf = svm.SVC(C=30, kernel='rbf', probability=True, degree=3,class_weight={0: 0.24})
    #clf = svm.SVC(C=45, kernel='rbf', probability=True, degree=3,class_weight={0: 0.15})
    clf = svm.SVC(C=10, kernel='rbf', probability=True, degree=3,class_weight={0: 0.2})
    #clf = svm.SVC(C=30, kernel='poly', probability=True, degree=3,class_weight={0: 0.3})
    clf = svm.SVC(C=100, kernel='poly', probability=True, degree=2,class_weight={0: 0.18})
    # Greater C: less missclassification
    # Smaller C: More missclassification         
    #
=======
    clf = svm.SVC(C=5, kernel='poly', probability=True, degree=1,class_weight={0: 2})
>>>>>>> pcs
    #clf = PCA(n_components=2)
    #clf = tree.DecisionTreeClassifier()
    #clf = svm.SVC(kernel='linear', probability=True,class_weight={1: 3})
    #clf = SGDClassifier(loss='hinge')
    
    #clf.fit(X_train, y_train)
    ppl = Pipeline(steps=[( 'rescaler', FunctionTransformer(rescale)), ('svc', clf)])


    ppl.fit(X,  y)

    #scoring = {'AUC': 'roc_auc', 'Accuracy': make_scorer(accuracy_score),\
        #"hinge":make_scorer(hinge_loss, greater_is_better=False),\
        #"balanced": make_scorer(balanced_accuracy_score), 'fp': make_scorer(fp)}

    
    #gs = GridSearchCV(svm.SVC(C=1, kernel='poly', degree=3, class_weight={1: 1.0}, probability=True),
                    #param_grid={'C': np.logspace(-1.5,0.5,15), 'class_weight':[{1: w} for w in np.logspace(-1,1,15)]},
                    #scoring=scoring, refit='AUC', return_train_score=True, cv=5)

    #gs.fit(X_train, y_train)
    #results = gs.cv_results_
    #print(gs.best_estimator_)
    #print(gs.best_params_)
    ###print(gs.scorer_)


    #clf = gs.best_estimator_

    
    ##clf.fit(X_train, y_train)
    b = ppl.predict(X)
    recovery(y, b)
    
    dump(ppl, 'svm.joblib') 
    exit()
    
    print()
    #kernel = 1.0 * RBF(1.0)
    #gpc = GaussianProcessClassifier(kernel=kernel, random_state=0)
    #gpc.fit(X, y)
    #print(gpc.score(X, y))
    #c = gpc.predict(X)
    #recovery(y, c)
    
    #exit()
    #-------------------------------
    #-------------------------------
    #-------------------------------
    multidim_visualization(clf, X, y, names={i:props[i] for i in range(len(props))})
    #-------------------------------
    #-------------------------------
    #-------------------------------

    #plt.show()

    print(80*"=")
    #fn = "major_catalog.fits"
    #ofn = "major_proba.fits"
    
    #fn = "test_catalog.fits"
    #fn = "stellar_test_sources.fits"
    #ofn = "stellar_test_proba.fits"
    fn = "../../ero_data/major_eFEDS.fits"
    ofn = "major_proba.fits"

    #fn = "../../ero_data/random_eFEDS.fits"
    #ofn = "random_proba.fits"
    
    
    
    Y, eid, idx = get_props(fn, category_column=None, prop_cols=props, name_col="srcID", with_index=True)
    print("Shape of property array", np.shape(Y))
    print("eid", eid)
    #print(np.sum(np.isfinite(Y[0])))
    #print(np.sum(np.isfinite(Y[1])))

    c = clf.predict(Y)
    
    #multidim_visualization(clf, Y, c, names={i:props[i] for i in range(len(props))})
    
    pp0 = clf.predict_proba(Y)[::,0]
    pp1 = clf.predict_proba(Y)[::,1]
    #print(pp0+pp1)
    #plt.scatter(c, pp1)
    #plt.show()
    print("# recovered", len(c)-np.sum(c), " from ",np.shape(pp1))
    #print(len(np.where(pp1>0.5)[0]))
    
    
    ff = pyfits.open(fn)
    cols = list(ff[1].columns)
    NN = len(ff[1].data[props[0]])
    pp_zero = np.zeros(NN)
    print("# entries: ",NN, len(idx))
    pred = np.zeros(NN)
    pred[:]=np.nan
    pred[idx] = c
    
    col0 = pyfits.Column(name="predicted", array=pred, format="J")
    cols.append(col0)

    pp_zero[idx] = pp1
    col1 = pyfits.Column(name="proba", array=pp_zero, format="E")
    cols.append(col1)
    xx = pyfits.BinTableHDU.from_columns(cols)
    hdu = pyfits.PrimaryHDU()    
    hdul = pyfits.HDUList([hdu, xx])
    hdul.writeto(ofn, overwrite=True)
    
    print(80*"=")

    fn = "../../ero_data/random_eFEDS.fits"
    ofn = "random_proba.fits"

    #fn = "../../ero_data/random_eFEDS.fits"
    #ofn = "random_proba.fits"
    
    
    
    Y, eid, idx = get_props(fn, category_column=None, prop_cols=props, name_col="srcID", with_index=True)
    print("Shape of property array", np.shape(Y))
    #print("eid", eid)
    #print(np.sum(np.isfinite(Y[0])))
    #print(np.sum(np.isfinite(Y[1])))

    c = clf.predict(Y)
    pp1 = clf.predict_proba(Y)[::,1]

    print("# recovered", len(c)-np.sum(c), " from ",np.shape(pp1))
    #print(len(np.where(pp1>0.5)[0]))
    
    
    ff = pyfits.open(fn)
    cols = list(ff[1].columns)
    NN = len(ff[1].data[props[0]])
    pp_zero = np.zeros(NN)
    print("# entries: ",NN, len(idx))
    pred = np.zeros(NN)
    pred[:]=np.nan
    pred[idx] = c
    
    col0 = pyfits.Column(name="predicted", array=pred, format="J")
    cols.append(col0)

    pp_zero[idx] = pp1
    col1 = pyfits.Column(name="proba", array=pp_zero, format="E")
    cols.append(col1)
    xx = pyfits.BinTableHDU.from_columns(cols)
    hdu = pyfits.PrimaryHDU()    
    hdul = pyfits.HDUList([hdu, xx])
    hdul.writeto(ofn, overwrite=True)
