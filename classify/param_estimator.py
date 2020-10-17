
from astropy.io import fits as pyfits
from astropy.coordinates import SkyCoord
from astropy import units as u
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
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import make_scorer
from sklearn.metrics import accuracy_score, hinge_loss, balanced_accuracy_score
from sklearn.metrics import confusion_matrix

def get_props(ifn, category_column="class", prop_cols=[], name_col=None, filter_column=None, filter_val=None, with_index=False):
    """
      Create np.array for use with scikit-learn
      
      Parameters
      ----------
       ifn : str, filename
       name_col : str, use this fits-column for the name/identifier for each entry
       with_index : booelan
       filter_column, filter_val : str, use-specific
         Filter column with column-name `filter_column` for `filter_val`
       
     Returns
     -------
       props : if `prop_cols` is populated
       props, category : if `prop_cols` is populated *and* `category_column` is not None
       
    """
    print("Reading ",ifn)
    ff = pyfits.open(ifn)
    props = []
    
    N = len(ff[1].data[prop_cols[0]])
    use = np.ones(N).astype(bool)
    
    if filter_column is not None:
        no_good = np.where(ff[1].data[filter_column] != filter_val)[0]
        print("Filtering for \'%s == %s\' leaves %i from %i entries in sample." %(str(filter_column), str(filter_val), N-len(no_good), N))
        use[no_good] = 0
        
    idx = np.arange(N)
    
    
    
    for pc in prop_cols:
        print("Reading col \'",pc,"\'")
        tmp = ff[1].data[pc].flatten()
        gi = np.where(np.isfinite(tmp)!=True)[0]
        use[gi] = 0
        #print(np.shape(tmp), np.sum(use))
        props.append(tmp)
        
    props = np.array(props)    
    #print(np.shape(props))
    props = props.T[use, ::]
    print("    final property array.shape: ", np.shape(props)) 
    
    
    
    if category_column is not None:
        cats = np.unique(ff[1].data[category_column][use])
        for cc in cats:
            print("category:", cc, len(np.where(ff[1].data[category_column][use] == cc)[0]))
        y = ff[1].data[category_column][use]
        if name_col is not None: 
            if with_index:
                return props, y, ff[1].data[name_col][use], idx[use]
            else: props, y, ff[1].data[name_col][use]
        else:
            if with_index:
                return props, y, idx[use]
            else:
                return props, y
    else:
        if name_col is not None: 
            if with_index:
                return props, ff[1].data[name_col][use], idx[use]
            else:
                return props, ff[1].data[name_col][use]
        if with_index:
                return props, idx[use]    
        return props


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
    props = ["bp_rp", "logFg","logFx", "pos","log_plx","skd"]
    #props = ["pos","skd"]
    #X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")
    #X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")
    X, y = get_props("../../ero_data/training_eFEDS.fits", prop_cols=props,category_column="category")
    y[y>0] = 1
    
    oo = open("params2.txt","w")
    oo.write("# C w i0 i1 i2 i3\n")
    oo.write("# i0: Stars as stellar recovered\n")
    oo.write("# i1: Others as others recovered\n")
    oo.write("# i2: Stars as others recovered\n")
    oo.write("# i3: Others as stars recovered\n")
    for c in [20,25,30,35]:#np.logspace(-1,1,10):
        for w in [0.18,0.2,0.22,0.24]:#np.logspace(-1,1,10):
            print("C=",c, "weight for class=0: ",w)
            clf = svm.SVC(C=c, kernel='rbf', probability=True, degree=3,class_weight={0: w})
            clf.fit(X,  y)
            b = clf.predict(X)
            i0,i1,i2,i3 = recovery(y, b)
            oo.write("%f %f %i %i %i %i\n" % (c, w, i0, i1, i2, i3))
            oo.flush()
            print()
    oo.close()            
            
            
