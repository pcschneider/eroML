from joblib import dump, load
from eroML.classify import get_props, multidim_visualization, recovery, rescale
import numpy as np
from astropy.io import fits as pyfits



if __name__ == "__main__":
    props = ["offset_sig", "expected_rnd","bp_rp", "logFg","logFx", "log_plx"]
    X, y = get_props("../../ero_data/training_eFEDS_clean.fits", prop_cols=props,category_column="category", pandas=True)
    y[y>0] = 1
    clf = load('svm.joblib') 
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
    
    
    
    Y, eid, idx = get_props(fn, category_column=None, prop_cols=props, name_col="srcID", with_index=True, pandas=True)
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
