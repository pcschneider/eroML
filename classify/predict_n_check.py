from joblib import dump, load
from eroML.classify import get_props, multidim_visualization, recovery, rescale
import numpy as np
from astropy.io import fits as pyfits



if __name__ == "__main__":
    clf = load('classify/svm_tmp.joblib') 
    props = clf.props
    
    rfn = "random4classify_eFEDS_HamStar.fits"
    Y = get_props(rfn, prop_cols=props, category_column=None, pandas=True)
    c = clf.predict(Y)
    Y = get_props(rfn, prop_cols=props+["NN"], category_column=None, pandas=True)
    
    ff = pyfits.open(rfn)
    fd = ff[1].data
    srcID = [x[0:7] for x in fd["srcID"]]
    #print(srcID)
    try:
        gi = np.where((Y["NN"] == 1) & (Y["FxFg"]<-1))[0]
    except:
        gi = np.where(Y["NN"] == 1)[0]
        
    
    si = np.where(c==0)[0]
    ui = np.unique(fd["original_srcID"][si])
    print("random stars: ",len(si), "(NN=1",len(np.where(c[gi]==0)[0]),",", len(ui),")")
    
    rfn = "major4classify_eFEDS_HamStar.fits"
    Y = get_props(rfn, prop_cols=props, category_column=None, pandas=True)
    c = clf.predict(Y)
    Y = get_props(rfn, prop_cols=props+["NN"], category_column=None, pandas=True)

    ff = pyfits.open(rfn)
    fd = ff[1].data

    si = fd["original_srcID"][np.where(c==0)[0]]
    print("real stars: ",len(np.where(c==0)[0]), "(",len(np.unique(si)),")")
    #multidim_visualization(clf, Y, c, names={i:props[i] for i in range(len(props))})

 
    fn = clf.filename 
    

    #props = ["offset_sig", "expected_rnd","bp_rp", "logFg","logFx", "log_plx"]
    #X, y = get_props("../../ero_data/training_eFEDS_clean.fits", prop_cols=props,category_column="category", pandas=True)
    X, y = get_props(fn, prop_cols=props,category_column="category", pandas=True)
    y[y>0] = 1
    
   #-------------------------------
    #-------------------------------
    multidim_visualization(clf, X, y, names={i:props[i] for i in range(len(props))})
    #-------------------------------
    #-------------------------------
    #-------------------------------

    #plt.show()

    print(80*"=")

    Y = get_props(rfn, prop_cols=props, category_column=None, pandas=True)
    c = clf.predict(Y)
    multidim_visualization(clf, Y, c, names={i:props[i] for i in range(len(props))})
    
    exit()
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
