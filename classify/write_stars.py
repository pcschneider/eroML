from joblib import dump, load
from eroML.classify import get_props, multidim_visualization, recovery, rescale
import numpy as np
from astropy.io import fits as pyfits


def write_one(clf, fn, ofn=None):
    props = clf.props
        
    
    #mfn = "random4classify_eFEDS.fits"
    Y = get_props(fn, prop_cols=props, category_column=None, pandas=True)
    c = clf.predict(Y)


    ff = pyfits.open(fn)
    fd = ff[1]
    column_formats = {col.name:col.format for col in ff[1].columns}


    cols = []
    for colname in fd.columns.names:        
        print(colname, column_formats[colname])
        col = pyfits.Column(name=colname , array=fd.data[colname], format=column_formats[colname])    
        cols.append(col)

    col = pyfits.Column(name="category", array=c, format="I")    
    cols.append(col)


    if ofn is not None:            
        hdu = pyfits.PrimaryHDU()    
        cc = pyfits.ColDefs(cols)    
        xx = pyfits.BinTableHDU.from_columns(cc)
        dst = len(xx.data[cc[0].name])
        hdul = pyfits.HDUList([hdu, xx])
        print("XXXXXX")
        if ofn is not None:
            hdul.writeto(ofn, overwrite=True)        
            print("datasets::prep_classify - Written ",dst," objects with ",len(cols)," properties to ",ofn)
    ff.close()   




clf = load('classify/svm2.joblib') 
write_one(clf, "major4classify_eFEDS.fits",     ofn = "major_eFEDS_classified.fits")
write_one(clf, "random4classify_eFEDS.fits",     ofn = "random_eFEDS_classified.fits")
