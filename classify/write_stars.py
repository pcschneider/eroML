from joblib import dump, load
from eroML.classify import get_props, multidim_visualization, recovery, rescale
import numpy as np
from astropy.io import fits as pyfits


clf = load('classify/svm.joblib') 
props = clf.props
    
mfn = "major4classify_eFEDS.fits"
Y = get_props(mfn, prop_cols=props, category_column=None, pandas=True)
c = clf.predict(Y)


ff = pyfits.open(mfn)
fd = ff[1]
column_formats = {col.name:col.format for col in ff[1].columns}

ofn = "major_eFEDS_classified.fits"
cols = []
for colname in fd.columns.names:        
    print(colname, column_formats[colname])
    col = pyfits.Column(name=colname , array=fd.data[colname], format=column_formats[colname])    
    cols.append(col)

col = pyfits.Column(name="category", array=c, format="I")    
cols.append(col)


        
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


