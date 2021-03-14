from astropy.io import fits as pyfits
import numpy as np
from eroML.ensemble import multi_fits_support

def gen_assocs(fn, ero_col, gaia_col, ero_mapper=lambda x: x, gaia_mapper=lambda x: x):
    ff = pyfits.open(fn)
    fd = ff[1].data
    col0 = ero_mapper(fd[ero_col])
    col1 = gaia_mapper(fd[gaia_col])
    assocs = [str(a).strip()+"&"+str(b).strip() for a, b in zip(col0, col1)]
    return assocs
    
def get_cols(fn, ext=1):
    ff = pyfits.open(fn)
    fd = ff[ext].data
    column_formats = {col.name:col.format for col in ff[ext].columns}
    ff.close()
    return column_formats

def get_props(fn):
    ff = pyfits.open(fn)
    ext = 1
    fd = ff[ext].data
    column_formats = {col.name:col.format for col in ff[ext].columns}
    ret = {}
    for col in column_formats.keys():
        ret[col] = fd[col]
    ff.close()
    return ret

if __name__== "__main__":
    ass0 = gen_assocs("eFEDS_HamStar.fits", "ero_ID", "gaia_ID")
    cf = get_cols("../ero_data/efeds_c001_V3_main_HamStar_internal2.fits")
    ass1 = gen_assocs("../ero_data/efeds_c001_V3_main_HamStar_internal2.fits", "ero_ID", "gaia_ID")
    ass2 = gen_assocs("../ero_data/nway2.fits", "ero_ID", "gaia_ID")    
    props = []
    for fn in ["eFEDS_HamStar.fits", "../ero_data/efeds_c001_V3_main_HamStar_internal2.fits", "../ero_data/nway2.fits"]:
        props.append(get_props(fn))
    
    mapper = {0:"SVM", 1:"Bayes", 2:"NWAY"}
    
    allass = [ass0, ass1, ass2]
    ass = np.unique(ass0+ass1+ass2)
    required = np.ones(len(ass))
    
    
    dct = {}
    for col in cf.keys():
        #print(col, cf[col])
        if "A" in cf[col]: dt = '<U32'
        elif "D" in cf[col]: dt = float
        elif "I" in cf[col]: dt = int
        dct[col] = np.empty(len(ass), dtype=dt)
        
    for m in mapper.values():
        dct[m] = np.zeros(len(ass))
        dct[m+"_ij"] = np.zeros(len(ass))
    
    #print(dct.keys())
    for j in range(len(allass)):
        print("Working on associations in ",mapper[j])
        a = allass[j]
        ii = np.in1d(ass, a) # ass in a
        print(ii, len(ii), np.sum(ii))
        setti = ii.astype(bool) & required.astype(bool)
        print(setti, np.sum(setti))
        
        xy, x_ind, y_ind = np.intersect1d(ass, a, return_indices=True)            
        if j < 2:
            dct[mapper[j]][x_ind] = props[j]["p_stellar"][y_ind]
            dct[mapper[j]+"_ij"][x_ind] = props[j]["p_ij"][y_ind]
        else:
            dct[mapper[j]][x_ind] = props[j]["p_any"][y_ind]
            dct[mapper[j]+"_ij"][x_ind] = props[j]["p_ij"][y_ind]
        
        
        if np.sum(setti) == 0: continue
        required[setti] = False
        idx = np.arange(len(ass))[setti]
        xy, x_ind, y_ind = np.intersect1d(ass[setti], a, return_indices=True)
        print(len(xy), x_ind, y_ind)
        print(ass[idx][x_ind[10]], props[j]["ero_ID"][y_ind[10]])
        for col in cf.keys():    
            dct[col][idx[x_ind]] = props[j][col][y_ind]
            
        
            #print(dct[col][y_ind])
        print()
    
    
    
    
    cols = []
    for col in cf.keys():
        c = pyfits.Column(name=col, array=dct[col], format=cf[col])    
        cols.append(c)
    for col in mapper.values():
        c = pyfits.Column(name=col, array=dct[col], format="D")    
        cols.append(c)
        c = pyfits.Column(name=col+"_ij", array=dct[col+"_ij"], format="D")    
        cols.append(c)
        
    ofn = "ero_master.fits"
    overwrite = True 
    verbose=1
    hdu = pyfits.PrimaryHDU()    
    cc = pyfits.ColDefs(cols)    
    xx = pyfits.BinTableHDU.from_columns(cc)
    dst = len(xx.data[cc[0].name])
    hdul = pyfits.HDUList([hdu, xx])
    print("XXXXXX")
    if ofn is not None:
        hdul.writeto(ofn, overwrite=overwrite)        
        if verbose>0:
            print("datasets::prep_classify - Written ",dst," objects with ",len(cols)," properties to ",ofn)
    #ff.close()  
    
    #print(dct)
    exit()
    #for a in 
    
    print(len(ass0),    len(ass1), len(np.intersect1d(ass0, ass1)))
    ii = np.in1d(ass0, ass1)
    print(len(ii)) # ass0 in ass1
    print(np.array(ass0)[~ii])
    #print(ass0[1], ass1[5])
    exit()
    dct = props2asso("../../eroML/major_eFEDS_classified.fits", assocs=ass0, props=["srcID", "srcID_NN"])
    print(len(ass0))
    print(dct)

