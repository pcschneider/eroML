import pandas as pd
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm


def rescale(X):
    X["offset_sig"] = np.log10(1.-norm.cdf(X["offset_sig"]))
    gi = np.where(-2*np.isfinite(X["offset_sig"])==False)[0]
    X["offset_sig"][gi] = 100
    return X

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


def get_props(ifn, category_column="class", prop_cols=[], name_col=None, filter_column=None, filter_val=None, with_index=False, pandas=False):
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
    
    if prop_cols == None:
        pp = []
        for c in ff[1].columns:
            if c.name == "srcID": continue
            if c.name == category_column: continue
            pp.append(c.name)
        prop_cols = pp
        print("Using prop_cols=",prop_cols)
    
    N = len(ff[1].data[prop_cols[0]])
    use = np.ones(N).astype(bool)
    
    if filter_column is not None:
        no_good = np.where(ff[1].data[filter_column] != filter_val)[0]
        print("Filtering for \'%s == %s\' leaves %i from %i entries in sample." %(str(filter_column), str(filter_val), N-len(no_good), N))
        use[no_good] = 0
        
    idx = np.arange(N)
    
    
    names = []
    formats = []
    rcols = []
    for pc in prop_cols:
        #print("Reading col \'",pc,"\'")
        rcols.append(pc)
        tmp = ff[1].data[pc].flatten()
        gi = np.where(np.isfinite(tmp)!=True)[0]
        use[gi] = 0
        #print(np.shape(tmp), np.sum(use))
        props.append(tmp.astype(float))
        fmt = tmp.dtype
        #print(pc, "format: ",fmt)
        #fmt = np.float64
        names.append(pc)
        formats.append(fmt)
        
    print("Read cols",", ".join(rcols))    
    props = np.array(props)    
    
    props = props.T[use, ::]
    
    if pandas:
        Y = pd.DataFrame(data=props, columns=names)
        props = Y
        props = Y[prop_cols]
    
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
