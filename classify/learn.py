
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

def multidim_visualization(clf, X, y, names=None, dims=[0,1]):    
    
    def draw_contours():
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ## create grid to evaluate model
        xx = np.linspace(xlim[0], xlim[1], 30)
        yy = np.linspace(ylim[0], ylim[1], 30)
        YY, XX = np.meshgrid(yy, xx)
        #print("lims: ",xlim, ylim)
        if np.shape(X)[1]==2:
            xy = np.vstack([XX.ravel(), YY.ravel()]).T
        else:
            
            tmp = [XX.ravel(), YY.ravel()]
            for i in range(np.shape(X)[1]-2):
                #zz = np.median(X[i+2,::])
                zz = sliders[i].val
                tmp.append(900*[zz])
            xy = np.vstack(tmp).T
            
        #Z = clf.decision_function(xy).reshape(XX.shape)
        Z = clf.predict_proba(xy)[::,0].reshape(XX.shape)
        ## plot decision boundary and margins
        global con
        #con = ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5,
                #linestyles=['--', '-', '--'])
        
        con = ax.contour(XX, YY, Z, colors='b', levels=[0.1, 0.5, 0.9], alpha=0.5,
                linestyles=['--', '-', '--'])
        return con
     
    def update_other_dim_values(val):
        global con        
        if con is not None:
            for coll in con.collections: 
                ax.collections.remove(coll) 
        con = draw_contours()
            
        fig.canvas.draw_idle()
    
    def value_range4dim(dim):
        mn, mx = np.min(X[::,dim]), np.max(X[::,dim])
        mean, med = np.mean(X[::,dim]), np.median(X[::,dim])
        return mn, med, mx
    
    def gen_slider(pos, dim, scaling=1.1):
        a,b,c = value_range4dim(dim)
        #print("pos", pos, dim, len(axsliders), names)
        axsliders[pos].clear()
        sliders[pos].__init__(axsliders[pos], names[dim], valmin=a/scaling, valmax=c*scaling, valinit=b)
        sliders[pos].on_changed(update_other_dim_values)
        
    def update_slider(old_dim, new_dim, scaling = 1.1):
        for i, sl in enumerate(sliders):
            tdim = sl.label.get_text()
            #print("old dim",old_dim, "new dim", new_dim, " names: ",tdim, names[new_dim])
            if tdim==names[new_dim]:
                #print("changing slider at pos",i)
                axsliders[i].clear()
                gen_slider(i, old_dim, scaling=scaling) 
                #a,b,c = value_range4dim(new_dim)
                #sliders[i].__init__(axsliders[i], names[new_dim], valmin=a/scaling, valmax=c*scaling, valinit=b)
        
    def update_displayed_dim0(valn):
        val = dim4name[valn]
        print("abc",valn, val)
        
        if dims[0] == int(val):
            print("Dim 0 unchanged (",dims[0],")", val)
            return
        elif dims[1] == int(val):
            print("Dim 0 cannot equal dim 1")
            radio0.set_active(dims[0])
            return                    
        
        update_slider(dims[0], int(val))
        dims[0] = int(val)
        #sca.set_offsets(np.array([X[::,dims[0]], X[::,dims[1]]]).T)        
        
        gi = np.where(y==0)[0]
        sca0.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
        gi = np.where(y>0)[0]
        sca1.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
        #plt.legend()
        
        fig.canvas.draw_idle()
        
    def update_displayed_dim1(valn):     
        val = dim4name[valn]
        #print(valn, val)
        if dims[1] == int(val):
            print("Dim 1 unchanged (",dims[1],")", val)
            return
        elif dims[0] == int(val):
            print("Dim 1 cannot equal dim 0")
            radio1.set_active(dims[1])
            return                    
        
        update_slider(dims[1], int(val))
        dims[1] = int(val)
        
        gi = np.where(y==0)[0]
        sca0.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
        gi = np.where(y>0)[0]
        sca1.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
        
        
        #sca.set_offsets(np.array([X[::,dims[0]], X[::,dims[1]]]).T)
        fig.canvas.draw_idle()

    
    sh = clf.shape_fit_
    if names is None:
        names = {i:"dim"+str(i) for i in range(sh[1])}
    dim4name = {name:i for i,name in enumerate(names.values())}    
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax.margins(x=0)
    
    #sca = ax.scatter(X[::,dims[0]], X[::,dims[1]], c=y, s=30, alpha=0.1, cmap=plt.cm.Paired)
    gi = np.where(y==0)[0]
    sca0 = ax.scatter(X[gi,dims[0]], X[gi,dims[1]], c=y[gi], s=30, alpha=0.7, vmin=0, vmax=1, cmap=plt.cm.Paired, label="Stars")
    gi = np.where(y>0)[0]
    sca1 = ax.scatter(X[gi,dims[0]], X[gi,dims[1]], c=y[gi], s=30, alpha=0.7, vmin=0, vmax=1, cmap=plt.cm.Paired, label="others")
    #plt.legend()

    axcolor = 'lightgoldenrodyellow'
    
    rax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)
    radio0 = RadioButtons(rax, list([str(names[i]) for i in range(sh[1])]), active=dims[0])
    radio0.on_clicked(update_displayed_dim0)

    ray = plt.axes([0.05, 0.5, 0.15, 0.15], facecolor=axcolor)
    radio1 = RadioButtons(ray, list([str(names[i]) for i in range(sh[1])]), active=dims[1])
    radio1.on_clicked(update_displayed_dim1)

    f0 = 1
    a0 = 1
    delta_f = 1
    ypos = 0.15
    sliders = []
    axsliders = []
    for i in range(sh[1]):
        if i==dims[0] or i==dims[1]: continue
        axslider = plt.axes([0.25, ypos, 0.65, 0.03], facecolor=axcolor)
        ypos-=0.05
        tmp = Slider(axslider, names[i], 0.1, 30.0, valinit=f0, valstep=delta_f)
        tmp.on_changed(update_other_dim_values)
        axsliders.append(axslider)
        sliders.append(tmp)
        gen_slider(len(axsliders)-1, i)
        
    draw_contours()

    
    #print(dir(sca))
    if sh[1]<2:
        raise Exception("multidim_visualization requires at least 2d training data (is "+str(sh)+").")
    #print(dir(clf))
    #print(np.shape(clf.support_vectors_))

    #sliders[1].label=Text(-0.02, 0.5, 'Arsch')
    #radio0.set_active(3)
    plt.show()
    
def get_props2(ifn):
    ff = pyfits.open(ifn)
    props = []
    
    gi = np.where(np.isfinite(ff[1].data["parallax"]) & (ff[1].data["parallax"]>0))[0]
    print("gi (parallax)",len(gi))
    Gmag = ff[1].data["Gmag"]
    Fg = 10**(-0.4* Gmag)*3.660e-08*720
    Fx = ff[1].data["ML_FLUX_0"]
    #props.append()
    props.append(ff[1].data["bp_rp"][gi])
    #
    print(np.shape(props))
    #

    props.append(np.log10(1000./ff[1].data["parallax"])[gi])

    tmp = ff[1].data["offset_sig"][gi].flatten()
    print("tmp",np.shape(tmp))
    props.append(tmp)
    print(np.shape(props))
    props.append(np.log10(Fx/Fg)[gi])    
    

    print(np.shape(props))
    
    #props.append(ff[1].data["ML_FLUX_0"])
 

    print(np.shape(props))
    props = np.array(props)
    #.T
    print(" x1",np.shape(props))
    sh = np.shape(props)
    gi = np.arange(sh[1])
    print("gi",np.shape(gi))
    for i in range(sh[0]):
        print("i ",i)
        finite = np.isfinite(props[i])
        print(len(finite), np.sum(finite))
        #fi = 
        #print(np.shape(fi), np.shape(gi))
        gi*=finite    
    print(" x2", np.shape(props))   
    
    props = props.T[gi, ::].T
    print(" x3", np.shape(props))   
    print()
    return props#[gi].T

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


if __name__ == "__main__":

    #props = ["bp_rp", "log_FxFg", "offset_sig", "log_distance"]
    #props = ["log_FG", "log_FX", "bp_rp", "offset_sig", "log_distance"]
    #X, y = get_props("../merged_training.fits", prop_cols=props)

    props = ["logFx","logFg","pos","log_plx","bp_rp"]
    props = ["pos", "logFxFg","bp_rp","log_plx"]
    
    #props = ["bp_rp", "logFg","logFx", "pos","log_plx"]
    
    X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")

    #clf = svm.SVC(class_weight={1: 3}, probability=True)
    clf = svm.SVC(C=5, kernel='poly', probability=True, degree=1,class_weight={0: 2})
    #clf = PCA(n_components=2)
    #clf = tree.DecisionTreeClassifier()
    #clf = svm.SVC(kernel='linear', probability=True,class_weight={1: 3})
    #clf = SGDClassifier(loss='hinge')
    clf.fit(X, y)
    b = clf.predict(X)
    recovery(y, b)
    print()
    #kernel = 1.0 * RBF(1.0)
    #gpc = GaussianProcessClassifier(kernel=kernel, random_state=0)
    #gpc.fit(X, y)
    #print(gpc.score(X, y))
    #c = gpc.predict(X)
    #recovery(y, c)
    
    #exit()
    multidim_visualization(clf, X, y, names={i:props[i] for i in range(len(props))})
    #exit()
    ##print(clf.support_vectors_)

    #ax = plt.gca()

    #plt.scatter(X[:, 0], X[:, 1], c=y, s=30, cmap=plt.cm.Paired)
    #xlim = ax.get_xlim()
    #ylim = ax.get_ylim()
    ### create grid to evaluate model
    #xx = np.linspace(xlim[0], xlim[1], 30)
    #yy = np.linspace(ylim[0], ylim[1], 30)

    #YY, XX = np.meshgrid(yy, xx)

    #tmp = [XX.ravel(), YY.ravel()]
    #for i in range(np.shape(X)[1]-2):
        #zz = np.median(X[i+2,::])
        #tmp.append(900*[zz])

    ##print(np.shape(X))
    ##print(XX, XX.ravel(), len(XX.ravel()))

    ##if np.shape(X)[1]==2:
        ##xy = np.vstack([XX.ravel(), YY.ravel()]).T
    ##else:
    #xy = np.vstack(tmp).T
    ##Z = clf.decision_function(xy).reshape(XX.shape)

    #Z = clf.predict_proba(xy)[::,0].reshape(XX.shape)
    #print(np.shape(Z))
    ### plot decision boundary and margins
    #ax.contour(XX, YY, Z, colors='k', levels=[-1, 0, 1], alpha=0.5,
            #linestyles=['--', '-', '--'])

    #ax.contour(XX, YY, Z, colors='b', levels=[0.1, 0.5, 0.9], alpha=0.5,
            #linestyles=['--', '-', '--'])

    ## plot support vectors
    ##ax.scatter(clf.support_vectors_[:, 0], clf.support_vectors_[:, 1], s=100,
            ##linewidth=1, facecolors='none', edgecolors='k')

    #b = clf.predict(X)
    #print(np.sum(b))
    #plt.show()

    print(80*"=")
    fn = "major_catalog.fits"
    ofn = "major_proba.fits"
    
    #fn = "test_catalog.fits"
    #fn = "stellar_test_sources.fits"
    #ofn = "stellar_test_proba.fits"
    fn = "random_matches.fits"
    ofn = "random_proba.fits"
    
    
    Y, eid, idx = get_props(fn, category_column=None, prop_cols=props, name_col="id_ero", with_index=True)
    print("Shape of property array", np.shape(Y))
    print("eid", eid)
    #print(np.sum(np.isfinite(Y[0])))
    #print(np.sum(np.isfinite(Y[1])))

    c = clf.predict(Y)
    
    multidim_visualization(clf, Y, c, names={i:props[i] for i in range(len(props))})
    
    pp0 = clf.predict_proba(Y)[::,0]
    pp1 = clf.predict_proba(Y)[::,1]
    #print(pp0+pp1)
    plt.scatter(c, pp1)
    plt.show()
    print("# recovered", np.sum(c), " from ",np.shape(pp1))
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
    exit()
    
