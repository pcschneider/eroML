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
            
        Z = clf.decision_function(xy).reshape(XX.shape)
        #Z = clf.predict_proba(xy)[::,0].reshape(XX.shape)
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
    
    
dd = np.genfromtxt("offs.dat", unpack=True)
#dd[1]*=100
X = np.transpose(dd[0:3])
y = dd[3]
y[y>0] = 1
#X, y = get_props("../merged_training.fits", prop_cols=props,category_column="category")

#clf = svm.SVC(class_weight={1: 3}, probability=True)
clf = svm.SVC(C=10, kernel='poly', degree=3,class_weight={0: 2}, gamma=0.01)
#clf = PCA(n_components=2)
#clf = tree.DecisionTreeClassifier()
#clf = svm.SVC(kernel='linear', probability=True,class_weight={1: 3})
#clf = SGDClassifier(loss='hinge')
clf.fit(X, y)
b = clf.predict(X)
#recovery(y, b)
print()
#kernel = 1.0 * RBF(1.0)
#gpc = GaussianProcessClassifier(kernel=kernel, random_state=0)
#gpc.fit(X, y)
#print(gpc.score(X, y))
#c = gpc.predict(X)
#recovery(y, c)
gi = np.where(b == y)[0]
print("Correct: ",len(gi), " of ",len(b))

gi = np.where((b == 0) & (y==0))[0]
print("star as star: ",len(gi), " of ",len(np.where(y==0)[0]))

gi = np.where((b == 1) & (y==1))[0]
print("random as random: ",len(gi), " of ",len(np.where(y==1)[0]))

gi = np.where((b == 0) & (y==1))[0]
print("random as star: ",len(gi), " of ",len(b))
gi = np.where((b == 1) & (y==0))[0]
print("star as random: ",len(gi), " of ",len(b))

#exit()
multidim_visualization(clf, X, y, names={1:"offset", 2:"sky_dens", 0:"sig"})
#exit()    
