import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.text import Text
import pandas as pd

def multidim_visualization(clf, X, y, names=None, dims=[0,1]):    
    def draw_contours():
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ## create grid to evaluate model
        xx = np.linspace(xlim[0], xlim[1], 100)
        yy = np.linspace(ylim[0], ylim[1], 150)
        YY, XX = np.meshgrid(yy, xx)
        #print("lims: ",xlim, ylim)
        
        tmp = []
        for i in range(np.shape(X)[1]):
            if dims[0] == i:
                tmp.append(XX.ravel())
            elif dims[1] == i:
                tmp.append(YY.ravel())
            else:
                for j in range(np.shape(X)[1] - 2):
                    if sliders[j].label.get_text() == names[i]:
                        zz = sliders[j].val
                        tmp.append(15000*[zz])  
                
        xy = np.vstack(tmp).T
        #print(X.columns)
        try:
            xy = pd.DataFrame(data=xy, columns=X.columns)
        except:
            pass
        Z = clf.predict(xy).reshape(XX.shape)
        #Z0 = clf.decision_function(xy).reshape(XX.shape)
        try:
            Z1 = clf.predict_proba(xy)[::,1].reshape(XX.shape)
        except:
            Z1 = Z
            
        ## plot decision boundary and margins
        global con
        #con = ax.contour(XX, YY, Z0, colors='k', levels=[-0.5, 0, .5], alpha=0.5,
                #linestyles=['--', '-', '--', ':'])
        
        
        #gg0 = np.where(Z==0)
        #llim0 = min(Z1[gg0])
        #gg1 = np.where(Z==1)
        #llim1 = min(Z1[gg1])
        #print("lim:",llim0, llim1)
        #print(np.shape(Z1))
        #print(xy)
        
        #print()
        #print(Z)
        #plt.show()
        #plt.imshow(Z.T, origin='lower')
        #plt.show()
        ax.imshow(Z.T, origin='lower', extent=(min(xx), max(xx), min(yy), max(yy)), aspect='auto')
        #plt.show()
        #plt.scatter(Z1.flatten(), Z.flatten())
        #plt.show()
        
        con = ax.contour(XX, YY, Z1, colors='r', levels=[0.3, 0.5, 0.7], alpha=0.5,
                linestyles=['-.', '-', '--',':'])
        return con
     
    def update_other_dim_values(val):
        global con        
        if con is not None:
            for coll in con.collections: 
                ax.collections.remove(coll) 
        con = draw_contours()
            
        fig.canvas.draw_idle()
    
    def value_range4dim(dim):
        try:
            mn, mx = np.min(X.iloc[::,dim]), np.max(X.iloc[::,dim])
            mean, med = np.mean(X.iloc[::,dim]), np.median(X.iloc[::,dim])
        except:
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
        
        try:
            gi = np.where(y>0)[0]
            sca0.set_offsets(np.array([X.iloc[gi,dims[0]], X.iloc[gi,dims[1]]]).T) 
            gi = np.where(y==0)[0]
            sca1.set_offsets(np.array([X.iloc[gi,dims[0]], X.iloc[gi,dims[1]]]).T)
        except:
            gi = np.where(y>0)[0]
            sca0.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
            gi = np.where(y==0)[0]
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
        
        try:
            gi = np.where(y>0)[0]
            sca0.set_offsets(np.array([X.iloc[gi,dims[0]], X.iloc[gi,dims[1]]]).T) 
            gi = np.where(y==0)[0]
            sca1.set_offsets(np.array([X.iloc[gi,dims[0]], X.iloc[gi,dims[1]]]).T) 
        except:
            gi = np.where(y>0)[0]
            sca0.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
            gi = np.where(y==0)[0]
            sca1.set_offsets(np.array([X[gi,dims[0]], X[gi,dims[1]]]).T) 
        
        #sca.set_offsets(np.array([X[::,dims[0]], X[::,dims[1]]]).T)
        fig.canvas.draw_idle()

    
    N_props = len(names)
    N_data = np.shape(X)
    if N_data[1] < N_props: N_props = N_data[1]
    if names is None:
        names = {i:"dim"+str(i) for i in range(N_props)}
    dim4name = {name:i for i,name in enumerate(names.values())}    
    fig, ax = plt.subplots()
    #if N_props == 2:
    bb = 0.1 if N_props==2 else 0.35
    plt.subplots_adjust(left=0.25, bottom=bb, top=0.98)
    ax.margins(x=0)
    
    #sca = ax.scatter(X[::,dims[0]], X[::,dims[1]], c=y, s=30, alpha=0.1, cmap=plt.cm.Paired)
    try:
        gi = np.where(y>0)[0]
        sca0 = ax.scatter(X.iloc[gi,dims[0]], X.iloc[gi,dims[1]], c=y[gi], s=30, alpha=0.7, vmin=0, vmax=1, cmap=plt.cm.Paired, label="Stars")
        gi = np.where(y==0)[0]
        sca1 = ax.scatter(X.iloc[gi,dims[0]], X.iloc[gi,dims[1]], c=y[gi], s=30, alpha=0.7, vmin=0, vmax=1, cmap=plt.cm.Paired, label="others")
    except:
        gi = np.where(y>0)[0]
        sca0 = ax.scatter(X[gi,dims[0]], X[gi,dims[1]], c=y[gi], s=30, alpha=0.7, vmin=0, vmax=1, cmap=plt.cm.Paired, label="Stars")
        gi = np.where(y==0)[0]
        sca1 = ax.scatter(X[gi,dims[0]], X[gi,dims[1]], c=y[gi], s=30, alpha=0.7, vmin=0, vmax=1, cmap=plt.cm.Paired, label="others")        
    #plt.legend()

    axcolor = 'lightgoldenrodyellow'
    
    rax = plt.axes([0.05, 0.7, 0.15, 0.15], facecolor=axcolor)
    radio0 = RadioButtons(rax, list([str(names[i]) for i in range(N_props)]), active=dims[0])
    radio0.on_clicked(update_displayed_dim0)

    ray = plt.axes([0.05, 0.5, 0.15, 0.15], facecolor=axcolor)
    radio1 = RadioButtons(ray, list([str(names[i]) for i in range(N_props)]), active=dims[1])
    radio1.on_clicked(update_displayed_dim1)

    f0 = 1
    a0 = 1
    delta_f = 1
    ypos = 0.25
    sliders = []
    axsliders = []
    for i in range(N_props):
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
    if N_props<2:
        raise Exception("multidim_visualization requires at least 2d training data (is "+str(sh)+").")
    #print(dir(clf))
    #print(np.shape(clf.support_vectors_))

    #sliders[1].label=Text(-0.02, 0.5, 'Arsch')
    #radio0.set_active(3)
    ax.plot([1,10] ,[4.3,7], color='k', ls='--', lw=2)

    plt.show()
    
