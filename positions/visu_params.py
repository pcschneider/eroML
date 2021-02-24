import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter

def make_image(x, y, z):    
    x = np.array(x)
    y = np.array(y)
    idx = np.lexsort((x,y))    
    z = np.array(z)
    sx = np.sort(x)
    sy = np.sort(y)
    #print(x, y, z)
    Nx, Ny = len(np.unique(x)),len(np.unique(y))
    dx = sx[Ny]-sx[0]
    dy = sy[Nx]-sy[0]
    print(Nx, Ny)
    print(dx, dy)
    im = z[idx].reshape((Ny,Nx))
    return im

def four_panels(dd, xindex=0, yindex=1, xlabel="C", ylabel="weight", zindices=[2,3,4,5],\
                ztitles=["Stars as stellar", "Others as others", "Stars as others",\
                         "Others as stars"]):
    """
    """
    print("First line: ",dd[xindex])
    x = np.sort(np.unique(dd[::,0]))
    y = np.unique(dd[::,yindex])
    print("x: ",x)
    print("y: ",y)
    xlog = np.sort(np.log10(np.unique(dd[::,xindex])))
    dlx = xlog[1]-xlog[0]
    dly = np.log10(y[2]) - np.log10(y[1])
    x0, x1 = min(x), max(x)
    y0, y1 = min(y), max(y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    print("d log x: ",dlx, " dx: ",dx)
    print("d log y: ",dly, " dy: ",dy)

    Nstars = dd[0,6]
    print("Nstars: ",Nstars)

    l4max = np.max(dd[::,4])
    l5max = np.max(dd[::,5])

    fig = plt.figure()
    

    for i, zi in enumerate(zindices):
        ax = fig.add_subplot(2,2,i+1)
        im = make_image(dd[::,xindex], dd[::,yindex], dd[::, zi])
        if i in [0,2,3]:
            axim = ax.imshow(im, aspect='auto', origin='lower',extent=(x0-dx/2, x1+dx/2, y0-dy/2, y1+dy/2))
        else:    
            axim = ax.imshow(im, aspect='auto', origin='lower',extent=(x0-dx/2, x1+dx/2, y0-dy/2, y1+dy/2))
        ax.set_title(ztitles[i]+" ("+str(i)+")")
        cb = fig.colorbar(axim, ax=ax)
        if i==1 or i==3:
            ax.set_yticklabels([])
        if i==0 or i==1:
            ax.set_xticklabels([])
        if i==2 or i==3:
            ax.set_xlabel(xlabel)
        if i==0 or i==2:
            ax.set_ylabel(ylabel)    
    plt.show()
    
Blues = plt.get_cmap('Blues')
Reds = plt.get_cmap('Reds')



dd = np.genfromtxt("params_scale.txt", unpack=False)
dd = np.genfromtxt("params2.txt", unpack=False)
    #oo.write("#  0 - C\n")
    #oo.write("#  1 - class weight for class 0 \n")
    #oo.write("#  2 - i0: Stars as stellar recovered\n")
    #oo.write("#  3 - i1: Others as others recovered\n")
    #oo.write("#  4 - i2: Stars as others recovered\n")
    #oo.write("#  5 - i3: Others as stars recovered\n")
    #oo.write("#  6 - Nstars: Stars in sample\n")
    #oo.write("#  7 - Ntotal: Sample size\n")


#four_panels(dd)
 

print("First line: ",dd[0])
x = np.sort(np.unique(dd[::,0]))
y = np.unique(dd[::,1])
print("x: ",x)
print("y: ",y)
xlog = np.sort(np.log10(np.unique(dd[::,0])))
dlx = xlog[1]-xlog[0]
dly = np.log10(y[2]) - np.log10(y[1])
x0, x1 = min(x), max(x)
y0, y1 = min(y), max(y)
dx = x[1] - x[0]
dy = y[1] - y[0]

print("d log x: ",dlx, " dx: ",dx)
print("d log y: ",dly, " dy: ",dy)

Nstars = dd[0,6]
print("Nstars: ",Nstars)

l4max = np.max(dd[::,4])
l5max = np.max(dd[::,5])

classified_as_star = dd[::,2] + dd[::,5]  - Nstars
imbalance = abs(dd[::,4] - dd[::,5])
wrong = dd[::,4] + dd[::,5]
wrong_im  = make_image(dd[::,0], dd[::,1], wrong)
imbalance_im = make_image(dd[::,0], dd[::,1], imbalance)
classified_as_star_im = make_image(dd[::,0], dd[::,1], classified_as_star)

xn = np.arange(len(x))
fig = plt.figure(figsize=(14,7))

ax0 = fig.add_subplot(1,2,1)
im = ax0.imshow(wrong_im, origin='lower', aspect='auto', extent=(-0.5, len(x)-0.5, y0-dy/2, y1+dy/2))
X, Y = np.meshgrid(xn, y)
#ax.contour(X, Y, gaussian_filter(imbalance_im, 0.2), levels=(100,200, 300), colors=['w','r','b'])
mwi = np.min(wrong_im)
CS = ax0.contour(X, Y, wrong_im, levels=(mwi+10, mwi+20, mwi+30), colors='k')
CS2 = ax0.contour(X, Y, imbalance_im, levels=(50, 100, 200), colors='r')
plt.colorbar(im, ax=ax0)
ax0.clabel(CS, inline=1, fontsize=10)
ax0.clabel(CS2, inline=1, fontsize=10)

#plt.show()
#fig = plt.figure(figsize=(14,7))

ax1 = fig.add_subplot(1,2,2)
im = ax1.imshow(imbalance_im, origin='lower', aspect='auto', extent=(-0.5, len(x)-  0.5, y0-dy/2, y1+dy/2))
#X, Y = np.meshgrid(x, y)
#print(x, y)
#ax.contour(X, Y, gaussian_filter(imbalance_im, 0.2), levels=(100,200, 300), colors=['w','r','b'])
CS = ax1.contour(X, Y, imbalance_im, levels=(50, 100, 200), colors='k')
CS2 = ax1.contour(X, Y, wrong_im, levels=(mwi+10, mwi+20, mwi+30), colors='r')
plt.colorbar(im, ax=ax1)
ax1.clabel(CS, inline=1, fontsize=10)
ax1.clabel(CS2, inline=1, fontsize=10)
ax1.set_xlabel("C")
ax0.set_xlabel("C")
ax0.set_ylabel("Class Weight (for class=0)")

for ax in [ax0, ax1]:
    xt = np.array(ax.get_xticks())
    gi = np.where((xt>=0) & (xt<len(x)))[0]
    print(xt, xt[gi], len(x))
    print(x[xt[gi].astype(int)])
    #ax.set_xticks(xt[gi])
    labels=[str(x[int(i)]) for i in xt[gi]] 
    print(labels)   
    #plt.xticks(ticks=xt[gi], labels=labels)
    ax.set_xticklabels(labels)
plt.show()

exit()
for line in dd:
    #print(line)
    xx = line[0]
    yy = line[1]
    #print(xx,yy)
    c1 = Blues(line[4]/l4max) # stars as other
    c2 = Reds(line[5]/l5max) # others as stars
    
    t1 = plt.Polygon(np.transpose([[xx-dx/2,xx-dx/2, xx+dx/2], [yy-dy/2,yy+dy/2,yy+dy/2]]), color=c1)
    plt.gca().add_patch(t1)
    
    t2 = plt.Polygon(np.transpose([[xx-dx/2,xx+dx/2, xx+dx/2], [yy-dy/2,yy-dy/2,yy+dy/2]]), color=c2)
    plt.annotate("%i" % line[4], xy=(xx-dx/4,yy+dy/4))
    plt.annotate("%i" % line[5], xy=(xx,yy-dy/4))
    plt.gca().add_patch(t2)

plt.xlim(x0, x1)
plt.ylim(y0, y1)
plt.xlabel("C")
plt.ylabel("class=0 weight")
plt.show()
