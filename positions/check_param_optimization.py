from joblib import dump, load
import numpy as np
import matplotlib.pyplot as plt
from eroML.classify import my_custom_loss_func

grid = load("grid.joblib")

def one_panel(x, y, z, ax=None, xlabel="x", ylabel="y"):
    if ax == None:
        ax=plt.gca()
    x = np.array(x)
    y = np.array(y)
    idx = np.lexsort((x,y))    
    z = np.array(z)
    sx = np.sort(x)
    sy = np.sort(y)
    print(ax, x, y, z)
    Nx, Ny = len(np.unique(x)),len(np.unique(y))
    dx = sx[Ny]-sx[0]
    dy = sy[Nx]-sy[0]
    print(Nx, Ny, dx)
    print(dx, dy)
    im = z[idx].reshape((Ny,Nx))
    ax.imshow(im, extent=(min(x)-dx/2,max(x)+dx/2, min(y)-dy/2, max(y)+dy/2), aspect='auto', origin='lower')
    sc = ax.scatter(x, y, c=z)
    plt.colorbar(sc)
    
def stellar_recorvery_fraction(cv, X, y):
    
    gi = np.where(y==0)[0]
    N_stars = len(gi)
    for pp in cv.cv_results_["params"]:
        clf = cv.best_estimator_
        print(clf.clf__C)
        b = clf.predict(X)
    
    N_recovered = np.sum(b[gi])
    return N_recovered / N_stars
    
x = grid.cv_results_["param_clf__C"]
y = [a[0] for a in grid.cv_results_["param_clf__class_weight"]]
#z = grid.cv_results_["mean_test_score"]
z = stellar_recorvery_fraction(grid, grid.X_, grid.y_)
print(z)
one_panel(x,y,z)
plt.show()
exit()


print(grid.cv_results_["rank_test_score"])
for i in np.argsort(grid.cv_results_["rank_test_score"]):
    print(grid.cv_results_["rank_test_score"][i], "->",i)
    #i = np.where(grid.cv_results_["rank_test_score"]==j+1)[0][0]
    print( i, grid.cv_results_["mean_test_score"][i]/-1187639277885.3333)
    print(grid.cv_results_["param_clf__class_weight"][i], grid.cv_results_["param_clf__C"][i])
    print()


#for k in grid.cv_results_.keys():
    #print(k, grid.cv_results_[k])
