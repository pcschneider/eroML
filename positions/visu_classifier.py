from eroML.classify import multidim_visualization
from joblib import dump, load

fn = "svc.joblib"
#fn = "svc_uniform.joblib"
#fn = "svc_linear.joblib"
#fn = "svc_LR.joblib"

print("Reading ",fn)

clf = load(fn)

multidim_visualization(clf, clf.X_, clf.y_, names={1:"offset", 2:"sky_dens", 0:"sig"})
