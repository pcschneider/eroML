from eroML.ensemble import fits_support
import numpy as np

@fits_support
def random_pos(e, random_category=2, category_col="category", abs_dist_col="match_dist", rel_dist_col="offset_sig"):
    abs_dist = e.to_array(abs_dist_col, array_type="array")
    rel_dist = e.to_array(rel_dist_col, array_type="array")
    cls = e.to_array(category_col, array_type="array")
    real_gi = np.where(cls != random_category)[0]
    rand_gi = np.where(cls == random_category)[0]
    print("Number of real sources: ", len(real_gi), " number of fake sources: ",len(rand_gi))
    return e
    #gi = np.where(e.
