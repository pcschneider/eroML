from eroML.ensemble import fits_support
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass


@fits_support
def random_pos(e, random_category=2, category_col="category", abs_dist_col="match_dist", rel_dist_col="offset_sig"):
    abs_dist = e.to_array(abs_dist_col, array_type="array")
    rel_dist = e.to_array(rel_dist_col, array_type="array")
    sig = e.to_array("RADEC_ERR", array_type="array")
    skd = e.to_array("eligible_sky_density", array_type="array") /  3600.
    cls = e.to_array(category_col, array_type="array")
    real_gi = np.where(cls != random_category)[0]
    rand_gi = np.where(cls == random_category)[0]
    print("Number of real sources: ", len(real_gi), " number of fake sources: ",len(rand_gi))
    
    sigma = np.zeros(len(e))
    sigma[real_gi] = sig[real_gi]
    sigma[rand_gi] = np.sqrt(1/(2*np.pi*skd[rand_gi]))
    
    try:
        plt.hist(sigma, bins=30, range=(0, 30))
        plt.hist(sigma[real_gi], bins=30, range=(0, 30), label="real")
        plt.hist(sigma[rand_gi], bins=30, range=(0, 30), label="random")
        plt.legend()
        plt.show()
    except:
        pass
    
    rnd = np.random.rand(len(sigma))
    rndx = np.sqrt(2*sigma**2 * (-np.log(1-rnd)))

    try:
        plt.hist(rndx, bins=30, range=(0,30))
        plt.hist(rndx[real_gi], bins=30, range=(0,30), label="real")
        plt.hist(rndx[rand_gi], bins=30, range=(0,30), label="random")
        plt.legend()
        plt.show()
    except Exception as ee:
        print(ee)
        pass
    
    e.add_col("fake_match_dist", rndx)
    off_sig = rndx / sig
    e.add_col("fake_offset_sig", rndx)
    
    return e
    #gi = np.where(e.
