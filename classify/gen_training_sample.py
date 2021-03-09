from eroML.tile import file4
from eroML.ensemble import from_fits, to_fits
import numpy as np
import copy
from eroML.utils import activity_filter
from eroML.positions import random4dens, gen_real_pos_offset, add_random2real, calc_sigma_from_RADEC_ERR
import matplotlib.pyplot as plt

def randomize_match_distances(ensemble, category_col="category", sigma=None, max_dist=60):
    """
    Parameters
    ----------
    ensemble : Ensemble instance
    sigma : array of float with length of Ensemble
        Use these sigmas instead of the entries in ensemble
    max_dist : float
        Return sources expected within ``max_dist'' (unit arcsec)       
    """
    dens_scaling=1.03
    
    cat = ensemble.to_array(category_col, array_type="array")
    srcIDs = np.array(ensemble.srcIDs())
    reali = np.where(cat==0)[0]
    randi = np.where(cat==2)[0]
    rIDs = srcIDs[reali]
    Nreal = len(reali)
    RADEC_ERR = ensemble.to_array("RADEC_ERR", array_type="array")
    SIG = calc_sigma_from_RADEC_ERR(RADEC_ERR)
    sk = ensemble.to_array("eligible_sky_density", array_type="array")
    
    real_offs = gen_real_pos_offset(sigma=SIG[reali])
    real_n_rand = add_random2real(real=real_offs, skd=sk[reali]* dens_scaling, max_dist=max_dist)
    # offs, dens, group, nth, class
    
    rpt = np.unique(real_n_rand[2], return_counts=True)[1]
    #print(rpt, len(rpt))
    real = copy.deepcopy(ensemble)
    real.keep(rIDs)
    arr = real.to_array(ensemble.known_cols)
    arr = arr.repeat(rpt)
    
    n_srcID = [sID if NN==1 else sID+"_"+str("NN%i" % NN) for sID, NN in zip(np.repeat(rIDs, rpt), real_n_rand[3])]
    arr["srcID"] = n_srcID
    print("unique",len(np.unique(arr["srcID"])))
    #print(n_srcID)
    
    print(np.shape(arr), type(arr))
    print(arr["srcID"]) 
    real.from_array(arr, clean=True)
    
    print("rr",len(real), np.shape(arr))
    real.set_col("NN", real_n_rand[3])
    real.set_col("match_dist", real_n_rand[0])
    real.set_col("category", real_n_rand[4])
    print(Nreal, len(real_offs))
    
 
    sk = sk[randi]
    rand_offs = random4dens(dens=sk* dens_scaling)
    rpt = np.unique(rand_offs[2], return_counts=True)[1]
    # offs, dens, group, nth
    print(len(np.unique(rand_offs[2])), len(sk))
    randomIDs = srcIDs[randi]
    urIDs = randomIDs[np.unique(rand_offs[2])]
    print("urIDs",len(urIDs))
    rand = copy.deepcopy(ensemble)
    rand.keep(urIDs)
    arr = rand.to_array(ensemble.known_cols)
    arr = arr.repeat(rpt)
    
    n_srcID = [sID if NN==1 else sID+"_"+str("NN%i" % NN) for sID, NN in zip(np.repeat(urIDs, rpt), rand_offs[3])]
    arr["srcID"] = n_srcID
    print("unique2",len(np.unique(arr["srcID"])))    
    rand.from_array(arr, clean=True)
    rand.set_col("NN", rand_offs[3])
    rand.set_col("match_dist", rand_offs[0])
    rand.set_col("category", np.ones(len(n_srcID)).astype(int)*2)
    
    
    
    
    rand.append(real)
    
    return rand
        
        


def filter_training(ensemble):
    """
    KEEP only those objects in ensemble that have stellar like properties
    """
    dd = ensemble.to_array(["bp_rp", "FxFg", "phot_g_mean_mag","match_dist", "NN"], array_type="dict")
    af = activity_filter(dd["bp_rp"], dd["FxFg"])
    
    
    gi = np.where((af==0) & (dd["phot_g_mean_mag"] > 5.5) & (dd["match_dist"] < 10) & (dd["NN"] == 1))[0]
    srcIDs = np.array(ensemble.srcIDs())
    #print(srcIDs[gi])
    training.keep(srcIDs[gi])
    
def add_random(training, major, rnd_factor=1, max_dist=60, randomize=[]):
    """
    
    Parameters
    ----------
    training, major : Ensembles
    rnd_factor : float
        The ratio between real and random sources (within ``max_dist``)
    max_dist : float
        Maximum considered match_dist (in arcsec)
    randomize : list
        Add "sky_density" or "RADEC_ERR" to assign randomized values
        also for these two properties
    """
    Nrnd = int(rnd_factor * len(training))
    training.add_col("category", np.zeros(len(training)).astype(int))
    print("Adding %i random sources." % Nrnd)
    rsrcIDs = training.srcIDs()
    #msrcIDs = major.srcIDs()
    mdd = major.to_array(["srcID", "NN"], array_type="dict")
    msrcIDs = mdd["srcID"][mdd["NN"] == 1]
    ri = np.setdiff1d(msrcIDs, rsrcIDs)
    si = np.random.choice(ri, size=Nrnd, replace=False)
    print("unique", len(np.unique(si)))
    rnd = copy.copy(major)
    rnd.keep(si)
    #if "category" in rnd.known_cols:
        #rnd.set_col("category", np.ones(len(rnd))*2)
    #else:
    rnd.add_col("category", np.ones(len(rnd)).astype(int)*2)
    print(len(rnd), len(si))
    rnd.append(training)
    return rnd
    
    
    

conf_file = "eFEDS_EDR3.ini"
training_ID_file = "svm_training_IDs.txt"

mfn = file4("major", cconfig=conf_file)
rfn = file4("random", cconfig=conf_file)
print("Using ",mfn, " and ",rfn)

tIDs = np.genfromtxt(training_ID_file, dtype=str)
major = from_fits(mfn)

gi = np.in1d(major.srcIDs(), tIDs)
training = copy.deepcopy(major)
training.keep(tIDs)
print(len(training))
filter_training(training)
print(len(training))
final = add_random(training, major, rnd_factor=12.5)
final = randomize_match_distances(final)
to_fits(final, "train.fits", overwrite=True)




