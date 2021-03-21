from eroML.tile import file4
from eroML.ensemble import from_fits, to_fits
import numpy as np
import copy
from eroML.utils import activity_filter
from eroML.utils.enrich import NN_Max
from eroML.positions import random4dens, gen_real_pos_offset, add_random2real, calc_sigma_from_RADEC_ERR
import matplotlib.pyplot as plt

def random_match_distances(ensemble, category_col="category", sigma=None, max_dist=60):
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
    
    print("Radnomizing match distances of %i sources, of which are %i real." % (len(ensemble), Nreal))
    
    RADEC_ERR = ensemble.to_array("RADEC_ERR", array_type="array")
    SIG = calc_sigma_from_RADEC_ERR(RADEC_ERR)
    sk = ensemble.to_array("eligible_sky_density", array_type="array")
    print("sigma: ",np.mean(SIG[reali]), np.median(SIG[reali]))
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
    NN_Max(rand)
    return rand
        
        
def random_props(train, random, category_col="category", keepNN1=True):
    """
    Parameters
    ----------
    train, random : Ensemble instances
    category_col : str
        Name of category column
    keepNN1 : boolean
        Keep properties of nearest neighbour association; otherwise randomize all associations
    """
    cat = train.to_array(category_col, array_type="array")
    NN = train.to_array("NN", array_type="array")
    srcIDs = np.array(train.srcIDs())
    realIDs = srcIDs[np.where(cat == 0)[0]]
    print("Real sources: ",len(realIDs))
    if keepNN1:
        gi = np.where((cat == 2) & (NN==1))[0]
        NN1_srcIDs = srcIDs[gi]
        gi = np.where((cat == 2) & (NN>1))[0]
    else:
        gi = np.where(cat == 2)[0]
        NN1_srcIDs = []

    # Identify associations to be randomized
    #    1) Occupied srcIDs (= real associations and (maybe) NN=1 sources)
    occupied_srcIDs = np.concatenate((realIDs, NN1_srcIDs))
    #    2) Select srcIDs 
    ii = np.in1d(srcIDs, occupied_srcIDs)
    srcIDs2rnd0 = srcIDs[~ii]
    indices2random = np.arange(len(srcIDs))[~ii]
    si = np.argsort(srcIDs[indices2random])
    indices2random = indices2random[si]
    
    #print("indices: ",indices2random, len(indices2random))
    #print("a",len(srcIDs2rnd0))
    #srcIDs2rnd1 = np.setdiff1d(srcIDs, occupied_srcIDs)
    #print("b",len(srcIDs2rnd1))
    #print("in", np.sum(np.in1d(srcIDs2rnd0, srcIDs2rnd1)))
    #print(" - ",np.sort(srcIDs2rnd0)[0:4], np.sort(srcIDs2rnd1)[0:4], np.sort(srcIDs[indices2random])[:4])
    #ii = np.in1d(srcIDs, occupied_srcIDs).nonzero()
    #print("ii", ii)
    #return
    #print(srcIDs[ii[0]], )
    osrcIDs2rnd = np.array([si[0:7] for si in srcIDs[indices2random]])
    print(osrcIDs2rnd[0:10])
    Nrnd = len(indices2random)
    print("Keeping properties for %i associations"  % len(occupied_srcIDs))
    print("Randomizing properties for %i associations. " % len(indices2random))
    
    mapper = np.zeros(len(indices2random)).astype(int)
    
    original_srcID = np.array([si[0:7] for si in random.to_array("original_srcID", array_type="array")])
    
    uoIDs, ioIDs, coIDs = np.unique(osrcIDs2rnd, return_inverse=True, return_counts=True)
    print("ioIDs: ",len(ioIDs))
    print(ioIDs)
    i0, i1 = 0, coIDs[0]
    for srnd, n in zip(uoIDs, coIDs):
        i1 = i0+n
        gi = np.where(srnd == original_srcID)[0]
        if n >= len(gi): rpl = True 
        else: rpl=False
        si = np.random.choice(gi, n, replace=rpl)
        mapper[i0:i1] = si
        
        #print(i0, i1, uoIDs[ioIDs[i0:i1]])
        #print("yyy",srcIDs[indices2random[i0:i1]], original_srcID[si])
        #print(srnd, n, original_srcID[gi])
        #print(si,"\n")
        i0+=n
    
    arr = train.to_array(train.known_cols)
    rarr= random.to_array(random.known_cols)
    
    for col in ["Fg", "phot_g_mean_mag", "parallax", "parallax_error", "bp_rp", "srcID_NN", "FxFg"]:
        print(col, indices2random, mapper)
        arr[col][indices2random] = rarr[col][mapper]
    #arr["FxFg"][indices2random][arr["FxFg"][indices2random] < 1e-5] = 0.1
    
        
    arr = train.from_array(arr, clean=True)    
    return train
        


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
        also for these two properties for the objects in the training sample
    """
    Nreal = len(training)
    Nmajor = len(major)
    Nrnd = int(rnd_factor * Nreal)
    training.add_col("category", np.zeros(len(training)).astype(int))
    print("Adding %i random sources." % Nrnd)
    
    if "RADEC_ERR" in randomize:
        radec_err = np.random.choice(major.to_array("RADEC_ERR", array_type="array"), Nreal)
        training.set_col("RADEC_ERR", radec_err)
        print("Randomizing RADEC_ERR for training objects.")
        
    if "sky_density" in randomize:
        sk = np.random.choice(major.to_array("eligible_sky_density", array_type="array"), Nreal)
        training.set_col("eligible_sky_density", sk)
        print("Randomizing sky_density for training objects.")
    
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
    
    
    

conf_file = "eFEDS_EDR3_HamStar.ini"
training_ID_file = "svm_training_IDs_HamStar.txt"

mfn = file4("major", cconfig=conf_file)
rfn = file4("random", cconfig=conf_file)
print("Using ",mfn, " and ",rfn)

tIDs = np.genfromtxt(training_ID_file, dtype=str)
major  = from_fits(mfn)
random = from_fits(rfn)

gi = np.in1d(major.srcIDs(), tIDs)
training = copy.deepcopy(major)
training.keep(tIDs)
sig = training.to_array("RADEC_ERR", array_type="array")
print("sigma0: ",np.mean(sig), np.median(sig), len(sig))

sig = major.to_array("RADEC_ERR", array_type="array")
print("sigma1: ",np.mean(sig), np.median(sig), len(sig))


print(len(training))
filter_training(training)
print(len(training))
#final = add_random(training, major, rnd_factor=12.5, randomize=["RADEC_ERR","sky_density"])
final = add_random(training, major, rnd_factor=35, randomize=["RADEC_ERR","sky_density"])
final = random_match_distances(final)
final = random_props(final, random)
to_fits(final, "train_HamStar.fits", overwrite=True)




