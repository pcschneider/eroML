from eroML.ensemble import from_fits
from eroML.tile import file4
from eroML.utils import activity_filter
import numpy as np
from eroML.ensemble import Ensemble

def training_filter(e, max_dist=1500, max_Lx=2e31 , verbose=1):
    """
    Check for distance, max Lx, and color-dependent activity filter
    
    Parameters
    ----------
    e : str, Ensemble, or dictionary
    max_dist : float
        Maximum allowed distance in pc
    max_Lx : float
        Maximum Lx in erg/s
        
    Returns
    -------
    good : array
        0 if bona fide star, else 1
    """
    if type(e) == type("xx"):
        if verbose>0:
            print("Reading ",e)
        e = from_fits(e)
        
    if type(e) == Ensemble:
        if verbose>0: print("Populating dictionary from Ensemble")
        dct = e.to_array(["parallax", "bp_rp", "FxFg", "Fx"])
    else:
        dct = x
    
    below = (1-activity_filter(dct["bp_rp"], dct["FxFg"])).astype(bool)
    cl = np.ones(len(dct["Fx"]))
    cl[below] = 0
    dist = 1000/dct["parallax"]
    Fx = dct["Fx"]
    Lx = 4*np.pi*(3.18e-2*dist)**2 * Fx * 1e40
    gi = np.where(((Lx > max_Lx) & (cl==0)) | (dist>max_dist))[0]
    cl[gi] = 1
    if verbose>0:
        gi = np.where(cl == 0)[0]
        print("training_filter - %i from %i are eligible training objects." % (len(gi), len(Fx)))
    return cl

if __name__ == "__main__":
    mfn = file4("major", cconfig="eFEDS_EDR3_HamStar.ini")
    #e = from_fits(mfn)
    cl = training_filter(mfn)
                      
