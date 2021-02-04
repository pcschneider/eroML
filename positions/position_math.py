import numpy as np
from scipy.special import factorial


def analytic_probability(match_dist=None, sigma=None, sky_density=None, ps=0.1, N0=1):
    """
    Calculate the analytic probability for a correct identification based
    on match distance, positional uncertainty (sigma), and local 
    sky density.
    
    Parameters
    ----------
    match_dist : float (or array of float)
        Distance between catalog entries in arcsec
    sigma : float (or array of float)
        Estimated positional uncertainty
    sky_density : float (or array of float)
        Density of sources in #objects/arcmin^2
    ps : float
        Stellar fraction (Bayes factor)
    N0 : int
        Number of sources in match catalog
        
    Returns
    -------
    probabilities : float (or array of float)
    """
    
    
    s = sigma / (3600 * 180/np.pi)
    md = match_dist / (3600 * 180/np.pi)

    tmp1 = 2/s**2*np.exp(-md**2 / (2* s**2))                     
    tmp2 = ps/(1-ps) * 1/N0 * tmp1            

    p_stellar = tmp2 / (1+tmp2)
        
    return p_stellar

def gen_real_pos_offset(N=None, sigma=1.):
    """
    Generate random match distances based on the positional uncertainty
    
    Parameters
    ----------
    N : int
        Numbrt of sources
    sigma : float (or array of float)
        positional uncertainty, must be of length N if array
        
    Returns
    -------
    offs : array
        The simulated positional offset (shape: (N))
    """
    
    if N is None:
        if type(sigma) == float:
            raise ValueError("If 'sigma' is float, you must provide 'N'")
        sigma = np.array(sigma)       
        N = len(sigma)
    
    rnd = np.random.rand(N)
    rndx = np.sqrt(2*sigma**2 * (-np.log(1-rnd)))
    return rndx

def gen_random_pos_offset(N=None, dens=1., NN=3):
    """
    Generate random match distances for up to the NN-th neighbour
    
    Parameters
    -----------
    dens : array of float
        Density (arcsec-2) of stars
    NN : Simulate sources up the NN-th neighbour
      
    Returns
    -------
    simulated data : tuple of arrays
        offs, dens, group, nth
    """
    if N is None:
        if type(dens) == float:
            raise ValueError("If 'dens' is float, you must provide 'N'")
        dens = np.array(dens)
        N = len(dens)
    
    if type(dens) != float:
        dens = np.array(dens)

    M = NN
    NN*=7
    max_dist = np.repeat(np.sqrt(NN/np.pi/dens), NN*2).reshape((N,NN*2))
    #print(np.shape(max_dist), max_dist)
    dd = np.repeat(dens,NN*2).reshape((N,NN*2))
    gg = np.repeat(np.arange(len(dens)),NN*2).reshape((N,NN*2))
    #print(N, np.shape(max_dist), max(max_dist))
    NNsimu = np.random.poisson(lam=np.array(len(dens)*[NN]))
    #print(np.shape(NNsimu), NNsimu)
    offs = max_dist*np.random.rand(N,NN*2)**0.5
    nth = np.tile(np.arange(2*NN)+1, N)
    #print("offs: ",np.shape(offs))
    for j, n in enumerate(NNsimu):
        offs[j][n::] = np.nan
        #print(offs[j])
    offs = np.sort(offs, axis=1).flatten()
    
    step = 2*NN
    for i in range(step-M):
        j = M+i
        
        offs[j::step] = np.nan
        
    #plt.hist(offs, range=(0, 2.3), bins=30)
    ##print(offs[0,:])
    #plt.show()
    #print(np.shape(offs))
    #nth = np.argsort(offs, axis=1)
    #offs = offs[nth]
    #nth = np.tile(np.arange(NN),N).flatten()[np.argsort(offs).flatten()]+1
    #print(nth[0,:])
    #print(offs[0,nth[0,:]-1])
    #offs = offs.flatten()
    #gi = np.where(nth.flatten() == 1)
    #print(offs)
    #print(nth)
    #print(offs[nth])
    #print(np.shape(offs))
    oo = offs.flatten()
    gi = np.where(np.isfinite(oo))[0]
    #print(len(gi), len(dens)*NN)
    if len(gi) < N * M: 
        print("Simulated only ",len(gi), " instead of ",M*N, " sources; retrying...")
        return gen_random_pos_offset(N=N, dens=dens, NN=M)
    return np.array([offs.flatten()[gi], dd.flatten()[gi], gg.flatten()[gi], nth.flatten()[gi]])

