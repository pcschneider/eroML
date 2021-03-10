import numpy as np
from scipy.special import factorial
from scipy.stats import uniform
from math import ceil
from scipy.stats import poisson


def generate_random_densities(N0, N1, dens0=0.1, dens1=1, dens_scaling='uniform'):
    """
    Generate random densities
    
    Parameters
    ----------
    N0, N1 : int
        Number of real and random sources
    dens0, dens1 : float
        Min and max sky densities (in units of TBD)
    dens_scaling : str
        Must be within ['uniform', 'proportional', 'uni_prop']
        
    Returns
    -------
    dens_real : array 
             Densities for real sources
    dens_rand : array 
             Densities for random sources        
    """
    
    print("Simulating ",N0+3*N1, " sources, with a real fraction of ",N0/N1)


    lc = (dens0/dens1)**2

    if dens_scaling == 'uniform':
        dens = uniform.rvs(size=N0+N1, loc=dens0, scale=(dens1-dens0)) # Uniform scaling in density
        dens_real = dens[0:N0]
        dens_rand = dens[N0::]
    elif dens_scaling == 'proportional':
        dens = uniform.rvs(size=N0+N1, loc=lc, scale=1-lc)**0.5 * dens1 # Scaling proportional to density    
        dens_real = dens[0:N0]
        dens_rand = dens[N0::]
    elif dens_scaling == 'uni_prop':
        dens_real = uniform.rvs(size=N0, loc=dens0, scale=(dens1-dens0))
        dens_rand = uniform.rvs(size=N1, loc=lc, scale=1-lc)**0.5 * dens1        
        
    return dens_real, dens_rand
        


def analytic_probability(match_dist=None, sigma=None, sky_density=None, ps=0.1):
    """
    Calculate the analytic probability for a correct identification based
    on match distance, positional uncertainty (sigma), and local 
    sky density. Considers only nearest neighbour.
    
    Parameters
    ----------
    match_dist : float (or array of float)
        Distance between catalog entries in arcsec
    sigma : float (or array of float)
        Estimated positional uncertainty (in arcsec)
    sky_density : float (or array of float)
        Density of sources in #objects/arcmin^2
    ps : float
        Stellar fraction (Bayes factor)
        
    Returns
    -------
    probabilities : float (or array of float)
    
    """
    N0 = 41253 * 3600 * sky_density
    # sky area: 41,253 deg^2
    
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


def add_random2real(real=None, skd=1, max_dist=60, verbose=1):
    """
    Create a tuple with real and random sources. Most importantly, calculate
    nearest neighbour (NN).
    
    Parameters
    ----------
    real : array of float 
        Offsets for presumably real sources
    skd : array of float (same length as 'real')
        Sky density (in units of stars/arcmin^2)
    max_dist : float
        Maximum match distance (in arcsec)
        
    Returns
    -------
    simulated data : tuple of arrays
        offs, dens, group, nth, class
    """
    N = len(real)
    print("Generating random sources for %i real sources." % N)
    offs, dens, group, nth = random4dens(skd, max_dist=max_dist) 
    Nrnd = len(offs)
    
    out_offs = np.zeros(N+Nrnd)
    out_dens = np.zeros(N+Nrnd)
    out_group = np.zeros(N+Nrnd)
    out_nth = np.zeros(N+Nrnd)
    out_cls = np.ones(N+Nrnd)
    
    i0, i1 = 0, None
    for i in range(N):
        left = np.searchsorted(group, i)
        right = np.searchsorted(group, i, side='right')
        Ngrp = right-left
        i1 = i0+Ngrp+1
        all_dists = np.concatenate(([real[i]], offs[left:right]))
        
        
        out_offs[i0:i1] = all_dists
        out_dens[i0:i1] = skd[i]
        out_group[i0:i1] = i
        out_nth[i0:i1] = np.argsort(np.argsort(all_dists))+1
        out_cls[i0] = 0
        out_cls[(i0+1):i1] = 2
        
        if verbose>1:
            print(left, right, "i",i, "NN max: ",Ngrp, "(i0, i1: ",i0,i1,')')
            if Ngrp>0:
                print(" grp: ",group[left], group[right-1], " - ", group[left:right], nth[left:right])
                print("   offs: ",offs[left:right], "real: ",real[i])
            print(out_offs[i0:i1], out_dens[i0:i1], out_group[i0:i1], out_nth[i0:i1], out_cls[i0:i1])
            print()
            
        
        i0 = i1
                   
    print("nth: ",out_nth)
    return out_offs, out_dens, out_group, out_nth, out_cls

def Nexp4dens(dens, max_dist=60):
    """
    Calculate number of expected sources up to ``max_dist``.
    
    Parameters
    ----------
    dens : array 
          Unit:  #source/arcmin^2
    max_dist : float
          Maximum distances (in arcsec)
    
    Returns
    -------
    N : float
        Number of expected random sources
    """
    Ni = dens/3600*max_dist**2*np.pi
    return np.sum(Ni)

def calc_sigma_from_RADEC_ERR(RADEC_ERR):
    """
    Calculate the sigma for the Rayleigh distribution
    """
    tmp = np.sqrt(((RADEC_ERR/1.15)**2 + 0.7**2)/2)
    #tmp = RADEC_ERR
    tmp = RADEC_ERR * 0.615 # = /sqrt(2)/1.15
    return tmp

    

def random4dens(dens, max_dist=60):
    """
    Generate random sources for the provided sky densities.
    
    The number of random sources will be proportional to the individual sky 
    densities, the expectation value being 
    
    .. math::
    
        N_i^{exp} = \pi*dens_i*max\_dist^2\,.
        
    Parameters
    ----------
    dens : array 
          Unit:  #source/arcmin^2
    max_dist : float
          Maximum distances (in arcsec)
    
    Returns
    -------
    simulated data : tuple of arrays
        offs, dens, group, nth

    """
        
    N = len(dens)    
    N_exp_max = ceil(max(dens)/3600*max_dist**2*np.pi)
    print("Maximum number of sources in match radius: %6.2f." % N_exp_max)
    #offs = max_dist*np.random.rand(N,NN*2)**0.5
    offs = max_dist*np.random.rand(N*2*N_exp_max)**0.5

    Nexp = dens/3600*max_dist**2*np.pi # /3600 to convert from arcmin^-2 to arcsec^-2 
    Nmeaus = poisson.rvs(Nexp, size=N)
    grp = np.arange(N)

    skdens = np.repeat(dens, Nmeaus)
    offs_idx = np.repeat(grp, Nmeaus)
    pos_off = offs[0:len(offs_idx)]
    sigout = [np.nan]*len(offs_idx)
    nth = np.zeros(len(offs_idx))
    j = 0
    print("calculating 'nth'")
    for i, g in enumerate(grp):
        gi = np.where(offs_idx == g)[0]
        n = np.argsort(np.argsort(pos_off[gi]))+1
        nth[j:j+len(gi)] = n
        #print(g, len(gi), Nmeaus[i], pos_off[gi], n)
        j+=len(gi)
    print("                done.")    
    print("Returning ",len(pos_off), " random sources.")    
    return pos_off, skdens, offs_idx, nth    


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

    oo = offs.flatten()
    gi = np.where(np.isfinite(oo))[0]
    #print(len(gi), len(dens)*NN)
    if len(gi) < N * M: 
        print("Simulated only ",len(gi), " instead of ",M*N, " sources; retrying...")
        return gen_random_pos_offset(N=N, dens=dens, NN=M)
    return np.array([offs.flatten()[gi], dd.flatten()[gi], gg.flatten()[gi], nth.flatten()[gi]])

