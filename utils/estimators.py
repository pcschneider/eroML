from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import PyAstronomy.funcFit as fuf


class Dist_model(fuf.OneDFit):
    """
    A model for the expected distance between real and random sources.
    """

    def __init__(self):
        fuf.OneDFit.__init__(self, ["fraction","sig","dens","N", "err_scaling"])
        self.Nsig = 1000
        self.sig = [1]
        self["sig"] = None
        
    def gen_sigs(self, sig_array):
        self.sig = np.random.choice(sig_array, self.Nsig)
        
    def evaluate(self, x):
        """
        Calculates and returns model according to the \
        current parameter values.

        Parameters:
        - `x` - Array specifying the positions at \
                which to evaluate the model.
        """
        
        # Real matches:
        s1 = self.sig*self["err_scaling"]
        s2 = np.repeat(np.array(s1)**2, len(x))
        #print(np.shape(s2))
        #print("s2", s2)
        #s2 = self["sig"]**2
        tmpx = np.tile(x, self.Nsig)
        tmp0 = tmpx/s2 * np.exp(-tmpx**2/(2*s2)) * self["fraction"]
        tmp0 = tmp0.reshape((self.Nsig, len(x)))
        #print(np.shape(tmp0))
        tmp0 = np.mean(tmp0, axis=0)
        #print(tmp0)
        
        ss2 = self["sig"]**2
        tmp00 = x/ss2 * np.exp(-x**2/(2*ss2)) * self["fraction"]
        #print("diff: ",np.array(tmp0 - tmp00)/tmp00)
        
        # Spurious (random) matches:
        
        Ng = self["dens"]
        tmp1 = 2*np.pi*x*Ng*np.exp(-np.pi*Ng*x**2) * (1-self["fraction"])
        
        return self["N"]*(tmp0+tmp1)



def NN_distribution(fn, Nsig=1000, ext=1, fkey="match_dist", errkey="RADEC_ERR", ofn=None, verbose=1):
    """
    Estimate number of real matches (as compared to random matches)
    
    The estimation is based on the 
    
    Example
    -------
    >>> from eroML.utils import NN_distribution 
    >>> N = NN_distribution("merged.fits", verbose=0)
    >>> abs(8663 - N)< 100
    True
    
    Parameters
    ----------
    fn : str
        Filename for the fits-file
    ext : int
        Extension in fits-file
    fkey : str
        Column name containing the match distance
    errkey : str
        Column name containing the positional error
    ofn : str (or None)
        Filename for diagnostic plots
        
    Returns
    -------
    estimated number of real matches : float
    
    """
    hist_range, hist_bins = (0, 80), 80
    dlim = 25 # assume that only random matches exist above this threshold
    disp = False if verbose<1 else True
        
    
    if verbose>0: print("Reading ",fn)
    ff = pyfits.open(fn)
    md = ff[ext].data[fkey] # match_distance
    N = len(md)
    
    hh = np.histogram(md, range=hist_range, bins=hist_bins)
    
    bin_x = 0.5*(hh[1][1::] + hh[1][0:-1])
    bin_w = np.mean(hh[1][1::] - hh[1][0:-1])
    bin_y = hh[0]
    if verbose>1: 
        print("Sources in hist: ",np.sum(hh[0]), " sources in fits-file: ",len(md))
        print("Bin width in hist: ",bin_w)
        
    mm = Dist_model()
    mm.Nsig = Nsig
    mm.gen_sigs(ff[ext].data[errkey])
    
    mm["N"] = N
    mm["sig"] = 2.5
    mm["dens"] = 4.5e-4
    mm["fraction"] = 0
    mm["err_scaling"] = 0.8
    mm.thaw(["dens","N"])
    #yerr = np.sqrt(bin_y)
    gi = np.where(bin_x > dlim)[0]

    mm.fit(bin_x[gi], bin_y[gi], maxfun=1e4, miniFunc="cash79", disp=disp)
    x = np.linspace(hist_range[0], hist_range[1],1000)

    if ofn is not None:
        fig = plt.figure(figsize=(8,11))
        ax0 = fig.add_subplot(311)
        ax0.hist(md, range=hist_range, bins=hist_bins)
        ax0.plot(x, mm.evaluate(x))
        
        ax = fig.add_subplot(312)
        #ax.plot(bin_x, bin_y, ls='steps-mid')
        ax.errorbar(bin_x, bin_y - mm.evaluate(bin_x), yerr=np.sqrt(bin_y))
        ax.errorbar(bin_x[gi], bin_y[gi] - mm.evaluate(bin_x[gi]), yerr=np.sqrt(bin_y[gi]), color='r')
        ax.axhline(y=0, color='k')
        ax.set_ylim(-100,100)
        
    #mm.parameterSummary()

    mm.freeze(["dens"])
    mm["fraction"] = 0.35
    mm["N"] = np.sum(hh[0])
    mm.thaw(["fraction", "err_scaling"])
    mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79", disp=disp)
    #mm.parameterSummary()

    mm.thaw(["dens"])
    mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79", disp=disp)
    
    if ofn is not None:
        ax = fig.add_subplot(313)
        #ax.plot(bin_x, bin_y, ls='steps-mid')
        ax0.plot(x, mm.evaluate(x))
        ax.errorbar(bin_x, bin_y - mm.evaluate(bin_x), yerr=np.sqrt(bin_y))
        #ax.errorbar(bin_x[gi], bin_y[gi] - mm.evaluate(bin_x[gi]), yerr=np.sqrt(bin_y[gi]), color='r')
        ax.axhline(y=0, color='k')
        ax.set_ylim(-100,100)
        plt.savefig(ofn)
        
    if verbose>0:
        mm.parameterSummary()
        print("Estimated number of sources: ",mm["N"], " instead of ", np.sum(hh[0]), " -> ",  np.sum(hh[0]) - mm["N"], "+-", np.sum(hh[0])**0.5,"(real-model)")
    
    return mm["N"] * mm["fraction"]


if __name__ == "__main__":
    
    import doctest
    #doctest.testmod()
    doctest.testmod()
    
  
