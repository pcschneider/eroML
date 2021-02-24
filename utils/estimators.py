from astropy.io import fits as pyfits
import numpy as np
try:
    import matplotlib.pyplot as plt
except:
    pass    
import PyAstronomy.funcFit as fuf
import astropy.units as u


def sky_dens4coordinates(coord, around=4,auto_adjust=True,  verbose=1):
    """
    Calculate the sky density using the provided coordinates
    
    Parameters
    ----------
    coord : SkyCoord instance
    around : float
        The default search area in arcmin (the search radius that is used may be smaller in high-density regions

    Returns
    -------
    stars per arcmin^2 : np.array
    """
    
    def split(cc, N):
        r = []
        step = len(coord)//N
        si = np.random.permutation(range(len(cc)))
        i0, i1 = 0, step
        for i in range(N):
            if i==N-1:
                i1 = len(coord)
            r.append(si[i0:i1])
            i0=i1
            i1+=step
        return r    
    
    if verbose>1:
        print("sky_dens:: Searching around ",around, "arcmin for ",len(coord)," objects.")
        print("sky_dens:: Coord range: (",min(coord.ra.degree), max(coord.ra.degree), " ; ", min(coord.dec.degree), max(coord.dec.degree),")")

    ra_range0 = abs(max(coord.ra.degree) - min(coord.ra.degree))
    ra_range1 = abs(max((coord.ra.degree  + 180) % 360)- min((coord.ra.degree + 180) % 360))
    #print("ra_range0, ra_range1",ra_range0, ra_range1)
    ra_range = ra_range0 if ra_range0<ra_range1 else ra_range1
    dec_range = abs(max(coord.dec.degree) - min(coord.dec.degree))
    sky_area = ra_range *dec_range * np.cos(np.nanmean(coord.dec.degree)/180*np.pi)
    #print("cos", np.cos(np.nanmedian(coord.dec.degree/180*np.pi)))
    sky_dens = len(coord)/sky_area/3600 # per arcmin^2
    if verbose>0: 
        print("sky_dens:: Sky area of Ensemble: ",sky_area, " (center: ",np.nanmedian(coord.ra.degree), np.nanmedian(coord.dec.degree),")")         
        print("sky_dens::    ",ra_range, dec_range)
        print("sky_dens::    Mean sky density: ",sky_dens," #stars/arcmin^2")
        print("sky_dens::    Mean number of stars in search radius: ",sky_dens*np.pi*around**2)
    
    N = len(coord)
    ret = np.zeros(N)
    ret[::] = np.nan
    
    if sky_dens > 3:
        #return None
        print("sky_dens:: splitting Ensemble...")
        cc = split(coord, 3)
        for c in cc:
            #print("ccc",c, len(c), min(c), max(c))
            ret[c] = sky_dens4coordinates(coord[c], around=around) * 3
        #return ret
    
    elif (sky_dens*np.pi*around**2 < 10) and auto_adjust:
        ta = (10/sky_dens/np.pi)**0.5
        print("sky_dens:: Adjusting 'search_around' to ",ta)
        return sky_dens4coordinates(coord, around=ta, auto_adjust=False)
        #idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, ta)
        #uni, cnt = np.unique(idxc, return_counts=True)
        ##print(uni)
        ##print(cnt, min(cnt), max(cnt), around, np.mean(cnt), len(cnt))
        #ret=np.array(cnt)/np.pi*(around/0.5)**2
        #print("sky_dens:: Using 0.5 arcmin search radius and extrapolating.", np.nanmean(ret), len(ret))        
    ##elif sky_dens > 3:
        ##ta = 0.7*u.arcmin
        ##idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, ta)
        ##uni, cnt = np.unique(idxc, return_counts=True)
        ##cnt=np.array(cnt)*(around/0.7)**2
        ##print("sky_dens:: Using 0.7 arcmin search radius and extrapolating.")
    else:
        idxc, idxcatalog, d2d, d3d = coord.search_around_sky(coord, around*u.arcmin)
        uni, cnt = np.unique(idxc, return_counts=True)
        #print(cnt)
        ret = (cnt-1)/np.pi/around**2
        #ret = cnt/np.pi/around**2 / 2.
        print("sky_dens:: Using ",around," arcmin search radius:", np.nanmean(ret))
        print("sky_dens::        median number of stars in search radius: ",np.median(cnt-1).astype(int))
        print("sky_dens::        mean number of stars in search radius: ",np.mean(cnt-1).astype(int))
        print("sky_dens::        number of unique densities: ",len(np.unique(ret)))
        
    #print(ret)
    print("sky_dens:: nanmean: ",np.nanmean(ret))
    return ret

class Dist_model(fuf.OneDFit):
    """
    A model for the expected distance between real and random sources.
    """

    def __init__(self):
        fuf.OneDFit.__init__(self, ["fraction","sig","dens","N", "err_scaling"])
        self.Nsig = 1000
        self.sig = [1]
        self["sig"] = 1.
        self["err_scaling"] = 1.
        self.gen_sigs([self["sig"]])
        
    def gen_sigs(self, sig_array):
        self.sig = np.random.choice(sig_array, self.Nsig)
        #mn = np.nanmedian(self.sig)
        ##st = np.std(self.sig)
        #self.sig-=mn
        #self.sig*=0.2
        #self.sig+=mn
        #plt.hist(np.abs(self.sig), range=(0,10), bins=30)
        #print("mn", mn)
        #plt.show()
        
    def evaluate(self, x):
        """
        Calculates and returns model according to the \
        current parameter values.

        Parameters
        ----------
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
        
    # Read file    
    if type(fn) == pyfits.HDUList:
        if verbose>0:
            print("utils.estimators.NN_distribution - fn is HDUList")
        ff = fn
    else:
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
    mm["sig"] = 3.7496337890625 * 0.75
    mm["dens"] = 4.5e-4
    mm["fraction"] = 0
    mm["err_scaling"] = 0.8
    mm.thaw(["dens","N"])
    #yerr = np.sqrt(bin_y)
    gi = np.where(bin_x > dlim)[0]

    mm.fit(bin_x[gi], bin_y[gi], maxfun=1e4, miniFunc="cash79", disp=disp)
    mm.parameterSummary()
    x = np.linspace(hist_range[0], hist_range[1],1000)

    if ofn is not None:
        fig = plt.figure(figsize=(8,11))
        ax0 = fig.add_subplot(311)
        ax0.hist(md, range=hist_range, bins=hist_bins)
        ax0.plot(x, mm.evaluate(x), ls='--')
        
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
    mm["err_scaling"]=0.75
    #mm.thaw(["fraction", "err_scaling"])
    mm.thaw(["fraction"])
    #mm.thaw(["fraction","sig"])
    mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79", disp=disp)
    #mm.parameterSummary()

    mm.thaw(["dens"])
    mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79", disp=disp)
    #mm["dens"] = 0.000441403 * 10
            
    if verbose>0:
        mm.parameterSummary()
        print("Estimated number of sources: ",mm["N"], " instead of ", np.sum(hh[0]), " -> ",  np.sum(hh[0]) - mm["N"], "+-", np.sum(hh[0])**0.5,"(real-model)")
    
    ff.close()
    if ofn is not None:
        ax = fig.add_subplot(313)
        #ax.plot(bin_x, bin_y, ls='steps-mid')
        ax0.plot(x, mm.evaluate(x), lw=2, color='r')
        ff, NN = mm["fraction"], mm["N"]
        mm["fraction"] = 1.
        mm["N"] = ff*NN
        ax0.plot(x, mm.evaluate(x))
        mm["fraction"] = 0
        mm["N"] = NN * (1-ff)
        
        ax0.plot(x, mm.evaluate(x), ls='-')

        mm["fraction"] = ff
        mm["N"] = NN
        ax.errorbar(bin_x, bin_y - mm.evaluate(bin_x), yerr=np.sqrt(bin_y))
        #ax.errorbar(bin_x[gi], bin_y[gi] - mm.evaluate(bin_x[gi]), yerr=np.sqrt(bin_y[gi]), color='r')
        ax.axhline(y=0, color='k')
        ax.set_ylim(-100,100)
        plt.savefig(ofn)
        plt.show()

    return mm["N"] * mm["fraction"]


if __name__ == "__main__":
    
    import doctest
    #doctest.testmod()
    doctest.testmod()
    
  

