from astropy.io import fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import PyAstronomy.funcFit as fuf
from scipy.stats import invgauss, chi2, lognorm, recipinvgauss,skewnorm,norminvgauss

fn = "../ero_data/merged_major.fits"
ff = pyfits.open(fn)


#hh = plt.hist(ff[1].data["RADEC_ERR"], range=(0,10), bins=30, label="eRO")
#mu = 5
#loc = 1
#scale = 3
##rr = invgauss(mu, loc=loc, scale=scale)
##err = gengamma.rvs(5, 1.3, size=int(sum(hh[0])))
##plt.show()

#gi = np.where(np.isfinite(ff[1].data["RADEC_ERR"]))[0]
#a, b, loc, scale = norminvgauss.fit(ff[1].data["RADEC_ERR"][gi])
#print(a,loc, scale, len(gi), len(ff[1].data["RADEC_ERR"]))

#err = norminvgauss.rvs(a,b, scale=scale, loc=loc, size=int(sum(hh[0])))
##rv = invgauss(mu)
##x = np.linspace(0, 20, 100)
#plt.hist(err, range=(0,10), bins=30, alpha=0.5, label="model")
#plt.legend()
#plt.show()



#exit()

class Dist_model(fuf.OneDFit):
    """
    Implements a straight line of the form y = "off" + x * "lin".
    """

    def __init__(self):
        fuf.OneDFit.__init__(self, ["fraction","sig","dens","N", "err_scaling"])

    def evaluate(self, x):
        """

        """
        
        # Real matches:
        Nsig = 1000 # int(self["N"])
        #sig = norminvgauss.rvs(a,b, scale=scale, loc=loc, size=Nsig)
        sig = np.random.choice(ff[1].data["RADEC_ERR"]*self["err_scaling"], Nsig)
        #sig = int(self["N"]) * [self["sig"]]
        s2 = np.repeat(np.array(sig)**2, len(x))
        #print(np.shape(s2))
        #print("s2", s2)
        #s2 = self["sig"]**2
        tmpx = np.tile(x, Nsig)
        tmp0 = tmpx/s2 * np.exp(-tmpx**2/(2*s2)) * self["fraction"]
        tmp0 = tmp0.reshape((Nsig, len(x)))
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


def likeli_r(r, sig):
    """ Likelihood for source at distance r given uncertainty sig """
    s2 = sig**2
    return r/s2 * np.exp(-r**2/(2*s2))

hh = plt.hist(ff[1].data["match_dist"], range=(0, 80), bins=80)
print("Sources in hist: ",np.sum(hh[0]), " sources in fits-file: ",len(ff[1].data["srcID"]))

bin_x = 0.5*(hh[1][1::] + hh[1][0:-1])
bin_w = np.mean(hh[1][1::] - hh[1][0:-1])
print("Bin width in hist: ",bin_w)
print()
bin_y = hh[0]
#print(bin_x, bin_y)
x = np.linspace(0,80,1000)

mm = Dist_model()
mm["N"] = len(ff[1].data["match_dist"])
mm["sig"] = 2.5
mm["dens"] = 4.5e-4
mm["fraction"] = 0
mm["err_scaling"] = 0.8
mm.thaw(["dens","N"])
yerr = np.sqrt(bin_y)
gi = np.where(bin_x > 20)[0]

mm.fit(bin_x[gi], bin_y[gi], maxfun=1e4, miniFunc="cash79")
plt.plot(x, mm.evaluate(x))
mm.parameterSummary()

#print(np.sum(

mm.freeze(["dens"])
mm["fraction"] = 0.35
mm["N"] = np.sum(hh[0])
mm.thaw(["sig","fraction", "err_scaling"])
mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79")
mm.parameterSummary()

mm.thaw(["dens"])
mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79")
mm.parameterSummary()


plt.plot(x, mm.evaluate(x))

plt.yscale("log")

plt.ylim(min(bin_y)*0.8, max(bin_y)*1.2)
plt.show()
