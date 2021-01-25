from eroML.utils import NN_distribution
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
import PyAstronomy.funcFit as fuf

class Dist_model2(fuf.OneDFit):
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



class Dist_model(fuf.OneDFit):
    """
    A model for the expected distance between real and random sources.
    """

    def __init__(self):
        fuf.OneDFit.__init__(self, ["fraction","sig_scale","dens_scale"])
        self["fraction"] = 0.07
        self["sig_scale"] = 1.
        self["dens_scale"] = 1.
        self.setRestriction({"fraction":[0,1]})
          
    def evaluate(self, x):
        M = len(x)
        xx = x.repeat(self.N).reshape((M,self.N))
        #yy = x/sigma**2*np.exp(-x**2/(2*sigma**2))
        
        sigma = self.sigs*self["sig_scale"]
        zz = xx/sigma**2*np.exp(-xx**2/(2*sigma**2))

        skd = self.dens*self["dens_scale"]
        sigma = 1/np.sqrt(2*np.pi*skd)
        yy =  xx/sigma**2*np.exp(-xx**2/(2*sigma**2))
        
        mm = (1-self["fraction"])*np.sum(yy, axis=1) + self["fraction"] * np.sum(zz, axis=1)
        return mm


fn_r = "../../ero_data/random_eFEDS.fits"
ffr = pyfits.open(fn_r)
ffrd = ffr[1].data

N_plot = 3
fig = plt.figure()
ax0 = fig.add_subplot(N_plot, 1,1)
fn = "../../ero_data/major_eFEDS.fits"   
#fn = '../../ero_data/training_eFEDS.fits'
fn = "../../ero_data/major_eFEDS_EDR3.fits"   
#../../ero_data/major_eFEDS_EDR3.fits
#fn = fn_r
Nsig = 1000
fkey="match_dist"
errkey="sigma_r"
#N = NN_distribution(, Nsig=1000, ext=1, fkey="match_dist", errkey="sigma_r", ofn="check.pdf", verbose=1)
ext = 1
ff = pyfits.open(fn)
ffd = ff[1].data
gi = np.where(ff[ext].data["NN"] == 1)[0]

md = ff[ext].data[fkey][gi]
N = len(md)
print("Found %i sources in file." % N)
mn_skd = np.mean(ff[ext].data["skd"][gi])
print("Mean sky density: ",mn_skd)
Nrnd = np.zeros(N)
x = np.linspace(0.5, 150.5, 151).repeat(N).reshape((151,N))
#y = np.zeros((N, len(x)))
#print(np.shape(x))
skd = ff[ext].data["skd"][gi] / 3600 / 10 * 1.07
#skd = np.ones(len(gi))
y = np.transpose(skd[:,None] * x.T *2* np.pi)
print(np.shape(y))

sigma = 1/np.sqrt(2*np.pi*skd)
xs = np.linspace(0,100,200)
yy = x/sigma**2*np.exp(-x**2/(2*sigma**2))
print(np.shape(yy))

sigma = ffd["sigma_r"][gi]*0.6
zz = x/sigma**2*np.exp(-x**2/(2*sigma**2))


gi = np.where(ffd["NN"] == 1)[0]
#bns = np.histogram(ffd["match_dist"][gi], range=(0, 150), bins=151)
bns = np.histogram(ffd["match_dist"][gi], range=(0, 100), bins=101)
#print(bns[1])
bnx = (bns[1][1:] + bns[1][0:-1])/2

bny = bns[0]
#print(bns)
print(len(bnx), len(bny))
ax0.bar(bnx, height=bny)
#plt.plot(x[:,0],np.sum(yy, axis=1)*0.95)
#plt.plot(x[:,0],np.sum(yy, axis=1)*0.925)
#plt.plot(x[:,0],np.sum(yy, axis=1)*0.90)
#plt.plot(x[:,0],np.sum(zz, axis=1)*0.1, ls ='--')

frac = 0.4
mm = (1-frac)*np.sum(yy, axis=1) + frac * np.sum(zz, axis=1)
print(len(mm))
print("model N: ",np.trapz(mm, x=x[:,0]))
print("data N: ",np.sum(bny))
print("N stars: ",frac*np.trapz(mm, x=x[:,0]))
#ax0.plot(x[:,0],mm, ls ='--', color='r', lw=2, label="Model")
ax0.plot(x[:,0], np.mean(y, axis=1)*len(gi), color='y', label="Random (all)")

mdl = Dist_model()
mdl.N = len(gi)
mdl.sigs = ffd["sigma_r"][gi]
mdl.dens = ffd["skd"][gi] / 3600 / 10 
mdl["sig_scale"] = 0.6
mdl["dens_scale"] = 1.5
xs = x[:,0]
mdl.thaw(["fraction"])
mdl.fit(bnx, bny)
mdl.parameterSummary()

mdl.thaw(["dens_scale", "sig_scale"])
mdl.fit(bnx, bny)
print("N star from model: ",len(gi)*mdl["fraction"])

mdl.fit(bnx, bny, yerr=1+np.sqrt(bny+0.75))
mdl.parameterSummary()
tmp = mdl["fraction"]
ax0.plot(xs, mdl.evaluate(xs), ls='-', color='r', lw=2, label="model")
mdl["fraction"] = 0
ax0.plot(xs, (1-tmp)*mdl.evaluate(xs), ls='--', color='r', lw=1, label="Random")

mdl["fraction"] = 7.14*1e-2
ax0.plot(xs, mdl.evaluate(xs), ls=':', color='r', lw=1, label="Juergen")
mdl["fraction"] = tmp

print("N star from model: ",len(gi)*mdl["fraction"])
ax0.set_xlim(0,100)
ax0.set_ylim(0, 650)
ax0.set_ylabel("N")
ax0.legend()

ax1 = fig.add_subplot(N_plot, 1,2, sharex=ax0)
ax1.plot(bnx, bny - mdl.evaluate(bnx))
ax1.plot([0,150],[0,0], color='k')
ax1.set_ylabel("N (data - model)")
ax1.set_xlim(0,100)

#ax1.set_xlabel("Match Distance (arcsec)")
ax1.grid(True)

if N_plot > 3:
    ax2 = fig.add_subplot(313, sharex=ax0)
    ax2.plot(bnx, np.cumsum(bny - mdl.evaluate(bnx)))
    ax2.plot([0,150],[0,0], color='k')
    ax2.set_ylabel("N cumsum(Data - Model)")
    ax2.set_xlim(0,100)
    ax2.set_ylim(-100,100)
    ax2.set_xlabel("Match Distance (arcsec)")
    #ax2.grid(True)

if N_plot > 2:
    ax2 = fig.add_subplot(313, sharex=ax0)
    frc = mdl["fraction"]
    #plt.plot(bnx, mdl.evaluate(bnx), lw=2)
    mdl["fraction"] = 0
    rand = mdl.evaluate(bnx) * (1-frc)
    mdl["fraction"] = 1.0
    real = mdl.evaluate(bnx) * frc
    ratio = np.cumsum(rand) / np.cumsum(real)
    frac = np.cumsum(real) / np.sum(real)
    #ax2.plot(bnx, real)
    #ax2.plot(bnx, rand, ls=':')
    ax2.plot(bnx, ratio, label="Contamination")
    ax2.plot(bnx, frac, label="Completeness")
    #ax2.plot([0,150],[0,0], color='k')
    ax2.set_ylabel("Fraction")
    ax2.set_xlim(0,100)
    #ax2.set_ylim(-100,100)
    ax2.legend()
    ax2.set_xlabel("Match Distance (arcsec)")
    ax2.grid(True)


mdl["fraction"] = 0

for xlim in [10,15,20]:
    ii = np.where(bnx <= xlim)[0]
    print("Data - model (<",xlim,"): ",np.sum(bny[ii] - (1-tmp)*mdl.evaluate(bnx[ii])))

plt.show()



exit()


##plt.plot(x, x**2


#mm = Dist_model()
#mm.Nsig = Nsig
#mm.gen_sigs(ff[ext].data[errkey])



#mm["N"] = N
#mm["sig"] = 3.7496337890625 * 0.75
#mm["dens"] = 4.5e-4
#mm["fraction"] = 0
#mm["err_scaling"] = 0.8
#mm.thaw(["dens","N"])
##yerr = np.sqrt(bin_y)
#gi = np.where(bin_x > dlim)[0]

#mm.fit(bin_x[gi], bin_y[gi], maxfun=1e4, miniFunc="cash79", disp=disp)
#mm.parameterSummary()
#x = np.linspace(hist_range[0], hist_range[1],1000)

#if ofn is not None:
    #fig = plt.figure(figsize=(8,11))
    #ax0 = fig.add_subplot(311)
    #ax0.hist(md, range=hist_range, bins=hist_bins)
    #ax0.plot(x, mm.evaluate(x), ls='--')
    
    #ax = fig.add_subplot(312)
    ##ax.plot(bin_x, bin_y, ls='steps-mid')
    #ax.errorbar(bin_x, bin_y - mm.evaluate(bin_x), yerr=np.sqrt(bin_y))
    #ax.errorbar(bin_x[gi], bin_y[gi] - mm.evaluate(bin_x[gi]), yerr=np.sqrt(bin_y[gi]), color='r')
    #ax.axhline(y=0, color='k')
    #ax.set_ylim(-100,100)
    
##mm.parameterSummary()

#mm.freeze(["dens"])
#mm["fraction"] = 0.35
#mm["N"] = np.sum(hh[0])
#mm["err_scaling"]=0.75
##mm.thaw(["fraction", "err_scaling"])
#mm.thaw(["fraction"])
##mm.thaw(["fraction","sig"])
#mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79", disp=disp)
##mm.parameterSummary()

#mm.thaw(["dens"])
#mm.fit(bin_x, bin_y, maxfun=1e4, miniFunc="cash79", disp=disp)
##mm["dens"] = 0.000441403 * 10
        
#if verbose>0:
    #mm.parameterSummary()
    #print("Estimated number of sources: ",mm["N"], " instead of ", np.sum(hh[0]), " -> ",  np.sum(hh[0]) - mm["N"], "+-", np.sum(hh[0])**0.5,"(real-model)")

#ff.close()
#if ofn is not None:
    #ax = fig.add_subplot(313)
    ##ax.plot(bin_x, bin_y, ls='steps-mid')
    #ax0.plot(x, mm.evaluate(x), lw=2, color='r')
    #ff, NN = mm["fraction"], mm["N"]
    #mm["fraction"] = 1.
    #mm["N"] = ff*NN
    #ax0.plot(x, mm.evaluate(x))
    #mm["fraction"] = 0
    #mm["N"] = NN * (1-ff)
    
    #ax0.plot(x, mm.evaluate(x), ls='-')

    #mm["fraction"] = ff
    #mm["N"] = NN
    #ax.errorbar(bin_x, bin_y - mm.evaluate(bin_x), yerr=np.sqrt(bin_y))
    ##ax.errorbar(bin_x[gi], bin_y[gi] - mm.evaluate(bin_x[gi]), yerr=np.sqrt(bin_y[gi]), color='r')
    #ax.axhline(y=0, color='k')
    #ax.set_ylim(-100,100)
    #plt.savefig(ofn)
    #plt.show()

#print("N: ",mm["N"] * mm["fraction"])
