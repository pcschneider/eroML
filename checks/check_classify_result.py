import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u

props = ["offset_sig","log_plx", "logFxFg","bp_rp", "pos", "log_sk"]


tfn = "../../ero_data/training_eFEDS.fits"
mfn = "../classify/major_proba.fits"
rfn = "../classify/random_proba.fits"

norm=False
for p in props:
    for what, fn in zip(["random", "training", "real"],[rfn, tfn, mfn]):
        ff = pyfits.open(fn)
        if what=="random":
            y, xx = np.histogram(ff[1].data[p], density=norm, bins=40)
            x = (xx[1:] + xx[0:-1])/2
            line, = plt.plot(x, y/10,label="random all", ls='--')
            line.set_drawstyle("steps-mid")
            gi = np.where(ff[1].data["predicted"] == 0)[0]
            print(what, len(gi))
            y, xx = np.histogram(ff[1].data[p][gi], density=norm, bins=xx)
            gi = np.where((ff[1].data["predicted"] == 0) & (ff[1].data["NN"]==1))[0]
            print("    NN==1: ", len(gi))
            x = (xx[1:] + xx[0:-1])/2
            line, = plt.plot(x, y/10,label=what, color=line.get_color())
            line.set_drawstyle("steps-mid")
        elif what=="training":
            gi = np.where(ff[1].data["category"] == 0)[0]
            print(what, len(gi))
            y, xx = np.histogram(ff[1].data[p][gi], density=norm, bins=xx)
            gi = np.where((ff[1].data["category"] == 0) & (ff[1].data["NN"]==1))[0]
            print("    NN==1: ", len(gi))
            x = (xx[1:] + xx[0:-1])/2
            line, = plt.plot(x, y,label=what)
            line.set_drawstyle("steps-mid")            
        elif what=="real":
            gi = np.where(ff[1].data["predicted"] == 0)[0]
            print(what, len(gi))
            y, xx = np.histogram(ff[1].data[p][gi], density=norm, bins=xx)
            gi = np.where((ff[1].data["predicted"] == 0) & (ff[1].data["NN"]==1))[0]
            print("    NN==1: ", len(gi))
            x = (xx[1:] + xx[0:-1])/2
            line, = plt.plot(x, y,label=what)
            line.set_drawstyle("steps-mid")
    plt.xlabel(p)    
    plt.legend()
    plt.show()
#exit()    

fn = "../classify/major_proba.fits"
ff = pyfits.open(fn)
gi = np.where(ff[1].data["predicted"]==0)[0]
print("Predicted eligible: ",len(gi))
srcIDs = np.unique(ff[1].data["srcID"][gi])
print(len(srcIDs))

#for i in gi:
    #NN = ff[1].data["NN"][i]
    #if NN>1:
        #osrcID = ff[1].data["srcID"][i][0:7]
        #print(ff[1].data["srcID"][i], osrcID)
        #si = np.where(ff[1].data["srcID"] == osrcID)[0]
        #print(ff[1].data["predicted"][si])
        #print()
fn = "../classify/random_proba.fits"
#fn = "../../ero_data/training_eFEDS.fits"
ff = pyfits.open(fn)

MM = 4
if "category" in ff[1].data.columns.names:
    gi = np.where((ff[1].data["category"]==0) & (ff[1].data["NN"]<MM) )[0]
    plt.hist(ff[1].data["pos"][gi])
    plt.xlabel("pos")
    plt.show()             
    gi = np.where((ff[1].data["category"]==0) & (ff[1].data["NN"]<MM)  & (ff[1].data["pos"]<110))[0]
else:
    gi = np.where((ff[1].data["predicted"]==0) & (ff[1].data["NN"]<MM)  & (ff[1].data["pos"]<110))[0]
    

dist = 1000/10**ff[1].data["log_plx"][gi]
iii = np.where(dist<1000)[0]
print(len(gi), len(iii))
gi = np.where((ff[1].data["predicted"]==0) & (ff[1].data["NN"]<MM)  & (ff[1].data["pos"]<110))[0]
dist = 1000/10**ff[1].data["log_plx"][gi]
iii = np.where(dist>1000)[0]
print(len(gi), len(iii))
gi = np.where((ff[1].data["predicted"]==0) & (ff[1].data["NN"]<MM)  & (ff[1].data["pos"]<35))[0]
dist = 1000/10**ff[1].data["log_plx"][gi]
iii = np.where(dist>1000)[0]
print(len(gi), len(iii))



nbins=20
logbins = np.logspace(np.log10(10),np.log10(5000),nbins)
plt.hist(dist, bins=logbins)
plt.xlabel("d (pc)")
plt.ylabel("N")
plt.xscale("log")
plt.show()


plt.scatter(ff[1].data["pos"], ff[1].data["skd"])
plt.scatter(ff[1].data["pos"][gi], ff[1].data["skd"][gi])
plt.xlabel("pos")
plt.ylabel("skd")
plt.show()


plt.scatter(ff[1].data["bp_rp"], ff[1].data["logFxFg"])
plt.scatter(ff[1].data["bp_rp"][gi], ff[1].data["logFxFg"][gi])
plt.xlabel("bp_rp")
plt.ylabel("logFxFg")
plt.show()

plt.scatter(ff[1].data["pos"], ff[1].data["logFxFg"])
plt.scatter(ff[1].data["pos"][gi], ff[1].data["logFxFg"][gi])
plt.xlabel("pos")
plt.ylabel("logFxFg")
plt.show()


plt.scatter(ff[1].data["pos"], ff[1].data["log_plx"])
plt.scatter(ff[1].data["pos"][gi], ff[1].data["log_plx"][gi])
plt.xlabel("pos")
plt.ylabel("log_plx")
plt.show()






