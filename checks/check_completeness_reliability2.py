import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits as pyfits
plt.rcParams.update({'font.size': 16})

def p_stellar(ids, probs, prior=0.075):
    uids = np.unique(ids)
    print(max(probs))
    r = np.zeros(len(ids))
    for ident in uids:
        tmp = np.in1d(ids, ident)
        p = probs[tmp]
        top = np.sum(p)
        bottom = top+(1-prior)
        r[tmp] = top/bottom
        r[tmp] = max(p)
        #print(p)
        
    return r#/max(r)
    
#def p_ij(ids, probs, prior=0.075):
    #uids = np.unique(ids)
    #print(max(probs))
    #r = np.zeros(len(ids))
    #for ident in uids:   
        #tmp = np.in1d(ids, ident)
        #p = probs[tmp]
    
def N4prob(fd, co=0.5):
        gi = np.where(fd["svm_prob"] > co)[0]
        ui = fd["original_srcID"][gi]
        return len(np.unique(ui))

bfn = "../ero_data/eFEDS_c001_main_ctp_star_1.0.fits"
bfn = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"
bff = pyfits.open(bfn)
bfd = bff[1].data

imfn = "checks/seb.png"
im = plt.imread(imfn)

fn = "major_eFEDS_classified_HamStar.fits"
ff = pyfits.open(fn)
fd = ff[1].data

rfn = "random_eFEDS_classified_HamStar.fits"
rff = pyfits.open(rfn)
rfd = rff[1].data

si = np.argsort(fd["svm_prob"])[::-1]
#plt.plot(fd["svm_prob"][si])
#plt.show()

a, b = [], []
cnt = 0
for p in fd["svm_prob"][si]:
    Nreal = N4prob(fd, p)
    Nrand = N4prob(rfd, p)
    a.append(Nreal)
    b.append(Nrand/10.)
    print(p)
    cnt+=1
    if p<0.1: break
a = np.array(a)
b = np.array(b)
    
#plt.imshow(im, extent=(1,0.2,0.4,1.0), aspect="auto")
    
pp = [ 1.98372460, -4.53855281,  3.53929735, -4.44042056e-03]
p = np.polyval(pp, fd["svm_prob"][si][0:len(a)])

fig = plt.figure(figsize=(8,4))
fig.subplots_adjust(top=0.98, right=0.98, left=0.08, bottom=0.15)

plt.plot(p/0.99, (a-b)/2060, label="SVC completeness", lw=2)
plt.plot(p/0.99, (a-b)/a, label="SVC reliability", lw=2)

bb = np.cumsum(1-p/0.99)
plt.plot(p/0.99, (a-bb)/a, ls=":")


bfn = "../ero_data/eFEDS_c001_main_ctp_star_1.0.fits"
bfn = "../ero_data/efeds_c001_V3_main_HamStar_internal.fits"
imfn = "checks/seb.png"
im = plt.imread(imfn)

#plt.imshow(im, extent=(1,0.2,0.4,1.0), aspect="auto")

bff = pyfits.open(bfn)
bfd = bff[1].data

uids = np.unique(bfd["ero_ID"])
eID = bfd["ero_ID"].tolist()
#print(uids)
ui = [eID.index(a) for a in uids]
print(bfn, " - unique eroIDs: ",len(ui), " catalog entries: ", len(bfd["ero_ID"]))

si = np.argsort(bfd["p_stellar"][ui])

pscaling = 0.88
pscaling = 1.0

spurious =  np.cumsum(1-bfd["p_stellar"][ui][si[::-1]]*pscaling)
missed = np.cumsum(bfd["p_stellar"][ui][si]*pscaling)
N = np.arange(len(si))

#print(bfd["p_stellar"][ui][si])
#print((N-spurious)/2060)
#print(max(spurious), max(missed), np.sum(bfd["p_stellar"][ui][si]))

#plt.plot(bfd["p_stellar"][ui][si], missed, label="missed")
#plt.plot(bfd["p_stellar"][ui][si[::-1]], spurious, label="spurious")
#plt.plot(bfd["p_stellar"][ui][si[::-1]], N, label="N")
plt.plot(bfd["p_stellar"][ui][si[::-1]], (N-spurious)/N, label="Bayes reliability", lw=2)
plt.plot(bfd["p_stellar"][ui][si][::-1], (N-spurious)/(N-spurious+missed[::-1]), label="Bayes completeness", lw=2)
#plt.plot(bfd["p_stellar"][ui][si][::-1], (N-spurious)/2092., label="completeness2")
plt.xlabel(r"$> p_{stellar}$")


plt.legend()
plt.xlim(0.97, 0.35)
plt.ylim(0.4,1.01)
#plt.plot()
plt.show()

pstar = p_stellar(fd["original_srcID"], fd["svm_prob"], prior = 0.075)

print(np.sum(np.unique(pstar)))


plt.scatter(fd["svm_prob"], pstar)
plt.show()

uprbs = np.unique(pstar)
si = np.argsort(uprbs)
plt.plot(uprbs[si])

spurious =  np.cumsum(1-uprbs[si[::-1]])
missed = np.cumsum(uprbs[si])[::-1]
N = np.arange(len(si))

print(max(missed))
#print(uprbs)

plt.show()


plt.plot(uprbs[si][::-1], missed, label="missed")
plt.plot(uprbs[si[::-1]], spurious, label="spurious")
plt.plot(uprbs[si[::-1]], N, label="N")
plt.legend()
plt.show()

plt.imshow(im, extent=(1,0.2,0.4,1.0), aspect="auto")

plt.plot(uprbs[si[::-1]], (N-spurious)/N, label="reliability")
plt.plot(uprbs[si][::-1], (N-spurious)/(N-spurious+missed), label="completeness")
plt.plot(uprbs[si][::-1], (N-spurious)/2060., label="completeness2")
plt.plot([0,1],[0.9,0.9])
plt.legend()
plt.show()

exit()

plt.imshow(im, extent=(1,0.2,0.4,1.0), aspect="auto")





uids = np.unique(bfd["ero_ID"]).tolist()
#print(uids)
ui = [uids.index(a) for a in uids]
print(bfn, " - unique eroIDs: ",len(ui), " catalog entries: ", len(bfd["ero_ID"]))

si = np.argsort(bfd["p_stellar"][ui])

pscaling = 0.88
pscaling = 1.0

spurious =  np.cumsum(1-bfd["p_stellar"][ui][si[::-1]]*pscaling)
missed = np.cumsum(bfd["p_stellar"][ui][si]*pscaling)
N = np.arange(len(si))

print(bfd["p_stellar"][ui][si])
print((N-spurious)/2060)
print(max(spurious), max(missed))

plt.plot(bfd["p_stellar"][ui][si], missed, label="missed")
plt.plot(bfd["p_stellar"][ui][si[::-1]], spurious, label="spurious")
plt.plot(bfd["p_stellar"][ui][si[::-1]], N, label="N")
#plt.plot(bfd["p_stellar"][ui][si[::-1]], (N-spurious)/N, label="reliability")
#plt.plot(bfd["p_stellar"][ui][si], (N-spurious)/(N-spurious+missed), label="completeness")
#plt.plot(bfd["p_stellar"][ui][si][::-1], (N-spurious)/2350., label="completeness2")


nn = np.where(fd["NN"] == 1)[0]
uids = fd["original_srcID"][nn]
print("my IDs: ",len(uids))
pp = [ 1.98372460, -4.53855281,  3.53929735, -4.44042056e-03]

svm_probs = np.polyval(pp, fd["svm_prob"])

plt.xlabel("p_stellar")
plt.ylabel("N")
plt.ylim(0.4, 1.05)
plt.xlim(1, 0.2)
plt.legend()
plt.show()

