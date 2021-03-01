from eroML.positions import analytic_probability
import matplotlib.pyplot as plt
import numpy as np

N = 1000
md = np.linspace(0.1, 15, N)
#md = np.linspace(2.999, 3.001, N)
sig = np.array(N*[3])
sig = np.linspace(2, 10, N)
skd = np.linspace(0.1, 10, N)
skd = np.linspace(0.9999, 1.0001, N)
#skd = np.array(N*[1])

#xv, yv = np.meshgrid(skd, md)
#pp = analytic_probability(match_dist=xv, sigma=sig, sky_density=yv)
#print(np.shape(pp))
#plt.imshow(pp, extent=(min(skd), max(skd), min(md), max(md)), aspect='auto', origin='lower')
#plt.xlabel("Sky density (stars/arcmin^2)")
#plt.ylabel("Match distance (arcsec)")
#plt.colorbar()
#plt.show()

xv, yv = np.meshgrid(sig, md)
pp = analytic_probability(match_dist=yv, sigma=xv, sky_density=0.8, ps=0.08)
plt.imshow(pp, extent=(min(sig), max(sig),min(md), max(md)), aspect='auto', origin='lower')
plt.gca().contour(pp, extent=(min(sig), max(sig),min(md), max(md)), origin='lower', colors='w')
plt.ylabel("Match distance (arcsec)")
plt.xlabel("Sigma (arcsec)")
cb = plt.colorbar()
cb.set_label("p_stellar")
plt.show()


#plt.plot(skd, pp)
#plt.show()
