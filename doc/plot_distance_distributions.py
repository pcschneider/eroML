import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0,10, 1000)
sig = 2.
y = x/sig**2 * np.exp(-x**2/(2*sig**2))

yy = 1 - np.exp(-sig**-2 * x**2 / 2)
#print(yy)
tmp = np.trapz(y, x=x)
print("sum under curve: ",tmp)

rnd = np.random.rand(1000)
print(rnd)

xx = np.sqrt(2*sig**2 * (-np.log(1-yy)))

rndx = np.sqrt(2*sig**2 * (-np.log(1-rnd)))

#print(xx)
plt.title("Distance distribution for real matches")
plt.plot(x, y, lw=2)

plt.hist(rndx, range=(0,10), bins=50, normed=True, color='g',alpha=0.5)
plt.xlabel("Distance (arcsec)")
plt.ylabel("pdf", color='b')


tx = plt.twinx()
tx.set_ylabel("cdf", color='r')


tx.plot(x, yy, lw=2, color='r')

#plt.plot(yy,xx)
plt.axvline(x=sig, color='0.5')
plt.xlim(0,10)
plt.ylim(0,1.05)
plt.savefig("distances_real_sources.png")
plt.show()

#\frac{x}{\sigma^2} e^{\frac{-x^2}{2 \sigma^2}}
