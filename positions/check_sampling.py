import numpy as np
import matplotlib.pyplot as plt

xx = np.linspace(0,12,100)
a = 3.2
plt.plot(xx, xx/72)
N=100000

rnd = 12*np.random.rand(N)**0.5
d = rnd
plt.hist(d, density=True, bins=30)

plt.show()


plt.plot(xx, 0.0017*xx**2)
N=100000

rnd = 12*np.random.rand(N)**(1/3)
d = rnd
plt.hist(d, density=True, bins=30)

plt.show()




plt.plot(xx, 0.00019*xx**3)
N=100000

rnd = 12*np.random.rand(N)**(1/4)
d = rnd
plt.hist(d, density=True, bins=30)

plt.show()
