import numpy as np
import matplotlib.pyplot as plt

dd = np.genfromtxt("offs.dat", unpack=True)

sk = np.unique(dd[1])
for s in sk:
    gi = np.where(dd[1] == s)[0]
    plt.hist(dd[0][gi], bins=60, range=(0,30),label=str(s))

plt.legend()
plt.show()
