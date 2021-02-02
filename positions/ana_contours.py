import numpy as np
import matplotlib.pyplot as plt
import PyAstronomy.funcFit as fuf



fn = "con.dat"
dd = np.genfromtxt(fn, unpack=True)
gi = np.where(dd[0] < 9.5)[0]

mm = fuf.PolyFit1d(3)
mm.thaw(["c0","c1","c2"])

plt.plot(dd[0], dd[1])
plt.plot(dd[0][gi], dd[1][gi])
mm.fit(dd[0][gi], dd[1][gi])
mm.parameterSummary()
x = np.linspace(0,14, 50)
plt.plot(x, mm.evaluate(x))
plt.show()
