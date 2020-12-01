import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import factorial 
from PyAstronomy import funcFit as fuf


class dModel(fuf.OneDFit):
    """
    Implements a straight line of the form y = "off" + x * "lin".
    """

    def __init__(self):
        fuf.OneDFit.__init__(self, ["r_scale", "y_scale", "eta", "r_shift"])
        self["r_scale"] = 1
        self["y_scale"] = 1
        self["eta"] = 1
        self["r_shift"] = 0
        
    def evaluate(self, x):
        """
        Calculates and returns model according to the \
        current parameter values.

        Parameters:
        - `x` - Array specifying the positions at \
                which to evaluate the model.
        """
        y = self["y_scale"] * 2*np.pi*eta*(self["r_scale"] * x + self["r_shift"]) * np.exp(-np.pi*self["eta"]*(self["r_scale"] * x + self["r_shift"])**2)
        y[y<0] = 0
        return y


r = np.linspace(0,30, 1000)
eta0 =0.01
mm = dModel()
mm["eta"] = eta0
mm.thaw(["r_scale", "y_scale", "r_shift"])

s = np.zeros(len(r))
for i in range(1,20):
    n=i
    eta = eta0
    pdf = 2*(np.pi*eta)**n*r**(2*n-1)/(factorial(n-1))*np.exp(-np.pi*eta*r**2)
    l, = plt.plot(r, pdf, label=str(n))
    mm.fit(r, pdf)
    mm.parameterSummary()
    #n = 1
    #eta = eta0/i
    #pdf = 2*(np.pi*eta)**n*r**(2*n-1)/(factorial(n-1))*np.exp(-np.pi*eta*r**2)
    plt.plot(r, mm.evaluate(r), color=l.get_color(), ls='--')
    
    s+=pdf
plt.plot(r,s/8, color='k', lw=2)    
plt.legend()    
plt.show()
