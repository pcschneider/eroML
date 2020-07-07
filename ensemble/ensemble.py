import copy
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from .astro_object import Astro_object

class Ensemble():
    """
    """
    def __init__(self, deepcopy=True):
        self.objects = {}
        self.uid_name = "uid"
        self.N = 0
        self.deepcopy = deepcopy
        
    def add_object(self, obj):
        uid = self.N + 1
        if self.deepcopy:
            self.objects[uid] = copy.deepcopy(obj)
        else:
            self.objects[uid] = obj
        self.N+=1
        
    def gen_columns(self, colnames=[]):
        pass
        
if __name__ == "__main__":
    
    import doctest
    #doctest.testmod()
    doctest.testmod(extraglobs={'a': Astro_object({"srcID":"a", "coord":SkyCoord(ra=150.0*u.degree, dec=20.0*u.degree), "pm":(10,-20)}, pm_name="pm")})
    
  
