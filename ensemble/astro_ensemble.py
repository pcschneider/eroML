import copy
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from .astro_object import Astro_Object

class Ensemble():
    """
    Holds a number of  `Astro_Objects`
    """
    def __init__(self, deepcopy=True):
        self.objects = {}
        self.uid_name = "uid"
        self.N = 0
        self.deepcopy = deepcopy
        
    def add_object(self, obj):
        """
        Add an object to `Ensemble`
        
        Parameters
        ----------
        obj : Astro_Object
        """
        #if type(obj) is not Astro_Object
        if type(obj) != Astro_Object:
            raise TypeError("`Ensemble`::add_object - `obj` must be `Astro_Object` instance!")
        
        uid = self.N + 1
        if self.deepcopy:
            self.objects[uid] = copy.deepcopy(obj)
        else:
            self.objects[uid] = obj
        self.N+=1
        
    def del_object(self, obj_name):
        """
        Remove object from `Ensemble`
        """
        pass
        
    def columns(self, colnames=(), make_array=True):
        """
        Usually used to get a numpy.array from given `colnames`.
        
        Returns
        -------
        array if `make_array` == True else dictionary
            array shape is (len(colnames), N) with N being the number of objects in Ensemble
        """
        dct = {}
        for c in colnames:
            dct[c] = []
            for o in self.objects.values():
                tmp = getattr(o,c)
                #print(o,c,tmp)
                dct[c].append(tmp)
        if make_array:
            return np.array(list([dct[c] for c in colnames]))
        else:
            return dct
        
if __name__ == "__main__":
    pass
    
    
  
