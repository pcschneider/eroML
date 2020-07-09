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
        self.mapper = {}
        self.uid_name = "uid"
        self._N = 0
        self.deepcopy = deepcopy
        
    def add_object(self, obj):
        """
        Add an object to `Ensemble`
        
        Parameters
        ----------
        obj : Astro_Object
        """
        if type(obj) != Astro_Object:
            raise TypeError("`Ensemble`::add_object - `obj` must be `Astro_Object` instance!")
        
        uid = self._N + 1
        if self.deepcopy:
            self.objects[uid] = copy.deepcopy(obj)
        else:
            self.objects[uid] = obj
        self.mapper[obj.srcID] = uid    
        self._N+=1

    def __len__(self):
        """
        Number of objects in `Ensemble`
        """
        return len(self.objects)
        
    def del_object(self, obj_name):
        """
        Remove object from `Ensemble`
        """
        if obj_name in self.mapper:
            del self.objects[self.mapper[obj_name]]
            del self.mapper[obj_name]
        else:
            raise IndexError(str("%s not in Ensemble." % obj_name))

    def skyCoords(self, srcIDs=None, epoch=2000, verbose=5):
        """
        The SkyCoords for the objects in the `Ensemble`
        
        Parameters
        ----------
        srcIDs : iterable (like list)
            Return only the coordinates for the objects in srcIDs
        epoch : float
            The epoch for the coordinates
        
        Returns
        -------
        Coordinates : SkyCoord instance
        """
        
        if srcIDs is None:
            srcIDs = self.mapper.keys()
        if verbose>3: print("astro_ensemble::Ensemble::skyCoords - Getting coords for ",srcIDs)    
        ra, dec = [], []
        for o in srcIDs:
            c = self.objects[self.mapper[o]].coord_tuple(ra_unit="degree", dec_unit="degree",epoch=epoch)
            ra.append(c[0])
            dec.append(c[1])
        return SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
        
    def array(self, colnames=(), make_array=True):
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
            print(c)
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
    
    
  
