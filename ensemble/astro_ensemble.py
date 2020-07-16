"""
A collection of astronomical objects, called an Ensemble.
"""
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
        self.known_cols = []
        
    def add_object(self, obj, auto_resolve=True):
        """
        Add an object to `Ensemble`
        
        Parameters
        ----------
        obj : Astro_Object
        auto_resolve : Boolean
            Append appropriate `_i` to `srcID` if `auto_resolve` == True
        
        Returns
        -------
        uid : ID for object as stored, needed to retrieve exactly this object 
        """
        def create_alternate_ID(name):
            for i in range(1, 100):
                an = str("%s_%i" % (name, i+1))
                if an not in self.mapper:
                    return an
        
        if type(obj) != Astro_Object:
            raise TypeError("`Ensemble`::add_object - `obj` must be `Astro_Object` instance!")
        
        if obj.srcID.strip() in self.mapper:
          if auto_resolve: 
              sID = create_alternate_ID(obj.srcID.strip())
          else:
              raise KeyError(str("`Ensemble`::add_object - An object with srcID=%s already in Ensemble." % obj.srcID.strip()))
        else:
          sID = obj.srcID
          
        uid = self._N + 1
        if self.deepcopy:
            self.objects[uid] = copy.deepcopy(obj)
        else:
            self.objects[uid] = obj
        self.mapper[obj.srcID.strip()] = uid    
        self._N+=1
        
        for c in obj.dct.keys():
            if c not in self.known_cols:
                self.known_cols.append(c)
        
        return sID

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
            del self.objects[self.mapper[obj_name.strip()]]
            del self.mapper[obj_name.strip()]
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
    
    def from_array(self, array, verbose=1, clean=True):
        """
        Standardized population of `Ensemble`. 
        
        In particular, the array shall contain the following columns:
          - srcID
          - RA in degree
          - Dec in degree
        Proper motion is optional, but if they exist, all the following three column must exist:
          - pm_RA in mas/year
          - pm_Dec in mas/year
          - ref_epoch, e.g., 2000.0
        Paralllax is option, but if it exists, the naming is:
          - plx in mas
          
        Parameters
        ----------
        array : ndarray 
            shape: NxM with N being the number of objects and M the number of properties
        """
        
        names = array.dtype.names
        self.known_cols = names
        for o in array:
            dct = {n:o[n] for n in names}
            #print(dct.keys())
            #print()
            tmp = Astro_Object(dct)
            self.add_object(tmp)
        return self    
    
    
    def __getitem__(self, name, verbose=1):
        """
        """
        return self.objects[self.mapper[name]]
    
    def _export(self, verbose=1):
        """
        """
        pass
    
    def array(self, colnames=(), array_type="recarray"):
        """
        Usually used to get a numpy.array from given `colnames`.
        
        Returns
        -------
        array if `array_type` == "recarray" else dictionary
            array shape is (len(colnames), N) with N being the number of objects in Ensemble
        """
        dct = {}
        for c in colnames:
            dct[c] = []
            for o in self.objects.values():
                tmp = o[c]
                #print(o,c,tmp)
                dct[c].append(tmp)
        if array_type=="recarray":
            cc = np.core.records.fromarrays([dct[n] for n in dct], names=",".join(dct.keys()))
            return cc #np.array(list([dct[c] for c in colnames]))
        elif array_type=="array":
            arr = []
            for c in dct.keys():
                arr.append(np.array(dct[c]).astype(float))
            return np.array(arr)    
        elif array_type=="dict":
            return dct
        else:
            raise LookupError("Ensemble::array - `array_type` must be in [recarray, array, dict], but is " +str(array_type))
        
         
        
if __name__ == "__main__":
    pass
    
    
  
