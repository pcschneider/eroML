"""
A collection of astronomical objects, called an Ensemble.
"""
import copy
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from .astro_object import Astro_Object
from collections import OrderedDict

class Ensemble():
    """
    Holds a number of  `Astro_Objects`
    
    Minimum requirement for an entry are
      - srcID (str)
      - coord (SkyCoord); internally, RA and Dec are stored in the array (in degree)
    If proper motion is present, the following columns must exist
      - pm_RA (float, unit mas/year)
      - pm_Dec (float, unit mas/year)
      - ref_epoch (float, decimal year)
    """
    def __init__(self, deepcopy=True):
        self.objects = OrderedDict()
        self.mapper = OrderedDict()
        self.uid_name = "uid"
        self._N = 0
        self.deepcopy = deepcopy
        self.known_cols = ["srcID","RA", "Dec"]
        
        self.array = np.recarray((0,),dtype=[("srcID",str), ("RA",float), ("Dec",float)])
        
    def add_one_object(self, obj, auto_resolve=True, verbose=1):
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
            raise TypeError("`Ensemble`::add_one_object - `obj` must be `Astro_Object` instance!")
        
        if type(obj.srcID)!=str:
            obj.srcID = str(obj.srcID).strip()
        if obj.srcID in self.mapper:
            if auto_resolve: 
                sID = create_alternate_ID(obj.srcID)
            else:
                raise KeyError(str("`Ensemble`::add_one_object - An object with srcID=%s already in Ensemble." % obj.srcID.strip()))
        else:
            sID = obj.srcID
          
        # Check if and which cols are not in Ensemble.array  
        object_cols = list(obj.dct.keys())
        non_array_cols = []
        known_none_cols = []
        for oc in object_cols:
            if oc not in self.known_cols:
                if oc == "coord": continue
                non_array_cols.append(oc)
                
        for kc in self.known_cols:
            if kc == "RA":
                continue
            elif kc == "Dec":
                continue
            elif kc not in object_cols:
                known_none_cols.append(kc)
                
        if verbose>3: print("Ensemble::add_one_object - non_array_cols: ",non_array_cols," known_none_cols: ",known_none_cols, " cols", self.known_cols)
        
        array_entry = []
        for c in self.known_cols:
            if c=="RA":
                array_entry.append(obj.dct["coord"].ra.degree)
            elif c=="Dec":
                array_entry.append(obj.dct["coord"].dec.degree)
            elif c in known_none_cols:
                array_entry.append(None)
            else:
                array_entry.append(obj.dct[c])
        array_entry = tuple(array_entry)
        self.array.resize(len(self)+1)
        self.array[-1] = array_entry
                
        if len(self) != self._N:
            raise IndexError("Ensemble::add_one_object - self.len should be equal to self._N!")
        
        if len(non_array_cols) > 0:    
          uid = self._N + 1
          if self.deepcopy:
              self.objects[uid] = copy.deepcopy(obj)
          else:
              self.objects[uid] = obj
          self._N+=1
          
        self.mapper[obj.srcID] = uid    
        
        #print(sID, self.mapper)
        return sID
    
    def add_col(self, array, name=None):
        """
        Add or overwrite column (if a column with this name already exists)
        
        Parameters
        ----------
        array : array, must be same length as len(Ensemble)
        name : str, Name of the column
        """
        for i,o in enumerate(self.mapper.values()):
            self.objects[o].dct[name] = array[i]
        if name not in self.known_cols: self.known_cols.append(name)
        
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
            N = self.mapper[obj_name.strip()]
            srcIDs = self.srcIDs()
            #for i in range(self._N - N):
                #self.mapper[srcIDs[i+N]]-=1
            del self.objects[self.mapper[obj_name.strip()]]
            del self.mapper[obj_name.strip()]
            #print(self.mapper)
            #self._N-=1
        else:
            raise IndexError(str("%s not in Ensemble." % obj_name))
        print("del",self.mapper)

    def keep(self, srcIDs, verbose=10):
        """
        Keep only given sources in Ensemble 
        
        Parameters:
        -----------
        srcIDs : array of str
        """
        n_mapper = OrderedDict()
        n_objects = OrderedDict()
        n_old = len(self)
        
        for i, si in enumerate(srcIDs):
            
            n_objects[i+1] = self.objects[self.mapper[si]]
            n_mapper[si] = i+1
        self.mapper = n_mapper
        self.objects = n_objects
        
        if verbose>1: print("Ensemble::keep - Keeping only ",len(self), " from ",n_old," objects.")
        
    def skyCoords(self, srcIDs=None, epoch=2000, verbose=1):
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
        N0 = len(self)
        N1 = len(array)
        
        
        if clean: 
            self.array = array
            self.known_cols = list(names)
            for i, srcID in enumerate(self.array["srcID"]):
                self.mapper[srcID.strip()] = i
        else:
            self.array.resize(N0+N1)    
            self.array[N0:N1] = array
            new_known = []
            for c in known_cols:
                if c in names:
                    new_known.append(c)
            self.known_cols = new_known
            for i, srcID in enumerate(array["srcID"]):
                self.mapper[srcID.strip()] = self._N+i+1
  
        #print(self.mapper)
            
            
        #for o in array:
            #dct = {n:o[n] for n in names}
            #ra, dec = o["RA"], o["Dec"]
            #coord = SkyCoord(ra, dec, unit=(u.degree, u.degree))
            #dct["coord"] = coord
            ##print(dct.keys())
            ##print()
            #tmp = Astro_Object(dct)
            #self.add_one_object(tmp)
            
        return self    
   
    def srcIDs(self):
        """
        The source IDs
        
        Returns
        -------
        srcIDs - list
        """
        return list(self.mapper.keys())
    
    def merge_add(self, other, conflict_resolution="append", col_postfix="_NN", criterium="distance", verbose=10, **kwargs):
        """
        Add properties of another Ensemble to this Ensemble. Also adds property `match_dist`
        
        Parameters
        ----------
        other : Ensemble
          The other Ensemble
        conflict_resolution : str
          How to handle columns with the same name. Options are ["append", "left", "right"]
            - append : append `col_postfix` to column-name and append column
            - left : use column in this Ensemble
            - right : use column in other Ensemble
        col_postfix : str
          Append this to column name in case `conflict_resolution` == 'append'
        criterium : str
          Must be "distance" at the moment.
          
          Additional `kwargs`: 
            - "NN" : int, nth nearest neighbor
            - "epoch" : float, epoch to match
        """
        if criterium != "distance":
            raise NotImplementedError("Ensemble::merge_add - `criterium` must be `distance` at the moment.")
        
        if "NN" in kwargs:
            NN = kwargs["NN"]
        else:
            NN = 1
        
        if "epoch" in kwargs:
            epoch = kwargs["epoch"]
        else:
            epoch=None
        
        coord0 = self.skyCoords(epoch=epoch)
        coord1 = other.skyCoords(epoch=epoch)
        oIDs = other.srcIDs()
        
        if verbose>1: print("Ensemble::merge_add - NN:",NN) 
        
        idx, d2d,d3d = coord0.match_to_catalog_sky(coord1, nthneighbor=NN) 
        cols = other.known_cols
        for c in cols:
            #print(c)
            if c=="coord": continue
            if c in self.known_cols:
                if conflict_resolution=="append":
                    ow = c+col_postfix
                    self.known_cols.append(ow)
                elif conflict_resolution=="right":    
                    ow = c
                elif conflict_resolution=="left":
                    continue
                for j, k in enumerate(self.mapper.keys()):
                    #print(j, idx[j], oIDs[idx[j]])
                    self.objects[self.mapper[k]].dct[ow] = other.objects[oIDs[idx[j]]][c]

        if NN>1:
            ow = "match_dist_"+str(NN)
        else:
            ow = "match_dist"

        self.known_cols.append(ow)
        for j, k in enumerate(self.mapper.keys()):
            #print(j, idx[j], oIDs[idx[j]])
            self.objects[self.mapper[k]].dct[ow] = d2d[j].arcsec
            
    def __getitem__(self, name, verbose=1):
        """
        """
        #print("Getting",name)
        #print("Name: ",name)
        #print(self.mapper.keys())
        idx = self.mapper[name]
        #print("idx: ",idx)
        #print("XX", self.objects.keys())
        #return 1
        if idx not in self.objects.keys():
            print("__getitem__ - Creating Astro_Object", name," -> ",idx)
            dct = {}
            for c in self.known_cols:
                if c=="RA": continue
                if c=="Dec": continue
                dct[c] = self.array[c][idx]
            coord = SkyCoord(ra=self.array["RA"][idx], dec=self.array["Dec"][idx], unit=(u.degree, u.degree))
            dct["coord"] = coord
            aa = Astro_Object(dct, id_name='srcID', ref_epoch=2000, coord_name='coord', pm_name=None)
            return aa
        else:    
            print("__getitem__ - returning Astro_Object", name," -> ",idx)
            return self.objects[idx]
    
    def to_array(self, colnames=(), array_type="recarray", verbose=1):
        """
        Usually used to get a numpy.array or np.recarray from given `colnames`.
        
        Returns
        -------
        array if `array_type` == "recarray", np.array if `array_type`=="array" else dictionary
            array shape is (len(colnames), N) with N being the number of objects in Ensemble
        """
        colnames = np.atleast_1d(colnames)
        
        non_array_cols = []
        for c in colnames:
            if c not in self.known_cols:
                non_array_cols.append(c)
        if verbose>5: print("to_array: non_array_cols:",non_array_cols)        
        if len(non_array_cols)==0:
            if array_type=="recarray":
                return copy.deepcopy(self.array)
            elif array_type=="array":
                arr = []
                for c in colnames:
                    arr.append(self.array[c].astype(float))
                if len(colnames)==1: return np.array(arr).flatten()    
                return np.array(arr)    
            elif array_type=="dict":
                dct = {}
                for c in colnames:
                    dct[c] = self.array[c]
                return dct    
            else:
                raise LookupError("Ensemble::array - `array_type` must be in [recarray, array, dict], but is " +str(array_type))
        
        else:
            return None
                
                
        if array_type=="recarray":
            cc = np.core.records.fromarrays([dct[n] for n in dct], names=",".join(dct.keys()))
            return cc #np.array(list([dct[c] for c in colnames]))
        elif array_type=="array":
            arr = []
            for c in dct.keys():
                arr.append(np.array(dct[c]).astype(float))
            if len(colnames)==1: return np.array(arr).flatten()    
            return np.array(arr)    
        elif array_type=="dict":
            return dct
        else:
            raise LookupError("Ensemble::array - `array_type` must be in [recarray, array, dict], but is " +str(array_type))
        
         
        
if __name__ == "__main__":
    pass
    
    
  
