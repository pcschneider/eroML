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
      
    The data structure is organized in an array, which contains the data available for all 
    objects and is kept up-to-date when adding or removing objects (may take time). Only if
    objects contain properties not generally available (check `known_cols`-attribute), 
    a copy of the object is stored and can be retrieved.
    With this approach, access to the `known_cols`-properties is fast.
    
    
    """
    def __init__(self, deepcopy=True):
        self.objects = OrderedDict()
        self.mapper = OrderedDict()
        self.row_mapper = OrderedDict()
        self.uid_name = "uid"
        self._N = 0
        self._Nrow = 0
        self.deepcopy = deepcopy
        self.known_cols = ["srcID","RA", "Dec"]
        
        self.array = np.recarray((0,),dtype=[("srcID",str), ("RA",float), ("Dec",float)])
     
    def filter_array_for_srcIDs(self, array, srcIDs, verbose=10):
        """
        Restrict array to those entries that match the given srcIDs
        
        Parameters
        ----------
        array - array, to be filtered
        srcIDs - list/np.array, the srcIDs
        
        Returns
        -------
        array, restricted to srcIDs
        """
        idx = list([self.row_mapper[s] for s in srcIDs])
        return array[idx]
        
    def add_one_object(self, obj, auto_resolve=True, verbose=1):
        """
        Add an object to `Ensemble`.
        
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
        NRows = len(self)
        self.array.resize(NRows+1)
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
        self.row_mapper[obj.srcID] = NRows
        
        #print(sID, self.mapper)
        return sID
    
    def add_col(self, colname, array):
        """
        Add or overwrite column (if a column with this name already exists)
        
        Parameters
        ----------
        array : array, must be same length as len(Ensemble)
        name : str, Name of the column
        """
        import numpy.lib.recfunctions as rfn
        if colname in self.known_cols: 
            raise ValueError(str("Ensemble::add_col - Column with name %s alrady known." % colname))
        #a = rfn.append_fields(a, 'USNG', np.empty(a.shape[0], dtype='|S100'), dtypes='|S100')
        dt = array.dtype
        self.array = rfn.append_fields(self.array, colname, array, dtypes=dt)
        self.known_cols.append(colname)
        
    def set_col(self, colname, array):
        """
        """
        if colname not in self.known_cols:
            raise IndexError(str("Ensemble::set_col - Column with name %s alrady known." % colname))
        if len(array) != len(self):
            raise ValueError(str("Ensemble::set_col - Column with name %s  should have %i instead of %i entries." % (colname, len(self), len(array))))
                             
        self.array[colname] = array
        
    def __len__(self):
        """
        Number of objects in `Ensemble`
        """
        return len(self.row_mapper)

    def shift_array(self, dN, n0=None, n1=None):
        """
        
        Shift (part of) the rows of the array.
        
        Parameters
        ----------
        dN : int
            The step-width, i.e., shift
        n0, n1 : int
            The start and end indices for shifting
        """
        
        N = len(self.array)
        
        if n0 is not None:
            i0 = n0+dN
            if n1 is not None:
                i1 = n0+dN, n1+dN
            else:
                n1 = N
                i1 = n1+dN
        else:
            i0=0
            i1 = N
            n0 = 0
            n1 = N
            
        self.array[i0:i1] = np.roll(self.array, dN)[n0:n1]        
        self.array.resize(N+dN)
            
    def rebuild_row_mapper(self, dN, n0=None):
        """
        Make sure that the row_mapper reflects the array.
        
        Parameters
        ----------
        n0 : int
            if provided, the array is assumed to correct up to uid `n`
        """
        if n0 is None:
            n0 = -1
        N = len(self.mapper)
        srcIDs = self.srcIDs()
        for si in srcIDs:
            if self.row_mapper[si] > n0:
                self.row_mapper[si]+=dN
        
    def del_object(self, obj_name, verbose=1):
        """
        Remove object from `Ensemble`
        """
        if verbose>5: print("Ensemble::del_object - Removing ",obj_name)
        if obj_name in self.mapper:
            N = len(self)            
            uid = self.mapper[obj_name.strip()]
            del self.mapper[obj_name.strip()]
            if uid in self.objects:
                del self.objects[uid]
            n0 = self.row_mapper[obj_name]
            if verbose>6:
                print("Ensemble::del_object - Removing ",obj_name," which has row index: ",n0, " and uid: ",uid)
            del self.row_mapper[obj_name]    
            self.shift_array(-1, n0=n0)
            self.rebuild_row_mapper(-1, n0=n0)
            if verbose>6:
                print("Ensemble::del_object - Now len(self)= ", len(self))
        else:
            raise IndexError(str("%s not in Ensemble." % obj_name))
        #print("del",self.mapper)

    def keep(self, srcIDs, verbose=10):
        """
        Keep only given sources in Ensemble 
        
        Parameters:
        -----------
        srcIDs : array of str
        """
        n_mapper = OrderedDict()
        n_row_mapper = OrderedDict()
        n_objects = OrderedDict()
        n_old = len(self)
        #print(self.mapper)
        narr = np.zeros(len(srcIDs), dtype=self.array.dtype)
        for i, si in enumerate(srcIDs):
            #print(i, si)
            narr[i] = self.array[self.row_mapper[si]]
            n_row_mapper[si] = i
            if si in self.mapper:
                n_objects[i+1] = self.objects[self.mapper[si]]
                n_mapper[si] = i+1
        self.mapper = n_mapper
        self.row_mapper = n_row_mapper
        self.objects = n_objects
        self.array = narr
        
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
        ra, dec = self.to_array(colnames=["RA","Dec"], array_type="array")
        if verbose>5: print("astro_ensemble::Ensemble::skyCoords - #coords",len(ra))    
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
                sid = str(srcID).strip()
                #self.mapper[sid] = i
                self.row_mapper[sid] = i
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
        return self    
   
    def srcIDs(self):
        """
        The source IDs
        
        Returns
        -------
        srcIDs - list
        """
        return list(self.row_mapper.keys())
    
    def merge_add(self, other, conflict_resolution="append", col_postfix="_NN", criterium="distance", verbose=1, **kwargs):
        """
        Add properties of another Ensemble to this Ensemble. Also adds property `match_dist`
        
        .. todo::
        
            Also add objects that are in other ensemble. Currently, only the other's array is considered.
        
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
        import numpy.lib.recfunctions as rfn

        if criterium != "distance":
            raise NotImplementedError("Ensemble::merge_add - `criterium` must be `distance` at the moment.")
        
        if "NN" in kwargs:
            NN = kwargs["NN"]
        else:
            NN = 1
        if verbose>1: print("Ensemble::merge_add - NN:",NN) 

        if "epoch" in kwargs:
            epoch = kwargs["epoch"]
        else:
            epoch=None
        
        coord0 = self.skyCoords(epoch=epoch)
        coord1 = other.skyCoords(epoch=epoch)
        oIDs = other.srcIDs()
        
        idx, d2d,d3d = coord0.match_to_catalog_sky(coord1, nthneighbor=NN) 
        
        cols = other.known_cols
        merge_cols = []
        ows = []
        
        # Identify the columns to be merged
        for c in cols:
            if c=="coord": continue
            if c == "srcID":
                cc = "srcID"+col_postfix
                oa = other.to_array(colnames=c)[idx]
                dt = oa.dtype[0]
                self.array = rfn.append_fields(self.array, cc, oa, dtypes=dt)
                self.known_cols.append(cc)                
                continue
        
            if c in self.known_cols:
                if conflict_resolution=="append":
                    ow = c+col_postfix
                    self.known_cols.append(ow)
                elif conflict_resolution=="right":    
                    ow = c
                elif conflict_resolution=="left":
                    continue
            else:
                ow = c
                self.known_cols.append(ow)                
            ows.append(ow)
            merge_cols.append(c)
        
        oa = other.to_array(colnames=merge_cols, array_type='recarray')[idx]
        
        # Build new merged array
        # First, the dtype
        dt = oa.dtype
        ca = []
        for i in range(len(ows)):
            ca.append(oa[merge_cols[i]])
        fdt = []
        for i in range(len(self.array.dtype)):
            fdt.append((self.array.dtype.names[i], str(self.array.dtype[i])))
        for i in range(len(dt)):
            fdt.append((ows[i], str(dt[i])))
        # Second, create and fill the new array
        xxx = np.zeros(len(self), dtype=fdt)        
        for j, c in enumerate(self.array.dtype.names):
            xxx[c] = self.array[c]        
        for j, c in enumerate(oa.dtype.names):
            xxx[ows[j]] = oa[c]
        # Last, overwrite existing array
        self.array = xxx
        
        
        if NN>1:
            ow = "match_dist_"+str(NN)
        else:
            ow = "match_dist"
        ows.append(ow)
        self.known_cols.append(ow)
        self.array = rfn.append_fields(self.array, ow, d2d.arcsec, dtypes=["f4"])
        
        #for ow in ows:
            #self.known_cols.append(ow)
            
        if verbose>5: print("Ensemble::merge_add - added cols: ",ows) 
        if verbose>3: print("Ensemble::merge_add - all known cols: ",self.known_cols)
        if verbose>0: print("Ensemble::merge_add - Added ", len(ows)+1, " columns to Ensemble.")
         
    def __getitem__(self, name, verbose=1):
        """
        """
        #print("Getting",name)
        #print("Name: ",name)
        #print(self.mapper.keys())
        row_idx = self.row_mapper[name]
        if name in self.mapper:
            idx = self.mapper[name]
        else:
            idx = None
        #print("idx: ",idx)
        #print("XX", self.objects.keys())
        #return 1
        if idx not in self.objects.keys():
            if verbose>5: print("__getitem__ - Creating Astro_Object", name," -> ",idx)
            dct = {}
            for c in self.known_cols:
                if c=="RA": continue
                if c=="Dec": continue
                dct[c] = self.array[c][row_idx]
            
            if idx is not None:
                for c in self.objects[idx].dct.keys():
                    #print("key: ",c)
                    dct[c] = self.objects[idx][c]
            
            coord = SkyCoord(ra=self.array["RA"][idx], dec=self.array["Dec"][idx], unit=(u.degree, u.degree))
            dct["coord"] = coord
            aa = Astro_Object(dct, id_name='srcID', ref_epoch=2000, coord_name='coord', pm_name=None)
            return aa
        else:    
            if verbose>5: print("__getitem__ - returning Astro_Object", name," -> ",idx)
            return self.objects[idx]
    
    def to_array(self, colnames=(), srcIDs=None, array_type="recarray", verbose=1):
        """
        Usually used to get a numpy.array or np.recarray from given `colnames`.
        
        Parameters
        ----------
        colnames : str, or list/tuple
            The column names from which to construct the array
        srcIDs : list (iterable in general)
            If not None: Return only values for given srcIDs
        array_type : str, in ["recarray", "array", "dict"]    
        
        Returns
        -------
        Content of self.array, the array shape is (len(colnames), N) with N being the number of objects in Ensemble : The return value depends on `array_type`, see above
          - recarray if `array_type` == "recarray"
          - np.array if `array_type` == "array" 
          - dictionary if `array_type` == "dict"
        """
        colnames = np.atleast_1d(colnames)
        
        non_array_cols = []
        for c in colnames:
            if c not in self.known_cols:
                non_array_cols.append(c)
        
        if verbose>4: print("Ensemble::to_array: Getting data for cols=",colnames)
        if verbose>5: print("Ensemble::to_array: non_array_cols:",non_array_cols)        
        
        if len(non_array_cols)==0:
            #if array_type=="recarray":
                #if srcIDs is not None:
                    #return self.filter_array_for_srcIDs(copy.deepcopy(self.array), srcIDs)
                #else:
                    #print("XX")
                    #aa = self.array[colnames]
                    #r = copy.copy(aa    )
                    #print("yy")
                    #return r
            if array_type=="array" or array_type=="recarray":
                arr = []
                dts = []
                #print("XXXXXXXXXXXXX")
                for c in colnames:
                    tmp = self.array[c]
                    #print("len(tmp)", len(tmp))
                    dt = (c, str(tmp.dtype))
                    if "i" not in dt[1] and "f" not in dt[1]: 
                        #print(c, type(tmp), dt)
                        tmp = np.array(tmp).astype(str)
                        dt = (c, '<U10')
                    elif "f" in dt[1]:
                        idx = np.isfinite(tmp)
                        #print(idx)
                        tmp[idx == False] = np.nan
                    elif "i" in dt[1]:
                        pass
                    
                    #if c=="astrometric_excess_noise": continue
                    #print(c, dt)
                    arr.append(np.array(tmp).astype(dt[1]))
                    #dt.name = c
                    dts.append(dt)
                #print()    
                if array_type=="recarray":
                    #print(dts, np.shape(arr))
                    if len(colnames)>1:
                        xxx = np.transpose(arr)
                        xxx = np.zeros(len(self), dtype=dts)
                        for i, c in enumerate(colnames):
                            xxx[c] = arr[i]
                        return xxx
                    return np.array(np.array(arr).flatten(), dtype=dts)
                
                if len(colnames)==1: 
                    if srcIDs is not None:
                        self.filter_array_for_srcIDs(np.array(arr), srcIDs).flatten()
                    else:
                        return np.array(arr).flatten()
                if srcIDs is not None:
                    return self.filter_array_for_srcIDs(np.array(arr), srcIDs)
                else:
                    return np.array(arr)
                
            elif array_type=="dict":
                dct = {}
                if srcIDs is not None:
                    idx = list([self.row_mapper[x] for x in srcIDs])
                else:
                    idx = range(len(self))
                for c in colnames:
                    dct[c] = self.array[c][idx]
                return dct    
            else:
                raise LookupError("Ensemble::array - `array_type` must be in [recarray, array, dict], but is " +str(array_type))
        
        else:
            raise LookupError("Ensemble::array -cannot find columns %s in array." % str(non_array_cols))
                
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
    
    
  
