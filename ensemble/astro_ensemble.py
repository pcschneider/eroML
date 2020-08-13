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


    def split(self, N, verbose=1):
        """
        Split Ensemble into `N` sub-Ensembles
        
        Not all (sub-) Ensembles will have the same size, the last (sub-) Ensemble will contain the ramining objects.
        
        Parameters
        ----------
        N : int
            Number of sub-Ensembles.
        """
        r = []
        step = len(self)//N
        si = np.random.permutation(range(len(self)))
        i0, i1 = 0, step
        
        for i in range(N):
            #print(i, " - ",i0, i1)
            if i==N-1:
                i1 = len(self)
            f = copy.deepcopy(self)
            sids = np.array(f.srcIDs())[si[i0:i1]]
            f.keep(sids, verbose=0)
            r.append(f)
            i0=i1
            i1+=step
            
        return r    


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
    

    def rename(self, name, newname):
        """
        Rename object `name` to `newname`
        
        Parameters
        ----------
        name, newname : str
        """
        if newname in self.row_mapper:
            raise ValueError("Ensemble::rename - Object with name " + str(newname) + " already in Ensemble.")
                
        ridx = self.row_mapper[name]
        self.array["srcID"][ridx] = newname
        del self.row_mapper[name]
        self.row_mapper[newname] = ridx
        
        if name in self.mapper:
            idx = self.mapper[name]
            self.objects[idx].srcID = newname
            del self.mapper[name]
            self.mapper[newname] = idx
            

    def add_col(self, colname, array, force=True):
        """
        Add or overwrite column (if a column with this name already exists)
        
        Parameters
        ----------
        array : array, must be same length as len(Ensemble)
        name : str, Name of the column
        force : boolean
            Overwrite extsting content
        """
        #import numpy.lib.recfunctions as rfn
        #print("Adding ",colname, self.known_cols)
        if colname in self.known_cols: 
            if force:
                return self.set_col(colname, array)
            else:
                raise ValueError(str("Ensemble::add_col - Column with name \'%s\' alrady known." % colname))
        #a = rfn.append_fields(a, 'USNG', np.empty(a.shape[0], dtype='|S100'), dtypes='|S100')
        #dt = array.dtype
        #self.array = rfn.append_fields(self.array, colname, array, dtypes=dt)
        #print(dir(self.array.dtype.names))
        dt = []
        for n in self.array.dtype.names:
            #print(n, self.array[n].dtype)
            dt.append((n, self.array[n].dtype))
        dt.append((colname, array.dtype))
        tmp_arr = np.zeros(len(self), dtype=dt)
        for j, c in enumerate(self.array.dtype.names):
            tmp_arr[c] = self.array[c]
        tmp_arr[colname] = array
        self.array = tmp_arr
        self.known_cols.append(colname)
        

    def set_col(self, colname, array):
        """
        Set values (entire array column) for one column
        """
        if colname not in self.known_cols:
            raise IndexError(str("Ensemble::set_col - Column with name %s not known." % colname))
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
        
        if "pm_RA" in self.known_cols:
            dtime = epoch - self.array["ref_epoch"]
            d_ra = self.array["pm_RA"] * dtime
            d_dec = self.array["pm_Dec"] * dtime
            
            #TODO: calc pos angle and add to coordinates
        
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
        

        oa_arr = other.to_array(colnames=merge_cols, array_type='recarray')
        oa = oa_arr[idx]
        
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
            ow = "match_dist"
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
            if verbose>5: print("__getitem__ - Creating Astro_Object", name," -> row-index ",row_idx)
            dct = {}
            for c in self.known_cols:
                if c=="RA": continue
                if c=="Dec": continue
                dct[c] = self.array[c][row_idx]
            
            if idx is not None:
                for c in self.objects[idx].dct.keys():
                    #print("key: ",c)
                    dct[c] = self.objects[idx][c]
            
            coord = SkyCoord(ra=self.array["RA"][row_idx], dec=self.array["Dec"][row_idx], unit=(u.degree, u.degree))
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
                        dt = (c, '<U36')
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
                            #print(i, "ccc",c, "(",len(self),np.shape(self.array),len(np.unique(self.srcIDs())),")")
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


    def append(self, other, postfix=None, cols="all", verbose=1):
        """
        Append another Ensemble
        
        Parameters
        ----------
        other : Ensemble
        postfitx : str
            append `postfix` to `srcID` if object with the same `srcID` already exists in `self`.
        cols : str or tuple
            How to handle columns existing in only one Ensemble
            - `all` : keep all columns and use np.nan for the entries in the other Ensemble
            - `same` : keep only columns existing in both Ensembles
            - `left` : keep all columns of left Ensemble
        """
        srcIDs1, srcIDs2 = self.srcIDs(), other.srcIDs()
        shared_cols = np.intersect1d(self.known_cols, other.known_cols)
        N0, N1 = len(self), len(other)
        
        if verbose>5:
            print("Ensemble::append - shared_cols: ",shared_cols)

        # Get columns for new array
        if type(cols) == type("all"):
            if cols == "all":
                outcols = np.unique(self.known_cols + other.known_cols).tolist()
            elif cols == "same":
                outcols = np.intersect1d(self.known_cols, other.known_cols).tolist()
            elif cols == "left":
                outcols = self.known_cols
            else:
                raise NotImplementedError(str("Ensemble::append - I don't understand cols=%s..." % cols))
        elif isinstance(cols, tuple) or isinstance(cols, list):
            outcols = cols
        else:
            raise NotImplementedError(str("Ensemble::append - I don't understand cols=%s..." % str(cols)))
        if verbose>5: print("Ensemble::append - Using cols=",outcols)

        # Build dtype
        dt = []
        for c in outcols:
            if c in self.known_cols:
                dtt = (c, self.array[c].dtype)
            elif c in other.known_cols:
                dtt = (c, other.array[c].dtype)
            else:
                print("Argh~!")
                exit()
            dt.append(dtt)
        if verbose > 2: print("Ensemble::append - Using dtype=", dt)   
        
        # Resolve ID conflicts:
        shared_ids = np.intersect1d(srcIDs1, srcIDs2)
        if verbose>5: print("Ensemble::append - shared_ids", shared_ids)
        if len(shared_ids)>0 and postfix is None:
            raise ValueError("Ensemble::append - Name conflicts, but no `postfix` given.")
            
        for si in shared_ids:
            other.rename(si, si+postfix)

        # Construct new array
        new_arr = np.zeros(N0+N1, dtype=dt)
        for c in outcols:
            left  = self.array[c] if c in self.known_cols else np.array(N0*[np.nan])
            right = other.array[c] if c in other.known_cols else np.array(N1*[np.nan])
                
            arr = np.hstack((left, right))
            new_arr[c] = arr

        # Assign new array to this Ensemble
        self.known_cols = outcols
        #print("outcols: ",type(outcols), outcols)
        self.array = new_arr
        
        for i, nn in enumerate(other.srcIDs()):
            self.row_mapper[nn] = i+N0


def fake_ensemble(N=10, random_pos=True, ID_prefix="src_", center=(0, 0), width=(1, 1), pm=False, randomIDs=False, ref_epoch=None, seed=None, random_cols=[], verbose=2):
    """
    Generate fake ensemble, either with random or uniform coordinates.
    
    The sky density will be approximately uniform. Only for large regions, the 
    RA coordinates will be skewed somewhat, because only `dec_center` is 
    considered.
    
    Example
    -------
    >>> from eroML.ensemble import random_ensemble
    >>> e = random_ensemble(center=(10,45), pm=10, seed=1)
    >>> e["src_0"].RA
    
    Parameters
    ----------
    N : int 
       Number of objects in Ensemble
    random_pos : boolean
       Use random coordinates within coordinate-box
    center : tuple (RA, Dec), float
       Coordinate center (in degree)
    ra_width, dec_width : float
       Width in degree (i.e., the region size will be width_RA x width_Dec)
    pm : float
       Add random proper motion (values are -pm...+pm in both, pm_RA and pm_Dec)
       Unit: mas/yr
    ref_epoch : float
       Epoch for pos (decimal-year)
    random_cols : list of str
       Include columns with names `random_cols` and random entries
    seed : int
       Control seed for numpy random (useful if a specific output is needed)
        
    Returns
    --------
    Ensemble
    """
    np.random.seed(seed=seed)
    
    if ref_epoch is None:
        ref_epoch = 2000.
    
    if randomIDs:
        ids = [ID_prefix+str(i) for i in np.random.randint(1000000, size=N)]
    else:
        ids = [ID_prefix+str(i) for i in np.arange(N)]
    
    if random_pos:
        dec = (np.random.rand(N)-0.5) * width[1] + center[1]
        cos_ra_min = (center[0]-width[0]) * np.cos(center[1]/180*np.pi)
        cos_ra_max = (center[0]+width[0]) * np.cos(center[1]/180*np.pi)
        ra = ((np.random.rand(N)-0.5) *  (cos_ra_max - cos_ra_min)) / 2 / np.cos(center[1]/180*np.pi) + center[0]
    else:
        dec = np.linspace(center[1]-width[1]/2, center[1]+width[1]/2, N)
        ra  = np.linspace(center[0]-width[0]/2, center[0]+width[0]/2, N)
    
    if pm:
        pm_ra, pm_dec = (np.random.rand(N)-0.5)*2*pm, (np.random.rand(N)-0.5)*2*pm
        dt = [("srcID", type("X"),16), ("RA", float), ("Dec", float), ("pm_RA", float), ("pm_Dec", float), ("ref_epoch", float)]
        rec = np.zeros(N, dtype=dt)
        rec["pm_RA"] = pm_ra
        rec["pm_Dec"] = pm_dec
        rec["ref_epoch"] = ref_epoch
    else:
        dt = [("srcID", '|S25'), ("RA", float), ("Dec", float)]
        rec = np.zeros(N, dtype=dt)
        
    print(max(ra), min(ra))
    print(max(dec), min(dec))    
        
    rec["srcID"] = ids
    rec["RA"] = ra
    rec["Dec"] = dec

    e = Ensemble()
    e.from_array(rec)
   
    for cc in random_cols:
        x = np.random.rand(N)
        e.add_col(cc, x)
   
    #if verbose>1: print(rec)    

    return e
    
    
    
if __name__ == "__main__":
    pass
    
    
  
