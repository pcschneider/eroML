import copy
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

class Astro_Object():
    """
    A basic and flexible class that holds the properties of an astrophysical object.
    The minimum requirement is that the object has a source ID (srcID). Most importantly,
    however, an astrophysical object also has a coordinate and,
    potentially, a proper motion. See example below.
    
    Example
    -------
    >>> a = Astro_Object({"srcID":"a","coord":SkyCoord(ra=150*u.degree, dec=45*u.degree),"pm":(10,-20)}, pm_name="pm")
    >>> a.srcID = "b"
    >>> print(a.srcID)
    b
    >>> round(a.coord_tuple(epoch=3000)[0], 5)
    150.00393
    
    """
    
    def __init__(self, dct, id_name = "srcID", ref_epoch=2000, coord_name="coord", pm_name=None, deepcopy=True):
        """
        
        The minimal object creation is        
        >>> a = Astro_Object({"srcID":"test"})
        >>> a.srcID
        'test'
    
        If *coordinates* are provided, they shall be a SkyCoord instance
        If *proper motion* values are provided, their unit is mas/year
    
    
        Parameters
        ----------
        dct : dictionary
        *_name : str
          describe the dictionary entries containing the relevant property
    
        """
        self.ref_epoch = ref_epoch
        self.coord_name = coord_name
        self.pm_name = pm_name
        self.id_name = id_name
        
        if id_name not in dct.keys():
            raise LookupError("Each Astro_Object must have an id")
        
        if deepcopy:
            self.dct = copy.deepcopy(dct)
        else:    
            self.dct = dct
     
    def __str__(self):
        return self.srcID
    
    @property
    def srcID(self):
        return self.dct[self.id_name]        
    
    @srcID.setter
    def srcID(self, value):
        self.dct[self.id_name] = value
             
    @property
    def coord(self):
        return self.dct[self.coord_name]        
    
    @property
    def ra(self):
        tmp = self.coord_tuple()
        return tmp[0]
    
    @property
    def dec(self):
        tmp = self.coord_tuple()
        return tmp[1]
    
    def coord_tuple(self, ra_unit="degree", dec_unit="degree", epoch=None):
        """
        
        Example
        --------
        >>> x = a.coord_tuple(ra_unit="hourangle", epoch=2100)
        >>> round(x[0], 7)
        10.0000197

        Returns
        -------
        (ra, dec) : tuple
        """
        c = self.coord4epoch(epoch=epoch)
        return getattr(c.ra, ra_unit),  getattr(c.dec, dec_unit)
        
    def coord4epoch(self, epoch=None):
        """
        Calculate the coordinate for a different epoch
        
        Returns
        -------
        new coordinate : SkyCoord
        """
        if epoch is None:
            epoch = self.ref_epoch
            dt = 0
        else:
            dt = float(epoch - self.ref_epoch)
            
        if self.pm_name is None:
            raise LookupError("Astro_Object has no entry for proper motion.")
        
        px  = self.dct[self.pm_name][0]/1000 * dt #/ np.cos(self.coord.dec.degree)
        py = self.dct[self.pm_name][1]/1000 * dt
        
        offset = np.sqrt(px**2+py**2)
        angle = np.arctan2(px, py) / np.pi*180.
    
        return self.dct[self.coord_name].directional_offset_by(angle*u.degree, offset*u.arcsec)

def from_Simbad(name, verbose=1):
    """
      Generate `Astro_Object` by querying Simbad
      
      Coordinates, and if available, proper motion as well as the parallax (plx) 
      are then available for the `Astro_Object` 
      
      Parameters
      ----------
      name : str
        Name as to be found in Simbad
        
      Returns
      -------
      Astro_Object
      
      Example
      -------
      >>> star = from_Simbad("AU Mic")
      >>> star.srcID
      'AU Mic'
    """
    from astroquery.simbad import Simbad
    Simbad.add_votable_fields('pmra')
    Simbad.add_votable_fields('pmdec')
    Simbad.add_votable_fields('plx')
    Simbad.add_votable_fields('ids')
    if verbose>5:
        print("Querying %s in Simbad..." % name)
    r = Simbad.query_object(name)
    ra, dec = r[0]["RA"], r[0]["DEC"]
    pmra, pmdec = r[0]["PMRA"]/1000, r[0]["PMDEC"]/1000
    plx = r[0]["PLX_VALUE"]
    if verbose>5:
        print("    Found ", len(r), "entries, using entry '0' with Coords=",ra, dec, ", pm=",pmra,pmdec," and plx=",plx)
    tmp = Astro_Object({"srcID":name, "coord":SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.degree)), "pm":(pmra,pmdec), "plx":plx}, pm_name="pm")
    return tmp
    

if __name__ == "__main__":
    
    import doctest
    #doctest.testmod()
    doctest.testmod(extraglobs={'a': Astro_Object({"srcID":"a", "coord":SkyCoord(ra=150.0*u.degree, dec=20.0*u.degree), "pm":(10,-50)}, pm_name="pm")})
    
  
