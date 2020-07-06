import copy
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

class Astro_object():
    """
    A basic and flexible class that holds the properties of an astrophysical object.
    The minimum requirement is that the object has a source ID (srcID). Most importantly,
    however, an astrophysical object also has a coordinate and,
    potentially, a proper motion. See example below.
    
    Example
    -------
    >>> a = Astro_object({"srcID":"a",\ 
                          "coord":SkyCoord(ra=150*u.degree, dec=45*u.degree),\
                          "pm":(10,-20)}, pm_name="pm")
    >>> a.srcID = "b"
    >>> print(a.srcID)
    b
    >>> round(a.coord_tuple(epoch=3000)[0], 5)
    150.00393
    
    """
    
    def __init__(self, dct, id_name = "srcID", coord_epoch=2000, coord_name="coord", pm_name=None, deepcopy=True):
        """
        
        The minimal object creation is        
        >>> a = Astro_object({"srcID":"test"})
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
        self.coord_epoch = coord_epoch
        self.coord_name = coord_name
        self.pm_name = pm_name
        self.id_name = id_name
        
        if id_name not in dct.keys():
            raise LookupError("Each AstroObject must have an id")
        
        if deepcopy:
            self.dct = copy.deepcopy(dct)
        else:    
            self.dct = dct
            
    @property
    def srcID(self):
        return self.dct[self.id_name]        
    
    @srcID.setter
    def srcID(self, value):
        self.dct[self.id_name] = value
             
    @property
    def coord(self):
        return self.dct[self.coord_name]        
    
    def coord_tuple(self, ra_unit="degree", dec_unit="degree", epoch=2000):
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
        
    def coord4epoch(self, epoch=2000):
        """
        Calculate the coordinate for a different epoch
        
        Returns
        -------
        new coordinate : SkyCoord
        """
        dt = float(epoch - self.coord_epoch)
        if self.pm_name is None:
            raise LookupError("AstroObject has no entry for proper motion.")
        
        px  = self.dct[self.pm_name][0]/1000 * dt
        py = self.dct[self.pm_name][1]/1000 * dt
        
        offset = np.sqrt(px**2+py**2)
        angle = np.arctan2(px, py) / np.pi*180.
    
        return self.dct[self.coord_name].directional_offset_by(angle*u.degree, offset*u.arcsec)
               
        
if __name__ == "__main__":
    
    import doctest
    #doctest.testmod()
    doctest.testmod(extraglobs={'a': Astro_object({"srcID":"a", "coord":SkyCoord(ra=150.0*u.degree, dec=20.0*u.degree), "pm":(10,-50)}, pm_name="pm")})
    
  
