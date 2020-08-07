import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as mpltPath
from astropy.io import fits as pyfits
import copy
from matplotlib.patches import PathPatch#, PointPatch

def iso_box(fn):
    """
      Return path from points in fn
      
      Parameters
      ----------
        iso_fn : str
            Filename containing the points for the iso-box
      
      Returns
      -------
        path, ppath : Path, PathPatch instances
      
    """
    dd = np.genfromtxt(fn, delimiter=',')
    xp, yp = dd[0], dd[1]   
    path = mpltPath.Path(np.array([xp, yp]).T, closed=True)
    ppath = PathPatch(path, fill=False)
    return path, ppath

def within_iso(points, iso_path, verbose=1):
    """
      Return boolean array indicating if the points are within the iso-box

      Parameters
      ----------
        points : array, shape: (N, 2)
          shall be color (BP-RP) and abs_mag (Gmag)
        iso_path : Path instance
    """
    if type(iso_path) == str:
        if verbose>0: print("iso_tools::within_iso - Reading iso-box from \'%s\'..." % iso_path)
        path, *x = iso_box(iso_path)
    return path.contains_points(points, radius=0)

def add_iso_column(ifn, ofn, iso_fn="aux_data/iso.dat", \
                   columns={"bp_rp":"bp_rp", "parallax":"parallax", "Gmag":"phot_g_mean_mag"},\
                   outcol="iso_compatible", extension=1, verbose=1, overwrite=False):
    """
      Add column to fits-file indicating if an object is within the iso-box.

      Parameters
      ----------
        ifn, ofn : str
          In- and out-filenames
        iso_fn : str
          Filename for the isochrone-box
        columns : dictionary
          Mapping of required properties to fits-file column names, must contain "bp_rp", "parallax","Gmag"
        extension : int
          Extension in fits-file
        outcol : str
          Name for fits-column indicating compatibility with isochrone
          
    """
    if verbose>0: print("iso_tools::add_iso_column - Reading ",ifn)
    ff = pyfits.open(ifn)
    color = ff[extension].data[columns["bp_rp"]]
    dist = 1000./ff[extension].data[columns["parallax"]]
    #psig = ff[1].data["parallax"] / ff[1].data["parallax_error"]
    mag = ff[extension].data[columns["Gmag"]]
    abs_mag = mag - 5*np.log10(dist)+5
    points = np.array([color, abs_mag])
    iso = within_iso(points.T, iso_fn)
    if verbose>0: print('iso_tools::add_iso_column - %i from %i are within iso-box (%5.2f%%).' % (np.sum(iso), len(iso), 100.*np.sum(iso)/len(iso)))
    
    colnames = [c.name for c in ff[extension].columns]     
    hdu = pyfits.PrimaryHDU()    
    col = pyfits.Column(name=outcol, array=iso, format="L")
    if outcol in colnames:
        ff[extension].data[outcol] = iso
    else:
        ff[extension].columns.add_col(col)
    hdu = copy.copy(ff[0]) 
    cols = ff[extension].columns
    hdx = pyfits.BinTableHDU.from_columns(cols)
    hdul = pyfits.HDUList([hdu, hdx])
    hdul.writeto(ofn, overwrite=overwrite)
    ff.close()
    if verbose>0: print('iso_tools::add_iso_column - Written: ',ofn)
    
def check_iso_box(iso_fn="iso.dat", parsec_fn="parsec_isochrones.dat",\
                  columns={"bp_rp":"bp_rp", "parallax":"parallax", "Gmag":"phot_g_mean_mag", "iso":"iso"},\
                  extension=1, data=None):
    """
      Plot isochron box
      
      Parameters
      ----------
        iso_fn : str
            Filename containing the points for the iso-box
        parsec_fn : str
            Filename for the PARSEC isochrones (http://stev.oapd.inaf.it/cgi-bin/cmd)
        data : str (optional)
            Fits-filename containing the data
        columns : dictionary
            Mapping of fits-column names (for \'data\')
        extension : int
            Fits-file extension (for \'data\')
    """
    # Iso-box:
    #dd = np.genfromtxt(iso_fn, delimiter=',')
    #xp, yp = dd[0], dd[1]   
    #path = mpltPath.Path(np.array([xp, yp]).T, closed=True)
    path, ppath = iso_box(iso_fn)
    plt.gca().add_patch(ppath)
    
    # Data (if provided)
    if data is not None:
        ff = pyfits.open(data)
        color = ff[extension].data[columns["bp_rp"]]
        dist = 1000./ff[extension].data[columns["parallax"]]    
        mag = ff[extension].data[columns["Gmag"]]
        abs_mag = mag - 5*np.log10(dist)+5
        points = np.array([color, abs_mag])
        plt.scatter(points[0],points[1], s=10, alpha=0.5)
        iso = within_iso(points.T, iso_fn)
        plt.scatter(points[0][iso],points[1][iso], color='k', s=13, alpha=0.2)


    # Isochrones
    dd = np.genfromtxt(parsec_fn, unpack=True)
    color = dd[-2] - dd[-1]
    mag = dd[-3]
    age = dd[2]
    mass = dd[5]
    plt.scatter(color, mag, c=age, vmin=5, vmax=10, s=5)

    # Within
    points = np.array([color, mag])
    iso = path.contains_points(points.T, radius=0)
    plt.scatter(color[iso], mag[iso], s=23, edgecolor='r', facecolor='', alpha=0.3)
    
    plt.xlabel("Bp-Rp")
    plt.ylabel("Gmag")
    plt.ylim(40,-15)
    plt.show()
    
if __name__ == "__main__":
    #check_iso_box()
    #add_iso_column("eFEDS_gaia.fits", "eFEDS_gaia_iso.fits", overwrite=True)
    check_iso_box(data="eFEDS_gaia_iso.fits")
    
