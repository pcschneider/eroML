from eroML.utils import *
from eroML.ensemble import from_fits
from .pixelize import add_healpix_col, extract_healpix_range
import glob

def loop(ero_fn=None, ofn=None, NSIDE=None, rID=1, skip=True):
    """
    Loop through tiles.
    
    Parameters
    ----------
    ero_fn : str
        filename for eROSITA sources list 
    ofn : str
        Source list with HPIX
    NSIDE : int (2**N)
        nside for healpix index calculation
    rID : int or str 
        run ID
    """

    print("using fn=",ero_fn, " ofn=",ofn," NSIDE=",NSIDE, "rID=",rID)
    #rID = 1
    hps = add_healpix_col(ero_fn, ofn=ofn, nside=NSIDE, overwrite=True)
    #hps = [500,501]
    #gi = np.where(hps==2814)[0][0]
    #print(gi)
    for i in hps:#[10:]:

        glob_str = "../ero_data/tile_hp"+str(i)+"_nside"+str(NSIDE)+"_rID"+str(rID)+"_training.fits"
        glob_str = "../ero_data/tile_hp"+str(i)+"_nside"+str(NSIDE)+"_rID"+str(rID)+"_gaia_sources.fits"
        fnames = glob.glob(glob_str)
        print(glob_str)
        print(fnames)
        if skip and len(fnames)>0: continue
        
        print()
        print()
        print(i)
        extract_healpix_range(ofn, "tmp_hp.fits", overwrite=True, min_index=i, max_index=i)
        t = Tile()
        t.populate_filenames("tmp_hp.fits", fn_prefix="../ero_data/tile_hp"+str(i)+"_nside"+str(NSIDE)+"_rID"+str(rID))
        t.prepare_data()
        #t.generate_sets()
        print()
        print()        


class Tile():
    """
    Orchestrates the analyzes of a sky region.
    
    Therer are four main objectives:
      1. Create the master catalog, which contains all the information needed to identify the stars
      2. Create a training catalog, based on well positionally matched sources
      3. Create a random catalog to assess to reliability of the sources
      4. Prepare the data for ML, i.e., convert/normalize columns etc.
    """
    def __init__(self, ero_fn=None, fn_prefix="Tile"):
        """
        Initialize Tile
        
        For parameters, see :func:`~eroML.tile.tile.Tile.populate_filenames`.
        """
        #self.e = None
        if ero_fn is not None:
            self.populate_filenames(ero_fn, fn_prefix=fn_prefix)            

    def populate_filenames(self, ero_fn, fn_prefix="Tile"):
        """
        Generate the filenames to be used later
        
        The following filenames will be available:
        
        - `ero_fn` : eROSITA source file
        - `gaia_fn` : File containing ALL Gaia sources
        - `major_fn` : Matched sources
        - `random_fn` : Shifted and matched sources
        - `training_fn` : subset of `major_fn` fullfilling certaing quality criteris
        - `training_random_fn` : training + random sources
        
        Parameters
        ----------
        ero_fn : str
            Filename for eROSITA source list
        fn_prefix : str
            String that will be prefixed when writing the corresponding files.
        """
        self.ero_fn = ero_fn
        self.gaia_fn = fn_prefix+"_gaia_sources.fits"
        self.major_fn = fn_prefix+"_major.fits"
        self.training_fn = fn_prefix+"_training.fits"
        self.training_random_fn = fn_prefix+"_training_random.fits"
        self.random_fn = fn_prefix+"_random.fits"
        
    def prepare_data(self, ero_fn=None, gaia_fn=None):
        """
        Add additional information to source lists.
        
        Parameters
        ----------
        ero/gaia_fn : str
            Filenames for eROSITA and Gaia source lists
        """
        
        if ero_fn is None: ero_fn = self.ero_fn 
        if gaia_fn is None: gaia_fn = self.gaia_fn     
        
        # First, download and enrich Gaia data
        get_gaia(ero_fn, gaia_fn, verbose=1)
        enrich_Gaia(gaia_fn)
        
        # Enrich eROSITA data
        e = enrich_eROSITA(ero_fn, mapper={"DETUID":"srcID", "DEC":"Dec"})
        
        
    def generate_sets(self, ero_fn=None, gaia_fn=None, training_rel_dist_cutoff=2, training_abs_dist_cutoff=3):       
        """
        Generate the relevant data sets.
        
        In particular, it generates
        
        - major set : Containing all matched sources (:func:`~eroML.utils.datasets.major_set`)
        - random set : Shift all source by a random amount and match  (:func:`~eroML.utils.datasets.random_set`)
        - training set : Best matching sources  (:func:`~eroML.utils.datasets.training_set`)
        - training+random : training set plus random source fullfilling the same criteria as the training set sources  (:func:`~eroML.utils.datasets.training_random_set`)
        
        Parameters
        ----------
        ero/gaia_fn : str
            Filenames for eROSITA and Gaia source lists
      
        """
        major_set(self.ero_fn, self.gaia_fn, self.major_fn)
        random_set(self.ero_fn, self.gaia_fn, self.random_fn)
        training_set(self.major_fn, self.training_fn, abs_dist_cutoff=training_abs_dist_cutoff, rel_dist_cutoff=training_rel_dist_cutoff)
        training_random_set(self.training_fn, self.random_fn, self.training_random_fn, abs_dist_cutoff=training_abs_dist_cutoff, rel_dist_cutoff=training_rel_dist_cutoff)
        
        return
        self.e = from_fits(ero_fn, mapper={"detUID":"srcID", "DEC":"Dec"})#, maxN=100)
        g = from_fits(gaia_fn, mapper={"source_id":"srcID", "ra":"RA", "dec":"Dec"})#, maxN=10)

        major_merged_catalog(ero_fn, gaia_fn)
        ff = to_fits(self.e)
        N_catalog = NN_distribution(ff, verbose=10)
        self.N_expected = N_catalog
        
