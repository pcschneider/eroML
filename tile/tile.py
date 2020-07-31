from ..utils.datasets import major_set, training_set, random_set
from ..utils.gaia_tools import gaia4ero, prepare_gaia
from ..utils.enrich import *
from ..utils.iso_tools import *
from ..ensemble import from_fits

class Tile():
    """
    Orchestrates the analyzes of a sky region.
    
    Therer are four main objectives:
      1. Create the master catalog, which contains all the information needed to identify the stars
      2. Create a training catalog, based on well positionally matched sources
      3. Create a random catalog to assess to reliability of the sources
      4. Prepare the data for ML, i.e., convert/normalize columns etc.
    """
    def __init__(self):
        """
        xxx
        """
        self.e = None
    def populate_filenames(self, ero_fn, fn_prefix="Tile"):
        """
        Generate the filenames to be used later
        """
        self.ero_fn = ero_fn
        self.gaia_fn = fn_prefix+"_gaia_sources.fits"
        self.major_fn = fn_prefix+"_major.fits"
        self.training_fn = fn_prefix+"_training.fits"
        self.random_fn = fn_prefix+"_random.fits"
        
    def prepare_data(self, ero_fn=None, gaia_fn=None):
        """
        """
        #return
        if ero_fn is None: ero_fn = self.ero_fn 
        if gaia_fn is None: gaia_fn = self.gaia_fn     
        
        # First, download and enrich Gaia data
        ####gaia4ero(ero_fn, ofn=gaia_fn, overwrite=True)
        prepare_gaia(ero_fn, gaia_fn, verbose=1)
        enrich_Gaia(gaia_fn)
        #g = from_fits(gaia_fn)
        add_iso_column(gaia_fn, gaia_fn, overwrite=True)        
        eligible(gaia_fn)
        sky_density(gaia_fn)
        sky_density(gaia_fn, filter_prop=None, out_col="sky_density")
        #exit()
        
        #e = from_fits(gaia_fn)
        # Enrich eROSITA data
        #print("fn",ero_fn)
        enrich_eROSITA(ero_fn, mapper={"detUID":"srcID", "DEC":"Dec"})
        
        
        
    def generate_sets(self, ero_fn=None, gaia_fn=None):       
        """
        """
        major_set(self.ero_fn, self.gaia_fn, self.major_fn)
        random_set(self.ero_fn, self.gaia_fn, self.random_fn)
        training_set(self.major_fn, self.training_fn)
        
        return
        self.e = from_fits(ero_fn, mapper={"detUID":"srcID", "DEC":"Dec"})#, maxN=100)
        g = from_fits(gaia_fn, mapper={"source_id":"srcID", "ra":"RA", "dec":"Dec"})#, maxN=10)

        major_merged_catalog(ero_fn, gaia_fn)
        ff = to_fits(self.e)
        N_catalog = NN_distribution(ff, verbose=10)
        self.N_expected = N_catalog
        
