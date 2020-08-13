.. eroML documentation master file, created by
   sphinx-quickstart on Fri Jul  3 16:29:59 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to eroML's documentation!
=================================
  
Run the stuff::

    python3.7 -c 'from eroML.utils import NN_distribution ; print(NN_distribution("merged.fits", ofn="x"))'
    python3.7 -c 'from eroML.utils.gaia_tools import gaia4ero ; gaia4ero("../eFEDS/SrcCat_V2T.fits", ofn="gaia.fits")'
  
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   learning.rst
   submodules.rst
   cheat_sheet.rst
      

Description of individual processing Steps
--------------------------------------------

1. Data preparation
    a. eROSITA data 
         - Annotate healpix, add FX
         - Split into individual tiles (=files)  
    b. Gaia data
         - Download Gaia (takes a long time)
         - Annotate Gaia (eligible, etc.)
         
2. Generate random eROSITA datasets        

3. (Position) Matching
    a. True eROSITA sources
    b. Random sources (N times)
    c. Enrich merged dataset(s)
    
4. Generate datasets
    a. Training
    b. Validation
    
5. Learn 

6. Match
              
      
Work Logic
-----------

Loop through Tiles 

  0. (method: :func:`~eroML.tile.tile.loop`)

  1. For each Tile: (method: :func:`~eroML.tile.tile.Tile.prepare_data`)
      a) Get Gaia sources
          - Get sky extent 
          - Download Gaia sources from archive
          - Convert Gaia data to fits-file
      b)  Prepare data (method: :func:
          - For Gaia, add columns: `Fg`, `iso_compatible`, `eligible`, `sky_density`, `sky_density_eligible`
          - For eROSITA, add columns: `Fx`       
          
  2. Generate data sets  (method: :func:`~eroML.tile.tile.Tile.generate_sets`)
      a) major set : Containing all matched sources (:func:`~eroML.utils.datasets.major_set`)
      b) random set : Shift all source by a random amount and match  (:func:`~eroML.utils.datasets.random_set`)
      c) training set : Best matching sources  (:func:`~eroML.utils.datasets.training_set`)
      d) training+random : training set plus random source fullfilling the same criteria as the training set sources (:func:`~eroML.utils.datasets.training_random_set`)
      
  3. Merge tiles (method: :func:``)
   
Entry points:


Todo
----

- Add quality criterium for NN_estimate
- Pixelize and merge


Testing
-------

Result of doctest:

.. doctest::
  :pyversion: > 3.8

Follow  `numpy/scipy <http://scipy.github.io/devdocs/dev/contributor/development_workflow.html>`_  rules.
  
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
