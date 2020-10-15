How to
======

Potentially useful entry points
--------------------------------

Once the Gaia data has been downloaded, one may

  0) re-run the enriching of the catalogs
  1) re-match the sources, i.e., produce new `master`, `training`, and `random`-sets

eFEDS example
-------------------

Required input Data
~~~~~~~~~~~~~~~~~~~~

The eFEDS source file::

  ../eFEDS/SrcCat_V2T.fits
  
The Gaia source tiles (expected number of files: 6305)::

  ../ero_data/Gaia_ID2_nside32_*.fits

How to download the Gaia tiles is described in 
:ref:`How to download Gaia tiles<How to download Gaia tiles>`.  
  
Generate required files
~~~~~~~~~~~~~~~~~~~~~~~~

We need to transform the eROSITA into the merged file containing all relevant
source properties. This can be done by setting::

  p37 run_eroML.py eFEDS_all.ini
  
This splits the eFEDS source file into individual source files according to 
their healpix index, matches the X-ray sources againt the Gaia catalog, 
generates random sources, and merges all into::

  merged_random_eFEDS.fits
  merged_training_eFEDS.fits
  merged_major_eFEDS.fits
  
These must be adapted for SVM by running.
