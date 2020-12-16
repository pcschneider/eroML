How to
======

Usage
-----
The software is controlled by ini-files, individual steps can be also selected
on the command line::

  p37 run_eroML.py -h



Default values are provided in::

  default.ini
  
  

Potentially useful entry points
--------------------------------

There are a few steps, that one may want change after the first run. 

Re-run the enriching of the catalogs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To apply new criteria for "eligible" source, create an ``ini`` file with the 
following keywords changed wrt to the original file::

  [Healpix]
  calculate=False

  [Gaia_Download]
  perform=False
  
Alternatively, one can perform only a single step, e.g.,::

  p37 run_eroML.py eFEDS_EDR3.ini --steps Xdata_preparation/enrich
  
However, this still requires to run all subsequent steps, too! Something like::

  p37 run_eroML.py eFEDS_EDR3.ini --steps Xdata_preparation/enrich Datasets/major Merging/shrink Merging/major
  
would create a new, merged major dataset for eFEDS, in this case::

  ../ero_data/merged_eFEDS_EDR3.fits
  
Re-match catalog
~~~~~~~~~~~~~~~~
To produce new `master`, `training`, and `random`-sets, run



eFEDS example
-------------------

As an example, we describe the processing of the eFEDS.

Required input Data
~~~~~~~~~~~~~~~~~~~~

The eFEDS source file::

  ../ero_data/eFEDS_c001_hs.fits
  
The Gaia source tiles (expected number of files: 6305)::

  ../ero_data/Gaia_ID2_nside32_*.fits

How to download the Gaia tiles is described in 
:ref:`How to download Gaia tiles<How to download Gaia tiles>`.  
  
Generate required files
~~~~~~~~~~~~~~~~~~~~~~~~

We need to transform the eROSITA files into the merged file containing all 
relevant source properties. This can be done by setting::

  p37 run_eroML.py eFEDS_all.ini
  
This splits the eFEDS source file into individual source files according to 
their healpix index, matches the X-ray sources againt the Gaia catalog, 
generates random sources, and merges all into::

  merged_random_eFEDS.fits
  merged_training_eFEDS.fits
  merged_major_eFEDS.fits
  
These must be adapted for SVM by running.

A new training sample can be constructed by running::

  p37 run_eroML.py eFEDS_tmp.ini  
  cd classify
  p37 clean_training.py 
  p37 learn.py
  p37 predict_n_check.py
  
  
  
