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

      
How to download Gaia tiles
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Gaia sources that need to be downloaded depend on the populated hpix (and 
NSIDE). It is therefore recommended to first calculate the hpix for a given 
source file::

  p37 run_eroML.py example.ini --steps Healpix/calculate

Then, the following command downloads the required Gaia tiles::  

  p37 run_eroML.py example.ini --steps Gaia_Download/perform
  
One can also combine both steps::

  p37 run_eroML.py example.ini --steps Healpix/calculate Gaia_Download/perform
  
Two parameters control how eroML performs the download:

  1) ``Gaia_Download/check_alternate``
      If this option is selected, Gaia-files with another runID can be used, 
      i.e., with ``Gaia_Download/directory`` ``/Gaia_*_nsideXX_YYY.fits`` where
      XX and YYY are  HPIX/nside and hpix-index, respectively. If such a file
      exists, a copy will be created.
      
      It may be advisable to run ``Enrich_Gaia/perform`` afterwards to ensure 
      that the correct filter-criteria are used.
      
  2) ``Gaia_Download/overwrite``
      Existing Gaia-files may be overwritten. If they already exist, no 
      download is initiated.
  
The most relevant entries in the ini-file that should point to the correct 
files are:

  - ``Sources/X_filename``
  - ``Sources/X_filename_hp`` 
  - ``General/runID``
  - ``Gaia_Download/directory``
  - ``Gaia_Download/prefix``

.. Relevant content of data sets
.. ------------------------------
.. 
.. Each data set has its specific, relevant columns:
.. 
..   - eROSITA source list (*Sources:ero_filename*)[``ero_filename``]



Re-run the enriching of the catalogs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To apply new criteria for "eligible" source, create an ``ini`` file with the 
following keywords changed wrt to the original file::

  [Healpix]
  calculate=False

  [Gaia_Download]
  perform=False
  
Alternatively, one can perform only a single processing step, e.g.,::

  p37 run_eroML.py eFEDS_EDR3.ini --steps Xdata_preparation/enrich
  p37 run_eroML.py eFEDS_EDR3.ini --steps Enrich_Gaia/perform

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

The Training sample
~~~~~~~~~~~~~~~~~~~

Construction of the training sample is a multi-step process. 

1) Estimate catalog fraction. This is done by running::
    
    p37 tools/estimate_catalog_N.py
  
2) Generate a training sample for the geometric classifier via::

    p37 positions/simu4major.py 2060 --conf eFEDS_EDR3.ini --ofn offs2.dat -o --rnd_factor=12.6

3) Train the SVM classifier based on the generated positions::

    p37 positions/generate_classifier.py 

4) Classify real associations::

    p37 positions/select_geometric_training_objects.py

5) Astrophysical training sample::

    p37 classify/gen_training_sample.py
    
6) Astrophysical screening (empirical Lx/Lbol screening, absolute Lx screening)::

    p37 classify/prepare.py
    
7) Train classifier::
    
    p37 classify/learn.py ; p37 classify/predict_n_check.py 

8) Identify stars::
    
    p37 classify/write_stars.py

9) Write HamStar-liek file::
    
    p37 tools/gen_HamStar_file.py
    
9) Compare with Sebastian::
    
    p37 classify/cmp.py
    

A new training sample can be constructed by running::

  p37 run_eroML.py eFEDS_tmp.ini  
  cd classify
  p37 clean_training.py 
  p37 learn.py
  p37 predict_n_check.py
  
  
