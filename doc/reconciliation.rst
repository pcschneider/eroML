Reconciliation
===============

Current input files
-------------------

The input file is::

  major_eFEDS_classified.fits
  
with the auxiliary files::

  ../ero_data/efeds_c001_V3_main_HamStar_internal2.fits
  ../../ero_data/nway.fits
  
HamStar-formatted files
------------------------------

Generate files akin to the HamStar format::

  p37 tools/gen_HamStar_file.py # -> eFEDS_HamStar.fits
  p37 tools/update_Sebs_HamStar_file.py # -> ../ero_data/efeds_c001_V3_main_HamStar_internal2.fits
  p37 tools/update_nway_HamStar_file.py # -> ../ero_data/nway2.fits
  
Master table
--------------

Generate master table (attention: filenames are hardcoded)::

  p37 reconciliation/gen_master_table.py # -> ero_master.fits
  fits2sqlite ero_master.fits ero_master.db
  
Comparison plot
--------------

Run::

  p37 reconciliation/compare_all.py
