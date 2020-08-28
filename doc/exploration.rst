Data Description and Data Sets
===============================


Individual data sets
--------------------
Names in brackets pertain to their corresponding name in the config-file (**idx** is the healpix index):
Names in square brackets give the identifier for `file4`.

  - eROSITA source list (*Sources:ero_filename*)[``ero_filename``]
      This file is provided by the catalog team and will not be changed
      
  - eROSITA source list with healpix indices (*Sources:ero_filename_hp*)[``ero_filename_hp``]
      
  - eROSITA tiles (*data_dir* / *prefix* _nside *nside* _ **idx**.fits)[``ero_tiles``]
      This file may be annotated, but shall always contain the full content. 
  - Gaia tiles (*data_dir*  / *prefix* _nside *nside* _**idx**.fits)[``gaia_tiles``]    
      This file may be annotated, but shall always contain the full content. 
  - For each tile, the following sets
      a) *major* [``major_tiles``]
      b) *random* [``random_tiles``]
      c) *training* [``training_tiles``]
    Plus the same data sets in the _small_ incarnation [``*_small``]  
  - Finally, the merged data sets
      a) *major* [``major``]
      b) *random* [``random``]
      c) *training* [``training``]
      
      
Data Processing Steps
--------------------------------------------

1. Data preparation
    a. eROSITA data 
         - Annotate healpix, add FX
         - Split into individual tiles (=files)  
    b. Gaia data
         - Download Gaia (takes a long time)
         - Annotate Gaia (eligible, etc.); this step is required to work on the full data 
         
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

Use a sub-sample
----------------------

Once, each eROSITA source has an asigned healpix (`ifn`), one can generate random hpix list::

  from eroML.tile import populated_hpix
  import numpy as np
  
  # The filename containing the healpix-annotated source list
  ifn = "../ero_data/eRASS1_hp.fits"
  ofn = "hpix_list.dat"
  N = 100
  
  hpix = populated_hpix(ifn)
  out = hpix[np.random.choice(range(len(hpix)), size=N, replace=False)]
  print(min(out), max(out))
  np.savetxt(ofn, out, fmt="%i")
  
Now, one can add the following::

  [Healpix]
  pix_file=hpix_list.dat
  
to the config (`.ini`) and only those pixels will be used.


Sky Density
------------
The sky density can be displayed by running::
  
  p37 tools/sky_density.py
