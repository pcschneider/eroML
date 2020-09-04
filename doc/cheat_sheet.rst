Useful Commands
=================


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



   
git
---

To get the latest git version in the `svm` directory, just run::

  git clone stgh325@pc13:/home/hslxrsrv3/stgh325/hdd/Xdata/ero/svm svm

The corresponding path should be then added to PYTHONPATH::

  export PYTHONPATH=/home/majestix/hdd/Xdata/ero/
  export PYTHONPATH=/hs/pc13/data/stgh325/Xdata/ero
  

Count lines::

   git ls-files | xargs wc -l

  
python
------


Run the stuff::

    python3.7 -c 'from eroML.utils import NN_distribution ; print(NN_distribution("merged.fits", ofn="x"))'
    python3.7 -c 'from eroML.utils.gaia_tools import gaia4ero ; gaia4ero("../eFEDS/SrcCat_V2T.fits", ofn="gaia.fits")'
  

Testing::

  p37 -c "import eroML ; eroML.test()"

or on pc13::

  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hs/pc13/data/stgh325/python/lib/
  python3.7 -c "import eroML ; eroML.test()"
  
  
::

    from importlib import reload
    reload(vo)

    
::

  echo "__pycache__" >> .gitignore 

sphinx
------

Install via pip, and soft-link sphinx-* etc into ``~/bin`` and then run
::

  sphinx-quickstart
  
  
  
