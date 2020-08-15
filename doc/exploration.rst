Exploration
============

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
