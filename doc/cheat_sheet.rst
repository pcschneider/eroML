Useful Commands
=================


   
git
---

To get the latest git version in the `svm` directory, just run::

  git clone stgh325@pc13:/home/hslxrsrv3/stgh325/hdd/Xdata/ero/svm svm

The corresponding path should be then added to PYTHONPATH::

  export PYTHONPATH=/home/majestix/hdd/Xdata/ero/
  export PYTHONPATH=/hs/pc13/data/stgh325/Xdata/ero
  

python
------

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
  
  
  
