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
      
   
Work Logic
-----------

Loop through Tiles 

  1. For each Tile: (method: :func:`~eroML.utils.gaia_tools.gaia4ero`)
      - Get sky extent 
      - Download Gaia sources from archive
      - Convert Gaia data to fits-file
  2. Prepare data 
   


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
