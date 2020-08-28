.. eroML documentation master file, created by
   sphinx-quickstart on Fri Jul  3 16:29:59 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to eroML's documentation!
=================================
  
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   exploration.rst
   training.rst 
   learning.rst
   submodules.rst
   cheat_sheet.rst
      

              
      
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

   
.. graphviz::

   digraph Flatland {
   
      a -> b -> c -> g; 
      a  [shape=polygon,sides=4]
      b  [shape=polygon,sides=5]
      c  [shape=polygon,sides=6]
   
      g [peripheries=3,color=yellow];
      s [shape=invtriangle,peripheries=1,color=red,style=filled];
      w  [shape=triangle,peripheries=1,color=blue,style=filled];
      
      }
   
   

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
