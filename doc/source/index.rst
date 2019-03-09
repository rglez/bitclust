.. bitclust documentation master file, created by
   sphinx-quickstart on Mon Jan 28 18:21:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

   
BitClust
========

Nowadays very long simulations are carried on routinely. Enhanced sampling
methods like metadynamics, REMD and accelerated dynamics allow scaping from
potential energy minima, returning trajectories that are conformationally sparsed
and where every cluster can be potentially important to detect and analyze.
Timescales in the order of micro and milliseconds are becoming familiar and
consequently, tools for analyzing such huge trajectories must be upgraded.

**BitClust** is a Python command line interface (CLI) conceived for fast
clustering of relatively long Molecular Dynamics trajectories following
Daura's algorithm [cite]. Retrieved clusters are roughly equivalent to those
obtained by VMD's internal command *measure cluster* but in a much faster way
(see benchmark section for details).

**BitClust** is built on the shoulders of two giants:

 *  `MDTraj <http://mdtraj.org/1.9.0/>`_ software that allows a very fast
    calculation of RMSD pairwise distances between all frames of trajectories in
    a parallelized fashion **and**

 * `bitarray <https://pypi.org/project/bitarray/>`_ 3rd party python library
   which offers a memory efficient data structure of bit vectors (bit arrays)
   and a set of bitwise operations that are the very heart of the clustering
   implementation it offers.


What **BitClust** offers is a classical tradeoff; RAM for speed. It is able to
calculate all the information needed to run a clustering job and then store it
in memory instead of recalculating it at every iteration. It is worth noting
that used memory have been deeply optimized by encoding similarity distances
as bits. This encoding result in a storage reduction as high as 16X compared to
algorithms that saves pairwise distances as single precision float values.


Citation
--------
If you make use of **BitClust** in your scientific work, please contribute with
a citation. The BibTeX reference is:


Licence
-------
**BitClust** is licensed under GNU General Public License v3.0.


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   sections/installation
   sections/arguments
   sections/usage_examples
   sections/benchmark
   sections/code_reference
   sections/changelog


Indices
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
