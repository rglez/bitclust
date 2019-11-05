.. bitclust documentation master file, created by
   sphinx-quickstart on Mon Jan 28 18:21:15 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BitClust's documentation
===================================

Description
-----------

**BitClust** is a Python command-line interface (CLI) conceived for fast
clustering of relatively long Molecular Dynamics trajectories following the
Daura's algorithm [1]_. Retrieved clusters are roughly equivalent to those
reported by **VMD's** internal command **measure cluster** but they are computed
in a much faster way (see benchmark section for more details).

What **BitClust** offers is a classical tradeoff; RAM for speed. It can
calculate all pairwise distances between frames to run a clustering job and
then store them in memory instead of recalculating them whenever a cluster is found.

It is worth noting that memory resources have been deeply optimized by encoding similarity distances
as bits (0 if the distance is less equal than a specified threshold, 1 otherwise).
This encoding result in a storage reduction of at least 32X when compared to similar
algorithms that save the same information as single-precision float values.


Main Dependencies
-----------------

**BitClust** is built on the shoulders of two giants:

 *  `MDTraj software <http://mdtraj.org/1.9.0/>`_  that allows a very fast
    calculation of RMSD pairwise distances between all frames of trajectories in
    a parallelized fashion **and**

 * `bitarray third-party python library <https://pypi.org/project/bitarray/>`_ 
   which offers a memory-efficient data structure of bit-vectors (bit arrays)
   and a set of bitwise operations that are the very heart of our clustering
   implementation.


Citation
--------
If you make use of **BitClust** in your scientific work, **BitCool** and `cite it ;) <https://doi.org/10.1021/acs.jcim.9b00828>`_


Licence
-------
**BitClust** is licensed under GNU General Public License v3.0.

.. toctree::
   :maxdepth: 2
   :hidden:

   sections/installation
   sections/usage
   sections/help
   sections/benchmark
   sections/changelog    

References
----------
.. [1] Daura, X.; van Gunsteren, W. F.; Jaun, B.; Mark, A. E.; Gademann, K.; Seebach, D. Peptide Folding: When Simulation Meets Experiment. Angew. Chemie Int. Ed. 1999, 38 (1/2), 236–240.

