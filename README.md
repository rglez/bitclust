# BitClust: Fast and memory efficient clustering of long Molecular Dynamics


# Home Page
-----------

BitClust´s latest documentation is available [here](https://bitclust.readthedocs.io/en/latest/) 


# Description
-------------

**BitClust** is a Python command line interface (CLI) conceived for fast
clustering of relatively long Molecular Dynamics trajectories following
Daura's algorithm [1]. Retrieved clusters are roughly equivalent to those
reported by **VMD's** internal command **measure cluster** but they are computed in a
much faster way (see benchmark section for more details).


# Motivation

Nowadays very long simulations are carried on routinely. Enhanced sampling
methods like metadynamics, REMD and accelerated dynamics allow escaping from
potential energy minima, returning trajectories that are conformationally sparsed
and where every cluster can be potentially important to detect and analyze. Improvements
on software designed to address this task is an important field of research.

**BitClust** offer is a classical tradeoff; RAM for speed. It is able to
calculate all pairwise distances between frames to run a clustering job and
then store them in memory instead of recalculating them whenever a cluster is found.
It is worth noting that used memory has been deeply optimized by encoding similarity distances
as bits (0 if the distance is less equal than a specified threshold, 1 otherwise).
This encoding result in a storage reduction as high as 16X compared to similar algorithms
that saves the same information as single precision float values.


# Main Dependencies

**BitClust** is built on the shoulders of two giants:

 *  [MDTraj software](http://mdtraj.org/1.9.0/)  that allows a very fast
    calculation of RMSD pairwise distances between all frames of trajectories in
    a parallelized fashion **and**

 * [bitarray third-party python library](https://pypi.org/project/bitarray/) 
   which offers a memory efficient data structure of bit vectors (bit arrays)
   and a set of bitwise operations that are the very heart of our clustering
   implementation.


# Citation

If you make use of **BitClust** in your scientific work, **BeCool** and cite it ;)

The BibTeX reference is:

.. todo::
  insert once published.


# Licence

**BitClust** is licensed under GNU General Public License v3.0.
  
  
# Reference

[1] Daura, X.; van Gunsteren, W. F.; Jaun, B.; Mark, A. E.; Gademann, K.; Seebach, D. Peptide Folding: When Simulation Meets Experiment. Angew. Chemie Int. Ed. 1999, 38 (1/2), 236–240.

