Benchmark
=========

A set of commonly used options for clustering MD simulations has been chosen for
performance comparison between BitClust and other clustering software. Run time
and memory consumption of each method is reported in Table 2 for the three trajectories
6K, 100K, and 500K described in the Computational Details section (having 6000, 100000 and
500000 frames respectively).

The clustering method selected for each software was as follows; **Daura** for BitClust and GROMACS (through the gromos option),
**quality threshold** for py-MS and VMD, **qt-like** for Wordom and **median-linkage** for TTClust. We would like to stress out
that despite their native denomination in their original software, the chosen algorithms correspond all to **Daura**, as it has been
recently reported [1]_

In the case of VMD, we decided to show the performance of processing five (VMD-5, the default value) and all (VMD-ALL) clusters
as a way to evaluate the usefulness of its implementation, which is specially conceived for preserving memory resources.


**Table 2**  Run time ([hh:]mm:ss) and memory peaks (GB) comparison of several clustering algorithms.
A crash was declared if the job consumed more than 64 GB of RAM or last for more than 60 h.
(For an extensive discussion of this benchmark, you can refer to the academic article [2]_)


+----------+----------+-------------+----------+-------------+----------+-------------+
| Software |      Trajectory 6K     |     Trajectory 100K    |     Trajectory 500K    |
+==========+==========+=============+==========+=============+==========+=============+
|          | Run time | Memory peak | Run time | Memory peak | Run time | Memory peak |
+----------+----------+-------------+----------+-------------+----------+-------------+
| BitClust | 00:04    |    0.15     | 01:25:28 |     9.41    | 06:00:08 |    33.84    |
+----------+----------+-------------+----------+-------------+----------+-------------+
| py-MS    | 00:04    |    0.15     | 01:46:02 |    61.54    | 03:03:49 |    CRASH    |
+----------+----------+-------------+----------+-------------+----------+-------------+
| Wordom   | 00:04    |    0.15     | 01:26:13 |    CRASH    | 02:54:23 |    CRASH    |
+----------+----------+-------------+----------+-------------+----------+-------------+
| VMD-5    | 00:04    |    0.15     | 04:59:22 |     6.97    | 29:34:43 |     4.13    |
+----------+----------+-------------+----------+-------------+----------+-------------+
| VMD-ALL  | 00:04    |    0.15     | 05:08:10 |     6.97    |  CRASH   |     4.14    |
+----------+----------+-------------+----------+-------------+----------+-------------+
| TTClust  | 00:04    |    0.15     | 00:01:59 |    CRASH    | 00:00:37 |    CRASH    |
+----------+----------+-------------+----------+-------------+----------+-------------+
| GROMACS  | 00:04    |    0.15     | 26:13:29 |    52.22    | 01:09:03 |    CRASH    |
+----------+----------+-------------+----------+-------------+----------+-------------+





References
----------
.. [1] Roy González-Alemán, David Hernández-Castillo, Julio Caballero, and Luis A. Montero-Cabrera. Journal of Chemical Information and Modeling. DOI: 10.1021/acs.jcim.9b00558 

.. [2] Roy González-Alemán, David Hernández-Castillo, Alejandro Rodríguez-Serradet, Julio Caballero, Erix W. Hernández-Rodríguez, and Luis Alberto Montero-Cabrera. Journal of Chemical Information and Modeling. DOI: 10.1021/acs.jcim.9b00828.


