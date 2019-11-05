Benchmark
=========

(For an extensive discussion of this benchmark, you can refer to the academic article [1]_)


A set of commonly used options for clustering MD simulations has been chosen for
performance comparison between BitClust and other clustering software. Run time
and memory consumption of each method is reported in Table 2 for the three trajectories
6K, 100K, and 500K described in the Computational Details section (having 6000, 100000 and
500000 frames respectively).

The clustering method selected for each software was as follows; **Daura** for BitClust and GROMACS (through the gromos option),
**quality threshold** for py-MS and VMD, **qt-like** for Wordom and **median-linkage** for TTClust. We would like to stress out
that despite their native denomination in their original software, the chosen algorithms correspond all to **Daura**, as it has been
recently reported [2]_

In the case of VMD, we decided to show the performance of processing five (VMD-5, the default value) and all (VMD-ALL) clusters
as a way to evaluate the usefulness of its implementation, which is specially conceived for preserving memory resources.


**Table 2**  Run time and memory performance comparison of several clustering algorithms.

==========  ================ ================  =================== ================   =================== ================
                        Trajectory 6K                         Trajectory 100K                         Trajectory 500K
Software    ---------------- ----------------  ------------------- ----------------   ------------------- ----------------
            Run time (mm:ss) Memory peak (GB)  Run time (hh:mm:ss) Memory peak (GB)   Run time (hh:mm:ss) Memory peak (GB)
==========  ================ ================  =================== ================   =================== ================
BitClust        00:04              0.15            01:25:28              9.41             06:00:08              33.84
py-MS           00:11              0.41            01:46:02             61.54           **03:03:49**          **crash**
Wordom          01:52              0.10          **01:26:13**         **crash**         **02:54:23**          **crash**
VMD-5           00:14              0.09            04:59:22              6.97             29:34:43               4.13
VMD-ALL         02:05              0.10            05:08:10              6.97           **60:00:00**           **4.14**
TTClust         01:30              1.16          **00:01:59**         **crash**           00:00:37            **crash**
GROMACS         02:10              0.16            26:13:29             52.22           **01:09:03**          **crash**
==========  ================ ================  =================== ================   =================== ================


References
----------
.. [1] Roy González-Alemán, David Hernández-Castillo, Alejandro Rodríguez-Serradet, Julio Caballero, Erix W. Hernández-Rodríguez, and Luis Alberto Montero-Cabrera. Journal of Chemical Information and Modeling. DOI: 10.1021/acs.jcim.9b00828.

.. [2] Roy González-Alemán, David Hernández-Castillo, Julio Caballero, and Luis A. Montero-Cabrera. Journal of Chemical Information and Modeling. DOI: 10.1021/acs.jcim.9b00558 

