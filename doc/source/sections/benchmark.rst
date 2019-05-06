Benchmark
=========

We have tested four different situations to demonstrate the superiority of
**BitClust** in terms of speed. The first of them refers to **BitClust**
processing all the clusters (denoted as BC-ALL). The other ones refer to VMD
and three possibilities: to process only one cluster, to process five clusters,
and to process all clusters of the trajectory (denoted as VMD-1, VMD-5, and
VMD-ALL respectively).

The only trajectory which permitted to evaluate all four possibilities in a
rational amount of time was 6K (Figure A). As it can be appreciated, the necessary time to
retrieve all clusters by **BitClust** is similar to time taken by VMD to retrieve
just the first cluster. This trend also holds for the bigger trajectories.

.. admonition :: **BitClust** (time & memory) benchmark
   
  **Figure A:** Elapsed time for retrieving 6K clusters: all of them using
  BitClust (BC-ALL) and one, five and all using VMD (VMD-1, VMD-5 and VMD-ALL
  respectively)

  **Figure B:** Runtime comparison for bigger trajectories. VMD would have taken
  an excessive amount of time to process all clusters and consequently this case
  was not measured. Instead, we calculated the percent of clustered frames by VMD
  after 24 hours to have an approximate idea of VMD behavior.

  **Figure C:** Time decomposition for two main processes associated to BitClust;
  the calculation of all vs. all RMSD and the clustering part. Most of the time
  is spent in the first task.

  **Figure D:** Memory decomposition for two main objects BitClust loads in RAM;
  the trajectory and the RMSD matrix.

.. figure :: /custom/benchmark.png
   :align: center
