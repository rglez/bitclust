Help
====

Basic help
----------
**BitClust** help is displayed in the console when typing **bitclust -h** ::

  $ bitclust -h 

  usage: bitclust [-h] [-top TOPOLOGY] [-traj TRAJECTORY] [-first FIRST]
                  [-last LAST] [-stride STRIDE] [-sel SELECTION]
                  [-cutoff CUTOFF] [-minsize MINSIZE] [-ref REFERENCE]
                  [-odir OUTDIR]

  BitClust: Fast & memory efficient clustering of long MD trajectories

  optional arguments:
    -h, --help        show this help message and exit
    -top TOPOLOGY     path to topology file (psf/pdb)
    -traj TRAJECTORY  path to trajectory file
    -first FIRST      first frame to analyze (starting from 0)
    -last LAST        last frame to analyze (starting from 0)
    -stride STRIDE    stride of frames to analyze
    -sel SELECTION    atom selection (MDTraj syntax)
    -cutoff CUTOFF    RMSD cutoff for pairwise comparisons in A
    -minsize MINSIZE  minimum number of frames inside returned clusters
    -ref REFERENCE    reference frame to align trajectory
    -odir OUTDIR      output directory to store analysis



Arguments details
-----------------

``-traj (str):`` This is the only argument that is **always** required. Valid
extensions for trajectories are ``.dcd``, ``.dtr``, ``.hdf5``, ``.xyz``, ``.binpos``,
``.netcdf``, ``.prmtop``, ``.lh5``, ``.pdb``, ``.trr``, ``.xtc``, ``.xml``,
``.arc``, ``.lammpstrj`` and ``.hoomdxml``.

``-top (str):`` If trajectory format includes topological information this
argument is not required. Otherwise, user must pass a path to a topology
file. Valid topology extensions are  ``.pdb``, ``.pdb.gz``,
``.h5``, ``.lh5``, ``.prmtop``, ``.parm7``, ``.prm7``, ``.psf``, ``.mol2``,
``.hoomdxml``, ``.gro``, ``.arc`` and ``.hdf5``.

``-first (int, default=0):`` First frame to analyze (starting count from 0)

``-last (int, default=-1):`` Last frame to analyze (starting count from 0). Value -1
indicates that last frame will be used.

``-stride (int, default=1):`` Stride of frames to analyze. Use this argument to
reduce trajectory size when performing exploratory analysis of cutoff value.

``-sel (str, default='all'):`` Atom selection. **BitClust** inherits ``MDtraj``
syntax selection which is very flexible. For a deeper insight please refer
to the ``MDTraj atom selection reference`` original documentation. Common cases
are listed at usage examples section. 
   
``-cutoff (int, default=1):`` RMSD cutoff for similarity measures given in Angstroms
(1 A = 0.1 nm).

``-minsize (int, default=2):`` Minimum number of frames inside returned clusters.
0 is not a meaningful value and 1 implies an unclustered frame (no other frame is
similar to it). Greater values of this parameter can speed up the algorithm.

``-ref (int, default=0):`` Reference frame to align trajectory.

``-odir (str, default="."):`` Output directory to store analysis. If not specified,
files and figs will be stored in current working directory.

