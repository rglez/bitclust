Usage Examples
==============


#usage: bitclust-cli.py [-h] [-top TOPOLOGY] [-traj TRAJECTORY] [-first FIRST]
#                       [-last LAST] [-stride STRIDE] [-sel SELECTION]
#                       [-cutoff CUTOFF] [-minsize MINSIZE] [-ref REFERENCE]
#                       [-odir OUTDIR]
#
#BitClust: Fast & memory efficient clustering of MD trajectories
#
#
#optional arguments:
#  -h, --help        show this help message and exit
#  -top TOPOLOGY     path to topology file (psf/pdb)
#  -traj TRAJECTORY  path to trajectory file
#  -first FIRST      first frame to analyze (starting from 0)
#  -last LAST        last frame to analyze (starting from 0)
#  -stride STRIDE    stride of frames to analyze
#  -sel SELECTION    atom selection (MDTraj syntax)
#  -cutoff CUTOFF    RMSD cutoff for pairwise comparisons in A
#  -minsize MINSIZE  minimum number of frames of returned cluster
#  -ref REFERENCE    reference frame to align trajectory
#  -odir OUTDIR      output directory to store analysis

