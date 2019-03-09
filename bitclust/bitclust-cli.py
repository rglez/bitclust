#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
.. note ::

  | **Author      :** Roy Gonzalez Aleman
  | **Contact     :** [roy_gonzalez@fq.uh.cu, roy.gonzalez.aleman@gmail.com]

**Description:**
  'BitClust: Fast & memory efficient clustering of MD trajectories'

'''
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
from subprocess import run
import matplotlib.pyplot as plt
from collections import OrderedDict

import mdtraj as md
from bitarray import bitarray as ba


# @profile
def parse_arguments():
    '''
    DESCRIPTION
    Parse all user arguments from the command line.

    Return:
        user_inputs (parser.argparse): namespace with user input arguments.
    '''

    # Initializing argparse ---------------------------------------------------
    desc = '\nBitClust: Fast & memory efficient clustering of MD trajectories'
    parser = argparse.ArgumentParser(description=desc,
                                     add_help=True,
                                     epilog='As simple as that ;)')
    # Arguments: loading trajectory -------------------------------------------
    parser.add_argument('-top', dest='topology', action='store',
                        help='path to topology file (psf/pdb)',
                        type=str, required=False)
    parser.add_argument('-traj', dest='trajectory', action='store',
                        help='path to trajectory file',
                        type=str)
    parser.add_argument('-first', dest='first',  action='store',
                        help='first frame to analyze (starting from 0)',
                        type=int, required=False, default=0)
    parser.add_argument('-last', dest='last', action='store',
                        help='last frame to analyze (starting from 0)',
                        type=int, required=False, default=-1)
    parser.add_argument('-stride', dest='stride', action='store',
                        help='stride of frames to analyze',
                        type=int, required=False, default=1)
    parser.add_argument('-sel', dest='selection', action='store',
                        help='atom selection (MDTraj syntax)',
                        type=str, required=False, default='all')
    # Arguments: clustering ---------------------------------------------------
    parser.add_argument('-cutoff', action='store', dest='cutoff',
                        help='RMSD cutoff for pairwise comparisons in A',
                        type=float, required=False, default=1.0)
    parser.add_argument('-minsize', action='store', dest='minsize',
                        help='minimum number of frames of returned cluster',
                        type=int, required=False, default=2)
    parser.add_argument('-ref', action='store', dest='reference',
                        help='reference frame to align trajectory',
                        type=int, required=False, default=0)
    # Arguments: analysis -----------------------------------------------------
    parser.add_argument('-odir', action='store', dest='outdir',
                        help='output directory to store analysis',
                        type=str, required=False, default='./')
    user_inputs = parser.parse_args()
    return user_inputs


# @profile
def load_trajectory(args):
    '''
    DESCRIPTION
    Loads trajectory file using MDTraj. If trajectory format is h5, lh5 or
    pdb, topology file is not required. Otherwise, you should specify a
    topology file.

    Arguments:
        args (argparse.Namespace): user input parameters parsed by argparse.
    Return:
        trajectory (mdtraj.Trajectory): trajectory object for further analysis.
    '''

    traj_file = args.trajectory
    traj_ext = traj_file.split('.')[-1]
    # Does trajectory file format need topology ? -----------------------------
    if traj_ext in ['h5', 'lh5', 'pdb']:
        trajectory = md.load(traj_file)
    else:
        trajectory = md.load(traj_file, top=args.topology)
    # Reduce RAM consumption by loading selected atoms only -------------------
    if args.selection != 'all':
        sel_indx = trajectory.topology.select(args.selection)
        sel_traj = trajectory.restrict_atoms(sel_indx)
        del trajectory
        trajectory = sel_traj[args.first:args.last:args.stride]
    else:
        trajectory = trajectory[args.first:args.last:args.stride]
    return trajectory


# @profile
def calc_rmsd_matrix(trajectory, args):
    '''
    DESCRIPTION
    Calculates RMSD square matrix using MDTraj. Pairwise similarity is saved in
    RAM as bits (dict of bitarrays), not floats.

    Args:
        trajectory (mdtraj.Trajectory): trajectory object to analyze.
        args (argparse.Namespace): user input parameters parsed by argparse.
    Return:
        matrix (collections.OrderedDict): dict of bitarrays
        degrees (collections.OrderedDict): dict of bitarrays´ lenghts.
    '''

    trajectory.center_coordinates()
    cutoff = np.float32(args.cutoff)/10  # transform A to nm (MDTraj coherence)
    matrix = OrderedDict()
    degrees = OrderedDict()
    to_explore = range(trajectory.n_frames)
    for i in to_explore:
        rmsd_ = md.rmsd(trajectory, trajectory, i, precentered=True)
#       closes = np.isclose(rmsd_, threshold, atol=1e-15, rtol=1e-15)
#       rounderrors = np.where(closes)[0]
        vector_np = np.less_equal(rmsd_, cutoff)
#       vector_np[rounderrors] = True
        matrix.update({i: ba(list(vector_np))})
        degrees.update({i: vector_np.sum()})
    return matrix, degrees


# @profile
def bitclusterize(matrix, degrees):
    '''
    DESCRIPTION
    Clusters bit matrix using bitwise operations.

    Args:
        matrix (collections.OrderedDict): dict of bitarrays
        degrees (collections.OrderedDict): dict of bitarrays lenghts.
    Return:
        clusters (numpy.ndarray): array of clusters ID.
        leaders (list) : list of clusters´ centers ID.
    '''
    # Memory allocation for clusters container --------------------------------
    clusters = np.empty(len(matrix), dtype='int32')
    clusters.fill(-1)
    # Declare all bits as available -------------------------------------------
    available_bits = ba(int(len(degrees)))
    available_bits.setall('1')
    # Start iterative switching -----------------------------------------------
    leaders = []
    clust_id = 0
    while True:
        # Find the biggest cluster --------------------------------------------
        leader = max(degrees.items(), key=lambda x: x[1])[0]
        biggest_cluster = matrix[leader] & available_bits
#        biggest_cluster = matrix[leader] # testing!!!
        available_bits = (available_bits ^ biggest_cluster) & available_bits
        leaders.append(leader)
        # Break 1: all candidates cluster have degree 1 (can´t clusterize) ----
        if (len(set(degrees.values())) == 1) and \
                (list(degrees.values())[0] == 1):
            return clusters, leaders
            break
        # Prune biggest cluster and its members from matrix -------------------
        discard = biggest_cluster.search(ba('1'))
        clusters[discard] = clust_id
        for element in discard:
            del degrees[element]
        # Assign next cluster ID ----------------------------------------------
        clust_id += 1
        # Update degrees of unclustered frames --------------------------------
        degrees_list = list(degrees.keys())
        for degree in degrees_list:
            degrees[degree] = ba.fast_hw_and(available_bits, matrix[degree])
        # Break 2: No more candidates available (empty matrix) ----------------
        if not degrees:
            return clusters, leaders
            break


# @profile
def calc_rms_vectors(trajectory, args):
    '''
    DESCRIPTION
    Calculates RMSD square matrix using MDTraj. Pairwise similarity is saved in
    RAM as bits (dict of bitarrays), not floats.

    Args:
        trajectory (mdtraj.Trajectory): trajectory object to analyze.
        args (argparse.Namespace): user input parameters parsed by argparse.
    Return:
        matrix (collections.OrderedDict): dict of bitarrays
        degrees (collections.OrderedDict): dict of bitarrays´ lenghts.
    '''

    trajectory.center_coordinates()
    reference = args.reference
    rmsd_ = md.rmsd(trajectory, trajectory, reference, precentered=True)
#    rmsf_ = md.rmsf(trajectory, trajectory, reference, precentered=True)
#    return rmsd_, rmsf_
    return rmsd_


def get_frames_stats(clusters, leaders):
    '''
    DESCRIPTION
    Gets plain text file "frames_statistics.txt" which contains as columns:
    frameID, clusterID and rmsd distance.

    Args:
        clusters
        leaders
        rmsd_
    Return:
        frames_df
    '''
    out_name = 'frames_statistics.txt'
    frames_df = pd.DataFrame(columns=['frame', 'cluster_id', 'rmsd'])
    frames_df['frame'] = range(len(clusters))
    frames_df['cluster_id'] = clusters
    frames_df['rmsd'] = rmsd_*10
    with open(out_name, 'wt') as on:
        frames_df.to_string(buf=on, index=False)
    return frames_df


def get_cluster_stats(clusters, leaders):
    '''
    DESCRIPTION
    Args:
    Return:
    '''
    out_name = 'cluster_statistics.txt'
    clusters_df = pd.DataFrame(columns=['cluster_id', 'size', 'leader',
                                        'percent'])
    clusters_df['cluster_id'] = range(-1, len(leaders)-1)
    leaders.insert(0, leaders[-1])
    leaders.pop()
    clusters_df['leader'] = leaders
    sizes = []

    if -1 in leaders:
        for x in range(-1, len(leaders)-1):
            sizes.append(len(np.where(clusters == x)[0]))
    else:
        for x in range(0, len(leaders)):
            sizes.append(len(np.where(clusters == x)[0]))

    clusters_df['size'] = sizes
    sum_ = sum(sizes)

    percents = [round(x/sum_*100, 4) for x in sizes]
    clusters_df['percent'] = percents
    with open(out_name, 'wt') as on:
        clusters_df.to_string(buf=on, index=False)
    return clusters_df


def generic_matplotlib():
    '''
    DESCRIPTION
    Some customizations of matplotlib.
    '''
    mpl.rc('figure', autolayout=True, figsize=[3.33, 2.5], dpi=300)
    mpl.rc('font', family='STIXGeneral')
    mpl.rc('lines', markersize=2)
    mpl.rc('mathtext', fontset='stix')

    mpl.rc('axes', titlesize=10, linewidth=1)
    mpl.rcParams['axes.labelweight'] = 'bold'

    mpl.rc('xtick', labelsize=10, direction='out', top=False)
    mpl.rc('xtick.major', top=False, )
    mpl.rc('xtick.minor', top=False, visible=False)

    mpl.rc('ytick', labelsize=10, direction='out', right=False)
    mpl.rc('ytick.major', right=True, )
    mpl.rc('ytick.minor', right=True, visible=False)


if __name__ == '__main__':
    # ---- Load user arguments ------------------------------------------------
    inputs = parse_arguments()
    sms = '\n\n ATTENTION !!! No trajectory passed.Run with -h for help.'
    assert inputs.trajectory, sms

    # ---- Stage 1: Load trajectory -------------------------------------------
    trajectory = load_trajectory(inputs)
    # ---- Stage 2: Calculate RMSD matrix -------------------------------------
    matrix, degrees = calc_rmsd_matrix(trajectory, inputs)
    # ---- Stage 3: Bitwise clustering ----------------------------------------
    clusters, leaders = bitclusterize(matrix, degrees)
    rmsd_ = calc_rms_vectors(trajectory, inputs)
    # --- Stage 4: Analysis ---------------------------------------------------
    os.mkdir('results')
    os.chdir('results')
    frames_stats = get_frames_stats(clusters, leaders)
    clusters_stats = get_cluster_stats(clusters, leaders)
    os.chdir('..')
    # =========================================================================
    # graph 1: rmsd_vs_reference (A)
    # =========================================================================
    generic_matplotlib()
    plt.ylabel('RMSD ($\AA$)', fontsize=10)
    plt.xlabel('Frame ID', fontsize=10)
    plt.plot(frames_stats.rmsd, lw=0.1, ls='-', color='navy', alpha=0.85)
    plt.savefig('rmsd_all_vs_reference', dpi=300, bbox_inches='tight',
                alpha=0.85)
    plt.close()
    # only for article publication
    run('mogrify -resize 1000x750 -density 300 -units PixelsPerInch rmsd_all_vs_reference.png', shell=True)

    # =========================================================================
    # graph 2: clusterlines
    # =========================================================================
    generic_matplotlib()
    plt.ylabel('Cluster ID', fontsize=10)
    plt.xlabel('Frame ID', fontsize=10)
    plt.plot(frames_stats.cluster_id, lw=0, marker='+', ms=1, color='navy',
             alpha=0.85)
    plt.savefig('clusters_lines', dpi=300, bbox_inches='tight')
    plt.close()
    run('mogrify -resize 1000x750 -density 300 -units PixelsPerInch clusters_lines.png', shell=True)

    # =========================================================================
    # graph 3: clustersizes
    # =========================================================================
    generic_matplotlib()
    plt.ylabel('Cluster size (%)', fontsize=10)
    plt.xlabel('Cluster ID', fontsize=10)
    plt.ylim(-1, clusters_stats['percent'].max()+5)
    plt.bar(clusters_stats.cluster_id,
            clusters_stats['percent'], color='navy', width=0.8, alpha=0.85)
    if -1 in leaders:
        plt.bar(clusters_stats.cluster_id[0], clusters_stats['percent'][0],
                color='r', width=0.8, label='Outliers', alpha=1)
    else:
        plt.bar(clusters_stats.cluster_id[0], clusters_stats['percent'][0],
                color='navy', width=0.8, label='Outliers', alpha=1)
    plt.legend(loc='upper right')
    plt.savefig('clusters_sizes', dpi=300, bbox_inches='tight')
    plt.close()
    run('mogrify -resize 1000x750 -density 300 -units PixelsPerInch clusters_sizes.png', shell=True)



    # =========================================================================
    # graph 4: colorbar
    # =========================================================================
    generic_matplotlib()
#    mpl.rc('axes', labelsize=12)
#    mpl.rc('font', family='monospace', size=12)
    mpl.rc('figure', autolayout=True, figsize=[6.197, 2.5], dpi=300)
    cmap = 'rainbow'
    data = np.asmatrix(frames_stats.cluster_id)
    plt.locator_params(axis='y', nbins=3)
    gs1 = mpl.gridspec.GridSpec(8, 1)
    gs1.update(hspace=0.08, wspace=0.01)
    ax1 = plt.subplot(gs1[:1, :1])
#    ax2 = plt.subplot(gs1[1:2, :1])
#    ax2.plot(range(0,1), color=cmap)
    plt.xlabel('Cluster ID', fontsize=8)
    im = plt.imshow(data, aspect='auto', interpolation='none', cmap=cmap)
#    cbar = plt.colorbar()
#    cbar.set_ticks([])
    plt.tight_layout()
    plt.tick_params(axis="y", which='both', left=False, right=False,
                    labelleft=False)
    plt.tick_params(axis="x", direction="out", which='both', top=False)
    colors_list = (im.cmap(im.norm(np.unique(data))))
    plt.savefig('clusters_colorbar', dpi=300, bbox_inches='tight')
    plt.close()
