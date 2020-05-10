#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
.. note ::

  | **Author      :** Roy Gonzalez Aleman
  | **Contact     :** [roy_gonzalez@fq.uh.cu, roy.gonzalez-aleman@u-psud.fr]

**Description:**
  'BitClust: Fast & memory efficient clustering of MD trajectories'

'''
import os
import sys
import shutil
import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import OrderedDict

import mdtraj as md
import bitarray.util as bau
from bitarray import bitarray as ba

def parse_arguments():
    '''
    DESCRIPTION
    Parse all user arguments from the command line.

    Return:
        user_inputs (parser.argparse): namespace with user input arguments.
    '''

    # Initializing argparse ---------------------------------------------------
    desc = '\nBitClust: Fast & memory efficient clustering of long MD trajectories'
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
                        type=int, required=False, default=None)
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
                        help='minimum number of frames inside returned clusters',
                        type=int, required=False, default=2)
    parser.add_argument('-max_clust', action='store', dest='max_clust',
                        help='maximum number of returned clusters',
                        type=int, required=False, default=np.inf)
    parser.add_argument('-ref', action='store', dest='reference',
                        help='reference frame to align trajectory',
                        type=int, required=False, default=0)

    # Arguments: analysis -----------------------------------------------------
    parser.add_argument('-numouts', action='store', dest='num_outs',
                        help='number of clusters and leaders to print as pdb',
                        type=int, required=False, default=5)
    parser.add_argument('-odir', action='store', dest='outdir',
                        help='output directory to store analysis',
                        type=str, required=False, default='./')

    user_inputs = parser.parse_args()
    return user_inputs


# Debuging jobs can start from here
# args = argparse.Namespace()
# args.topology = 'examples/aligned_tau.pdb'
# args.trajectory = 'examples/aligned_original_tau_6K.dcd'
# args.selection = 'all'
# args.cutoff = 4
# args.minsize = 2
# args.reference = 0
# args.outdir = './out'
# args.first = 0
# args.last = 100
# args.stride = 1

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
    orig_traj_idx = range(trajectory.n_frames)
    # Reduce RAM consumption by loading selected atoms only -------------------
    if args.selection != 'all':
        try:
            sel_indx = trajectory.topology.select(args.selection)
        except ValueError:
            print('Specified selection is invalid')
            sys.exit()
        if sel_indx.size == 0:
            print('Specified selection in your system corresponds to no atoms')
            sys.exit()
        trajectory = trajectory.atom_slice(sel_indx)[args.first: args.last:
                                                     args.stride]
    else:
        trajectory = trajectory[args.first:args.last:args.stride]
    trajectory_idx = orig_traj_idx[args.first:args.last:args.stride]
    # Center coordinates of loaded trajectory ---------------------------------
    trajectory.center_coordinates()
    return trajectory, trajectory_idx


def np_to_bitarray(np_array, N):
    '''
    DESCRIPTION
    Converts a numpy array to a bitarray using the fastest way. It creates a
    bitarray of N bits and sets to 1 only those indices that coincides with
    the integers in the numpy array.

    Arguments:
        np_array (numpy.array): numpy array.
        N (int): size of the desired bitarray.
    Return:
        bitarr (bitarray): a bitarray of N bits having as 1 only those indices
        that coincides with the integers present in the numpy array.
    '''
    zero_arr = np.zeros(N, dtype=np.bool)
    zero_arr[np_array] = 1
    bitarr = ba()
    bitarr.pack(zero_arr.tobytes())
    return bitarr


def calc_rmsd_matrix(trajectory, args):
    '''
    DESCRIPTION
    Calculates RMSD square matrix using MDTraj. Pairwise similarity is saved in
    RAM as bits (dict of bitarrays), not floats.

    Args:
        trajectory (mdtraj.Trajectory): trajectory object to analyze.
        args (argparse.Namespace): user input parameters parsed by argparse.
    Return:
        matrix (list): list of bitarrays
        degrees (collections.OrderedDict): dict of bitarrays´ lenghts.
    '''
    trajectory.center_coordinates()
    N = trajectory.n_frames
    cutoff = np.float32(args.cutoff)/10  # transform A to nm (MDTraj coherence)
    matrix = []
    degrees = OrderedDict()
    to_explore = range(N)
    for i in to_explore:
        rmsd_ = md.rmsd(trajectory, trajectory, i, precentered=True)
        vector_np = np.less_equal(rmsd_, cutoff)
        bitarr = np_to_bitarray(vector_np, N)
        matrix.append(bitarr)
        degrees.update({i: bitarr.count()})
    return matrix, degrees


def bitclusterize(matrix, degrees, args):
    '''
    DESCRIPTION
    Clusters the bit matrix using bitwise operations.

    Args:
        matrix (list): list of bitarrays
        degrees (collections.OrderedDict): dict of bitarrays lenghts.
    Return:
        clusters (numpy.ndarray): array of clusters ID.
        leaders (list) : list of clusters´ centers ID.
    '''
    degrees = np.asarray([degrees[x] for x in degrees])
    # Memory allocation for clusters container --------------------------------
    clusters = np.empty(len(matrix), dtype='int32')
    clusters.fill(-1)
    # Declare all bits as available -------------------------------------------
    available_bits = ba(int(len(degrees)))
    available_bits.setall('1')
    # Start iterative switching -----------------------------------------------
    leaders = []
    clust_id = 0
    ncluster = -1
    while True:
        ncluster += 1
        # Break 0: break if max number of cluster was reached -----------------
        if ncluster > args.max_clust:
            break
        # Find the biggest cluster --------------------------------------------
        leader = degrees.argmax()
        # Break 1: all candidates cluster have degree 1 (can´t clusterize) ----
        if degrees.sum() == np.nonzero(degrees)[0].size:
            return clusters, leaders
            break
        biggest_cluster = matrix[leader] & available_bits
        biggest_cluster_list = np.frombuffer(biggest_cluster.unpack(),
                                             dtype=np.bool)
        # Break 2: all candidates cluster have degree < minsize ---------------
        if biggest_cluster_list.sum() < args.minsize:
            return clusters, leaders
            break
        # Break 3: No more candidates available (empty matrix) ----------------
        if degrees.sum() == 0:
            return clusters, leaders
            break
        degrees[biggest_cluster_list] = 0
        available_bits = (available_bits ^ biggest_cluster) & available_bits
        if biggest_cluster.count() <= 1:
            leaders.append(-1)
            return clusters, leaders
            break
        else:
            leaders.append(leader)
            # Assign next cluster ID ------------------------------------------
            clusters[biggest_cluster_list] = clust_id
            clust_id += 1
        # Update degrees of unclustered frames --------------------------------
        for degree in available_bits.itersearch(ba('1')):
            # degrees[degree] = ba.fast_hw_and(available_bits, matrix[degree])
            degrees[degree] = bau.count_and(available_bits, matrix[degree])



def calc_rms_vectors(trajectory, args):
    '''
    DESCRIPTION
    Calculates RMSD of all frames vs. referenece.
    Args:
        trajectory (mdtraj.Trajectory): trajectory object to analyze.
        args (argparse.Namespace): user input parameters parsed by argparse.
    Return:
        rmsd_ (list): list containing the RMSD of all frames vs. reference.
    '''
    trajectory.center_coordinates()
    reference = args.reference
    rmsd_ = md.rmsd(trajectory, trajectory, reference, precentered=True)
    return rmsd_


def get_frames_stats(clusters, rmsd_):
    '''
    DESCRIPTION
    Gets plain text file "frames_statistics.txt" which contains as columns:
    frameID, clusterID and RMSD distance with respect to reference.

    Args:
        clusters (numpy.ndarray): array of clusters ID.
        rmsd_ (list): list containing the RMSD of all frames vs. reference.

    Return:
        frames_df (pandas.DataFrame): dataframe with frames_statistics info.
    '''
    frames_df = pd.DataFrame(columns=['frame', 'cluster_id', 'rmsd'])
    frames_df['frame'] = range(len(clusters))
    frames_df['cluster_id'] = clusters
    frames_df['rmsd'] = rmsd_*10
    return frames_df


def get_cluster_stats(clusters, leaders):
    '''
    DESCRIPTION
    Gets plain text file "cluster_statistics.txt" which contains as columns:
    clusterID, cluster_size, cluster_leader and cluster percentage from
    trajectory.

    Args:
        clusters (numpy.ndarray): array of clusters ID.
        leaders (list) : list of clusters´ centers ID.

    Return:
        clusters_df (pandas.DataFrame): dataframe with cluster_statistics info.
    '''
    clusters_df = pd.DataFrame(columns=['cluster_id', 'size', 'leader',
                                        'percent'])

    # treatment in case of outliers
    if len(np.where(clusters == -1)[0]):
        leaders = leaders + [-1]
        clusters_df['leader'] = leaders
        clusters_df['cluster_id'] = list(range(0, len(leaders) - 1)) + [-1]
        sizes = []
        for x in clusters_df.cluster_id:
            sizes.append(len(np.where(clusters == x)[0]))
        clusters_df['size'] = sizes
    # treatment in case of non outliers
    else:
        clusters_df['leader'] = leaders
        clusters_df['cluster_id'] = list(range(0, len(leaders)))
        sizes = []
        for x in clusters_df.cluster_id:
            sizes.append(len(np.where(clusters == x)[0]))
        clusters_df['size'] = sizes

    sum_ = clusters_df['size'].sum()
    percents = [round(x/sum_*100, 4) for x in clusters_df['size']]
    clusters_df['percent'] = percents
    return clusters_df


def generic_matplotlib():
    '''
    DESCRIPTION
    Some customizations of matplotlib.
    '''
    mpl.rc('figure', autolayout=True, figsize=[7, 5], dpi=300)
    mpl.rc('font', family='STIXGeneral')
    mpl.rc('lines', markersize=2)
    mpl.rc('mathtext', fontset='stix')

    mpl.rc('axes', titlesize=10, linewidth=1)
    mpl.rcParams['axes.labelweight'] = 'bold'

    mpl.rc('xtick', labelsize=12, direction='out', top=False)
    mpl.rc('xtick.major', top=False, )
    mpl.rc('xtick.minor', top=False, visible=False)

    mpl.rc('ytick', labelsize=12, direction='out', right=False)
    mpl.rc('ytick.major', right=True, )
    mpl.rc('ytick.minor', right=True, visible=False)


def main():

    # ---- Load user arguments ------------------------------------------------
    inputs = parse_arguments()
    sms = '\n\n ATTENTION !!!! No trajectory passed.Run with -h for help.'
    assert inputs.trajectory, sms

    # ---- Stage 1: Load trajectory -------------------------------------------
    trajectory, trajectory_idx = load_trajectory(inputs)
    # ---- Stage 2: Calculate RMSD matrix -------------------------------------
    matrix, degrees = calc_rmsd_matrix(trajectory, inputs)
    # ---- Stage 3: Bitwise clustering ----------------------------------------
    clusters, leaders = bitclusterize(matrix, degrees, inputs)
    rmsd_ = calc_rms_vectors(trajectory, inputs)

    if inputs.outdir == './':
        pass
    else:
        try:
            os.makedirs(inputs.outdir)
            os.chdir(inputs.outdir)
        except FileExistsError:
            shutil.rmtree(inputs.outdir)
            os.makedirs(inputs.outdir)
            os.chdir(inputs.outdir)

    # ---- Stage 4: Output clusters/leaders -----------------------------------
    nouts = inputs.num_outs
    for c in range(nouts):
        clust_out = np.where(clusters == c)[0]
        if clust_out.size == 0:
            break
        traj_out = trajectory[clust_out]
        clust_name = 'cluster_{}.pdb'.format(c)
        traj_out.save_pdb(clust_name)

        leader_out = trajectory[leaders[c]]
        lead_name = 'leader_{}.pdb'.format(c)
        leader_out.save_pdb(lead_name)

    # # --- Stage 5: Analysis -------------------------------------------------
    # log files
    frames_stats = get_frames_stats(clusters, rmsd_)
    clusters_stats = get_cluster_stats(clusters, leaders)
    # log files conversion to original indices
    frames_stats.frame = trajectory_idx
    if clusters_stats.leader.iloc[-1] == -1:
        clusters_stats.leader[:-1] = np.asarray(trajectory_idx)[clusters_stats.leader[:-1]]
    else:
        clusters_stats.leader = np.asarray(trajectory_idx)[clusters_stats.leader]
    # writing log files in original indices
    with open('frames_statistics.txt', 'wt') as on:
        frames_stats.to_string(buf=on, index=False)
    with open('cluster_statistics.txt', 'wt') as on:
        clusters_stats.to_string(buf=on, index=False)

    # graph 1: rmsd_vs_reference (A)
    generic_matplotlib()
    mpl.rc('figure', autolayout=True, figsize=[7, 3.5], dpi=300)
    mpl.rc('xtick', labelsize=18, direction='out', top=False)
    mpl.rc('ytick', labelsize=18, direction='out', right=False)
    plt.ylabel('RMSD ($\AA$)', fontsize=18)
    plt.xlabel('Frame ID', fontsize=18)
    colors = ['r', 'orange', 'green', 'blue', 'purple']
    plt.scatter(frames_stats.frame, frames_stats.rmsd, marker='x', s=8,
                color='gray', alpha=1, label='RMSD vs. reference')
    plt.savefig('rmsd_all_vs_reference', dpi=300,
                alpha=0.85)
    plt.close()

    # graph 2: clusterlines
    generic_matplotlib()
    plt.ylabel('Cluster ID', fontsize=14)
    plt.xlabel('Frame ID', fontsize=14)
    plt.plot(frames_stats.frame, frames_stats.cluster_id, lw=0, marker='+',
             ms=1, color='k', alpha=0.85)
    plt.savefig('clusters_lines', dpi=300, bbox_inches='tight')
    plt.close()

    # graph 3: clustersizes
    generic_matplotlib()
    plt.ylabel('Cluster size (%)', fontsize=14)
    plt.xlabel('Cluster ID', fontsize=14)
    plt.ylim(-1, clusters_stats['percent'].max()+5)
    plt.bar(clusters_stats.cluster_id,
            clusters_stats['percent'], color='k', width=0.8, alpha=0.85)
    plt.bar(clusters_stats.cluster_id.iloc[-1],
            clusters_stats['percent'].iloc[-1],
            color='r', width=0.8, label='Unclustered', alpha=1)
    plt.legend(loc='upper right')
    plt.savefig('clusters_sizes', dpi=300, bbox_inches='tight')
    plt.close()

    # graph 4: coloRMSD colored by clusters
    generic_matplotlib()
    mpl.rc('figure', autolayout=True, figsize=[7, 3.5], dpi=300)
    mpl.rc('xtick', labelsize=18, direction='out', top=False)
    mpl.rc('ytick', labelsize=18, direction='out', right=False)
    plt.ylabel('RMSD ($\AA$)', fontsize=18)
    plt.xlabel('Frame ID', fontsize=18)
    colors = ['r', 'orange', 'green', 'blue', 'purple']
    plt.scatter(frames_stats.frame, frames_stats.rmsd, marker='x', s=8,
                color='gray', alpha=0.3, label='RMSD')
    for i, c in enumerate(colors):
        frames = frames_stats[frames_stats.cluster_id == i].frame
        rms = frames_stats[frames_stats.cluster_id == i].rmsd
        plt.scatter(frames, rms, color=c, marker='+', s=8, alpha=0.5,
                    label='Clust-{}'.format(i+1))
        plt.legend(bbox_to_anchor=(0.95, 1), markerscale=5,
                   fontsize='xx-large')
    plt.savefig('rmsd_all_vs_reference_colored', dpi=300)
    plt.close()

    print('\n\n\nNORMAL TERMINATION ;)')

if __name__ == '__main__':
    main()
