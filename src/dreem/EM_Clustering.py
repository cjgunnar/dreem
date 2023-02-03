#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# The MIT License (MIT)
# Copyright (c) <2019> <The Whitehead Institute for Biomedical Research>

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
Created on Tue Jul 30 2019

@author: harish

Step 2 of the DREEM pipeline: EM Clustering of bit vectors.
"""
import argparse
import time
import os
import BitVector_Functions
import EM_Plots
import EM_CombineRuns
import Run_EMJobs
import EM_Files


def EM_Clustering():
    """
    """
    print('Starting EM clustering...')

    for ref in refs_seq:  # Each seq in the ref genome

        start_time = time.time()

        # bvfile_basename = '{}_{}_{}_{}'.format(sample_name, ref, START, END)
        # outplot_dir = outfiles_dir + bvfile_basename + '/'
        # if not os.path.exists(outplot_dir):
        #     os.makedirs(outplot_dir)
        # else:  # Folder exists
        #     if os.path.exists(outplot_dir + 'log.txt'):  # Log file exists
        #         print('EM Clustering already done for', bvfile_basename)
        #         return

        wind_size = int(END) - int(START)
        norm_bases = int((wind_size * NORM_PERC_BASES) / 100)

        # Read the bit vector file and do the filtering
        input_file = bitvector_file
        X = EM_Files.Load_BitVectors(input_file, INFO_THRESH, SIG_THRESH,
                                     inc_TG)

        if X is None:
            print("Too little viable bitvectors!")
            return

        K = 1  # Number of clusters
        cur_BIC = float('inf')  # Initialize BIC
        BIC_failed = False  # While test is not passed
        while not BIC_failed and K <= MAX_K:
            print('Working on K =', K)

            RUNS = NUM_RUNS if K != 1 else 1  # Only 1 Run for K=1
            ITS = MIN_ITS if K != 1 else 10  # Only 10 iters for K=1

            for run in range(1, RUNS + 1):
                print('Run number:', run)
                Run_EMJobs.Run_EMJob(X, ITS, INFO_THRESH,
                                         CONV_CUTOFF, SIG_THRESH,
                                         outfiles_dir, K, CPUS, run)


            # Processing of results from the EM runs
            EM_CombineRuns.Post_Process(K, RUNS,
                                        cur_BIC, norm_bases, struct,
                                        ref_file, outfiles_dir)

            # Check BIC
            latest_BIC = EM_CombineRuns.Collect_BestBIC(K,
                                                        outfiles_dir)
            if latest_BIC > cur_BIC:  # BIC test has failed
                BIC_failed = True
            cur_BIC = latest_BIC  # Update BIC

            K += 1  # Move on to next K

        end_time = time.time()
        time_taken = round((end_time - start_time) / 60, 2)
        print('Time taken:', time_taken, 'mins')

        # Write params to log file
        EM_Plots.Log_File(sample_name, NUM_RUNS, MIN_ITS,
                          CONV_CUTOFF, INFO_THRESH, SIG_THRESH, inc_TG,
                          norm_bases, K - 2, time_taken, log_file_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='EM Clustering')
    parser.add_argument('sample_name', help='Name of sample')
    parser.add_argument('bitvector', help='path to the bitvector.txt file')
    parser.add_argument('ref_name', help='Name of reference genome')
    parser.add_argument('ref_file', help='path to .fasta reference file')
    parser.add_argument('start', type=int, help='Start coord for bit vector')
    parser.add_argument('end', type=int, help='End coord for bit vector')
    parser.add_argument('out', help='output folder for cluster data')

    parser.add_argument('--min-its', type=int, help='Min iterations per EM run, default=300', default=300)
    parser.add_argument('--info-thresh', help='Threshold for informative bits, default=0.05', default=0.05)
    parser.add_argument('--conv-cutoff', help='Diff in log like for convergence, default=0.5', default=0.5)
    parser.add_argument('--num-runs', type=int, help='Number of EM runs per K, default=10', default=10)
    parser.add_argument('--max-k', type=int, help='Max K to work on, default=3', default=3)
    parser.add_argument('--cpus', type=int, help='Num processors to use, default=1', default=1)
    parser.add_argument('--norm-perc-bases', type=int, help='Perc of bases to use for norm, default=10', default=10)
    parser.add_argument('--sig-thresh', help='Signal threshold, default=0.0005', default=0.005)
    parser.add_argument('--inc-TG', action='store_true', help='Include signal from Ts & Gs, default=off')
    parser.add_argument('--struct', action='store_true', help='Secondary structure prediction, default=off')
    parser.add_argument('--log', help='path and name of log file')

    args = parser.parse_args()
    sample_name = args.sample_name
    ref_name = args.ref_name
    START = args.start
    END = args.end
    MIN_ITS = args.min_its
    INFO_THRESH = float(args.info_thresh)
    CONV_CUTOFF = float(args.conv_cutoff)
    NUM_RUNS = int(args.num_runs)
    MAX_K = int(args.max_k)
    CPUS = int(args.cpus)
    NORM_PERC_BASES = int(args.norm_perc_bases)
    inc_TG = args.inc_TG
    SIG_THRESH = float(args.sig_thresh)
    struct = args.struct
    bitvector_file = args.bitvector
    log_file_name = args.log

    ref_file = args.ref_file
    refs_seq = BitVector_Functions.Parse_FastaFile(ref_file)  # Ref seqs
    outfiles_dir = args.out

    if not os.path.exists(outfiles_dir):
        os.makedirs(outfiles_dir)

    EM_Clustering()
