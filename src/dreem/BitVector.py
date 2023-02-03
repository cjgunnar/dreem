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

Step 1 of the DREEM pipeline: Conversion of reads to bit vectors.

Output: Text file containing all the bit vectors. One bit vector is created
per read pair. This and other files, such as read coverage and pop avg plots,
are created in separate output subdirectories.
IMPORTANT: Assumption of SAME start and end coords for each seq in ref genome.

Calls Bit_Vector_Functions.py and Bit_Vector_Outputs_New.py.
"""
import os
import argparse
import time
import datetime
import BitVector_Functions
import BitVector_Outputs


def Bit_Vectors():
    """
    Create bit vectors for a sample based on the ref seq
    """
    start_time = time.time()

    # Initialize plotting variables
    mod_bases, mut_bases, delmut_bases = {}, {}, {}
    info_bases, cov_bases = {}, {}
    files, num_reads = {}, {}
    for ref in refs_seq:  # Each seq in the ref genome file
        ref_seq = refs_seq[ref]
        num_reads[ref] = 0
        mod_bases[ref], mut_bases[ref], delmut_bases[ref] = {}, {}, {}
        info_bases[ref], cov_bases[ref] = {}, {}
        for base in bases:
            mod_bases[ref][base] = {}
            for pos in range(start, end + 1):
                mod_bases[ref][base][pos] = 0
        for pos in range(start, end + 1):
            mut_bases[ref][pos], delmut_bases[ref][pos] = 0, 0
            info_bases[ref][pos], cov_bases[ref][pos] = 0, 0
        # Write header lines to output text file
        # this is a problem if refs_seq is more than 1 as it
        # will overwrite this file
        files[ref] = open(bitvector_file, 'w')
        files[ref].write('@ref' + '\t' + ref_name + ';' + ref + '\t' +
                         ref_seq[start - 1:end] + '\n')
        files[ref].write('@coordinates:length' + '\t' + str(start) + ',' +
                         str(end) + ':' + str(end - start + 1) + '\n')
        files[ref].write('Query_name\tBit_vector\tN_Mutations\n')

    # Compute Bit Vectors
    print('Computing bit vectors...')
    Process_SamFile(sam_file, paired, refs_seq, start, end,
                    cov_bases, info_bases, mod_bases,
                    mut_bases, delmut_bases, num_reads, files)

    print('Writing to the output file...')

    for ref in refs_seq:
        files[ref].close()

    end_time = time.time()
    time_taken = str(round((end_time - start_time) / 60, 2))

    if plot_dir:
        print("Creating plots...")
        # Write to output file
        BitVector_Outputs.writeOutputFiles(sample_name, num_reads,
                                           bitvector_file, plot_dir,
                                           refs_seq, start, end, mod_bases,
                                           mut_bases, delmut_bases,
                                           info_bases, cov_bases)

    if log_file_name:
        # Write to log file
        log_file = open(log_file_name, 'w')
        log_file.write('Sample: ' + sample_name + '\n')
        log_file.write('Reference file: ' + ref_file + '\n')
        log_file.write('Reference genomes: ' + str(list(refs_seq.keys())) + '\n')
        log_file.write('Num surrounding bases for del: ' + str(SUR_BASES) + '\n')
        log_file.write('ASCII+' + str(encoding))
        log_file.write('Q score cutoff: ' + str(QSCORE_CUTOFF) + '\n')
        log_file.write('Time taken: ' + str(time_taken) + ' mins\n')
        log_file.write('Finished at: ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M") + '\n')
        log_file.close()

    print('Finished creating bit vectors.')


def Process_SamFile(sam_file, paired, refs_seq, start, end,
                    cov_bases, info_bases, mod_bases,
                    mut_bases, delmut_bases, num_reads, files):
    """
    Read SAM file and generate bit vectors.
    """
    ignore_lines = len(refs_seq.keys()) + 2
    sam_fileobj = open(sam_file, 'r')
    for line_index in range(ignore_lines):  # Ignore header lines
        sam_fileobj.readline()
    while True:
        try:
            if paired:
                line1, line2 = next(sam_fileobj), next(sam_fileobj)
                line1, line2 = line1.strip().split(), line2.strip().split()
                mate1 = BitVector_Functions.Mate(line1)
                mate2 = BitVector_Functions.Mate(line2)
                assert mate1.PNEXT == mate2.POS and \
                    mate1.RNAME == mate2.RNAME and mate1.RNEXT == "="
                assert mate1.QNAME == mate2.QNAME and mate1.MAPQ == mate2.MAPQ
                GenerateBitVector_Paired(mate1, mate2, refs_seq,
                                         phred_qscore, cov_bases, info_bases,
                                         mod_bases, mut_bases, delmut_bases,
                                         num_reads, files)
            else:
                line = next(sam_fileobj)
                line = line.strip().split()
                mate = BitVector_Functions.Mate(line)
                GenerateBitVector_Single(mate, refs_seq, phred_qscore,
                                         cov_bases, info_bases, mod_bases,
                                         mut_bases, delmut_bases, num_reads,
                                         files)
        except StopIteration:
            break
    sam_fileobj.close()


def GenerateBitVector_Paired(mate1, mate2, refs_seq, phred_qscore,
                             cov_bases, info_bases, mod_bases, mut_bases,
                             delmut_bases, num_reads, files):
    """
    Create a bitvector for paired end sequencing.
    """
    bitvector_mate1 = Convert_Read(mate1, refs_seq, phred_qscore)
    bitvector_mate2 = Convert_Read(mate2, refs_seq, phred_qscore)
    bit_vector = Combine_Mates(bitvector_mate1, bitvector_mate2)
    Plotting_Variables(mate1.QNAME, mate1.RNAME, bit_vector, start, end,
                       cov_bases, info_bases, mod_bases, mut_bases,
                       delmut_bases, num_reads, files)


def GenerateBitVector_Single(mate, refs_seq, phred_qscore,
                             cov_bases, info_bases, mod_bases,
                             mut_bases, delmut_bases, num_reads, files):
    """
    Create a bitvector for single end sequencing.
    """
    bit_vector = Convert_Read(mate, refs_seq, phred_qscore)
    Plotting_Variables(mate.QNAME, mate.RNAME, bit_vector, start, end,
                       cov_bases, info_bases, mod_bases, mut_bases,
                       delmut_bases, num_reads, files)


def Convert_Read(mate, refs_seq, phred_qscore):
    """
    Convert a read's sequence to a bit vector of 0s & 1s and substituted bases
    Args:
        mate (Mate): Read
        refs_seq (dict): Sequences of the ref genomes in the file
        phred_qscore (dict): Qual score - ASCII symbol mapping
    Returns:
        bitvector_mate (dict): Bitvector. Format: d[pos] = bit
    """
    bitvector_mate = {}  # Mapping of read to 0s and 1s
    read_seq = mate.SEQ  # Sequence of the read
    ref_seq = refs_seq[mate.RNAME]  # Sequence of the ref genome
    q_scores = mate.QUAL  # Qual scores of the bases in the read
    i = mate.POS  # Pos in the ref sequence
    j = 0  # Pos in the read sequence
    CIGAR_Ops = BitVector_Functions.Parse_CIGAR(mate.CIGAR)
    op_index = 0
    while op_index < len(CIGAR_Ops):  # Each CIGAR operation
        op = CIGAR_Ops[op_index]
        desc, length = op[1], int(op[0])

        if desc == 'M':  # Match or mismatch
            for k in range(length):  # Each base
                if phred_qscore[q_scores[j]] >= QSCORE_CUTOFF:
                    bitvector_mate[i] = read_seq[j] \
                        if read_seq[j] != ref_seq[i - 1] else nomut_bit
                else:  # < Qscore cutoff
                    bitvector_mate[i] = ambig_info
                i += 1  # Update ref index
                j += 1  # Update read index

        elif desc == 'D':  # Deletion
            for k in range(length - 1):  # All bases except the 3' end
                bitvector_mate[i] = ambig_info
                i += 1  # Update ref index
            ambig = BitVector_Functions.Calc_Ambig_Reads(ref_seq, i, length,
                                                         SUR_BASES)
            bitvector_mate[i] = ambig_info if ambig else del_bit
            i += 1  # Update ref index

        elif desc == 'I':  # Insertion
            j += length  # Update read index

        elif desc == 'S':  # Soft clipping
            j += length  # Update read index
            if op_index == len(CIGAR_Ops) - 1:  # Soft clipped at the end
                for k in range(length):
                    bitvector_mate[i] = miss_info
                    i += 1  # Update ref index
        else:
            print('Unknown CIGAR op encountered.')
            return ''

        op_index += 1
    return bitvector_mate


def Combine_Mates(bitvector_mate1, bitvector_mate2):
    """
    Combine bit vectors from mate 1 and mate 2 into a single read's bit vector.
    0 has preference. Ambig info does not. Diff muts in the two mates are
    counted as ambiguous info.
    Args:
        bitvector_mate1 (dict): Bit vector from Mate 1
        bitvector_mate2 (dict): Bit vector from Mate 2
    Returns:
        bit_vector (dict): Bitvector. Format: d[pos] = bit
    """
    bit_vector = {}
    for (pos, bit) in bitvector_mate1.items():  # Bits in mate 1
        bit_vector[pos] = bit
    for (pos, bit) in bitvector_mate2.items():  # Bits in mate2
        if pos not in bitvector_mate1:  # Not present in mate 1
            bit_vector[pos] = bit  # Add to bit vector
        else:  # Overlap in mates
            mate1_bit = bitvector_mate1[pos]
            mate2_bit = bitvector_mate2[pos]
            bits = set([mate1_bit, mate2_bit])
            if len(bits) == 1:  # Both mates have same bit
                bit_vector[pos] = mate1_bit
            else:  # More than one bit
                if nomut_bit in bits:  # 0 in one mate
                    bit_vector[pos] = nomut_bit  # Add 0
                elif ambig_info in bits:  # Ambig info in one mate
                    other_bit = list(bits - set(ambig_info))[0]
                    bit_vector[pos] = other_bit  # Add other bit
                elif mate1_bit in bases and mate2_bit in bases:
                    if mate1_bit != mate2_bit:  # Diff muts on both mates
                        bit_vector[pos] = ambig_info
    return bit_vector


def Plotting_Variables(q_name, ref, bit_vector, start, end, cov_bases,
                       info_bases, mod_bases, mut_bases, delmut_bases,
                       num_reads, files):
    """
    Create final bit vector in relevant coordinates and all the
    variables needed for plotting
    Args:
        q_name (string): Query name of read
        ref (string): Name of ref genome
        bit_vector (dict): Bit vector from the mate/mates
        start (int): start position in reference (1-based)
        end (int): end position in reference (1-based)
    """
    # Create bit vector in relevant coordinates
    num_reads[ref] += 1  # Add a read to the count
    bit_string = ''
    for pos in range(start, end + 1):  # Each pos in coords of interest
        if pos not in bit_vector:  # Pos not covered by the read
            read_bit = miss_info
        else:
            read_bit = bit_vector[pos]
            cov_bases[ref][pos] += 1
            if read_bit != ambig_info:
                info_bases[ref][pos] += 1
            if read_bit in bases:  # Mutation
                mod_bases[ref][read_bit][pos] += 1
                mut_bases[ref][pos] += 1
                delmut_bases[ref][pos] += 1
            elif read_bit == del_bit:  # Deletion
                delmut_bases[ref][pos] += 1
        bit_string += read_bit
    # Write bit vector to output text file
    n_mutations = str(float(sum(bit.isalpha() for bit in bit_string)))
    if not bit_string.count('.') == len(bit_string):  # Not all '.'
        files[ref].write(q_name + '\t' + bit_string + '\t' + n_mutations + '\n')


def main():
    parser = argparse.ArgumentParser(description='Creation of bit vectors')
    parser.add_argument('sample_name', help='Name of sample')
    parser.add_argument('sample_sam', help='.sam file of sample')
    parser.add_argument('ref_name', help='Name of reference genome')
    parser.add_argument('ref_file', help='Fasta file of reference')
    parser.add_argument('start', type=int, help='Start pos in ref genome (1-based)')
    parser.add_argument('end', type=int, help='Start pos in ref genome (1-based)')
    parser.add_argument('out', help='bitvector file output (include .txt)')
    parser.add_argument('plots', help='directory for plots and log for describing bitvector. If nothing is provided '
                                      'the plots will not be generated', nargs='?')
    parser.add_argument('--sur-bases', type=int, help='Bases surrounding a deletion, default=10', default=10)
    parser.add_argument('--ascii-64', help='Use ASCII+64 char to Q score map instead of the default ASCII+33',
                        action='store_true')
    parser.add_argument('--qscore-cutoff', type=int, help='Qscore cutoff for a valid base, default=20', default=20)
    parser.add_argument('--paired', help='Paired-end sequencing, default=off', action='store_true')
    parser.add_argument('-l', '--log', help='location of log file to write to, default=no log')
    args = parser.parse_args()
    sample_name = args.sample_name
    ref_name = args.ref_name
    start = args.start
    end = args.end
    SUR_BASES = args.sur_bases
    QSCORE_CUTOFF = args.qscore_cutoff
    bitvector_file = args.out
    plot_dir = args.plots
    paired = args.paired
    log_file_name = args.log

    # Symbols to represent a read as a bit vector
    miss_info, ambig_info = '.', '?'
    nomut_bit, del_bit = '0', '1'
    bases = ['A', 'T', 'G', 'C']

    # output paths
    bitvector_root = os.path.split(bitvector_file)[0]
    if not os.path.exists(bitvector_root):
        os.makedirs(bitvector_root)
    if plot_dir and not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    # Generate input variables from input file
    ref_file = args.ref_file
    refs_seq = BitVector_Functions.Parse_FastaFile(ref_file)  # Ref seqs

    encoding = 33 if not args.ascii_64 else 64

    # generate map of ascii characters to qscore
    # where '!' (ascii=33) = 0 in ASCII+33 encoding
    # see https://people.duke.edu/~ccc14/duke-hts-2018/bioinformatics/quality_scores.html
    # not sure why the max is 41 but that qscore represents 0.008% probability
    # so they probably just don't get more precise than that
    max_qscore = 41
    phred_qscore = dict((chr(i + encoding), i) for i in range(0, max_qscore + 1))

    sam_file = args.sample_sam

    Bit_Vectors()


if __name__ == '__main__':
    main()
