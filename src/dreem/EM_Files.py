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

Read the bitvector text file and create object X
Read the FASTA file to get the ref genome sequence
Filtering of the bit vectors is done here using various criteria
Changing of . and ? to 0 is done here
"""
import EM_Class
import EM_Functions
import numpy as np


def Load_BitVectors(bv_file, INFO_THRESH, SIG_THRESH, inc_TG):
    """
    """
    bases = ['A', 'T', 'G', 'C']
    bit_strings, mut_popavg, n_discard = [], {}, 0
    f, f1, f2, f3, f4 = 0, 0, 0, 0, 0

    bv_fileobj = open(bv_file)
    bvfile_contents = bv_fileobj.readlines()
    bv_fileobj.close()

    first_line = bvfile_contents[0]
    first_line_split = first_line.strip().split()
    ref_info, seq = first_line_split[1], first_line_split[2]
    ref_file, ref = ref_info.split(';')[0], ref_info.split(';')[1]

    second_line = bvfile_contents[1]
    second_line_split = second_line.strip().split()
    indices = second_line_split[1].split(':')[0]

    l = len(bvfile_contents[3].strip().split()[1])  # Len of 1st bit string
    nmuts_min = int(round(0.1 * l))
    nmuts_thresh = max(nmuts_min, EM_Functions.calc_nmuts_thresh(bv_file))
    print('Mutations threshold:', nmuts_thresh)

    for i in range(3, len(bvfile_contents)):
        f += 1
        line = bvfile_contents[i].strip().split()
        bit_string = line[1]
        n_mut = float(line[2])

        # Replace bases with 1
        for base in bases:
            bit_string = bit_string.replace(base, '1')

        # Filter 1 - Number of mutations
        if n_mut > nmuts_thresh:
            n_discard += 1
            f1 += 1
            continue

        # Filter 2 - Fraction of informative bits
        if (bit_string.count('.') + bit_string.count('?') +
           bit_string.count('N')) >= INFO_THRESH * len(bit_string):
            n_discard += 1
            f2 += 1
            continue

        # Filter 3 - Distance between mutations
        if not EM_Functions.is_distmuts_valid(bit_string):
            n_discard += 1
            f3 += 1
            continue

        # Filter 4 - Bits surrounding mutations
        if not EM_Functions.is_surmuts_valid(bit_string):
            n_discard += 1
            f4 += 1
            continue

        bit_strings.append(bit_string)

    print('Total bit vectors:', f)
    print('Bit vectors removed because of too many mutations: ', f1)
    print('Bit vectors removed because of too few informative bits: ', f2)
    print('Bit vectors removed because of mutations close by: ', f3)
    print('Bit vectors removed because of no info around mutations: ', f4)

    minimum_bitvectors = 30
    if len(bit_strings) < minimum_bitvectors:
        return None

    D = len(bit_strings[0])
    thresh_pos = []  # Positions below signal threshold
    for d in range(D):  # Each position of interest in the genome
        bits_list = [bs[d] for bs in bit_strings]  # List of bits at that pos
        noinfo_count = bits_list.count('.') + bits_list.count('?') + \
            bits_list.count('N')
        info_count = len(bits_list) - noinfo_count  # Num of informative bits
        try:
            mut_prob = bits_list.count('1') / info_count
        except ZeroDivisionError:
            mut_prob = 0
        if mut_prob < SIG_THRESH:
            mut_prob = 0
            thresh_pos.append(d)
        mut_popavg[d] = mut_prob

    for i in range(len(bit_strings)):  # Change . and ? to 0, noise to 0
        bit_string = bit_strings[i]
        bit_string = bit_string.replace('?', '0')
        bit_string = bit_string.replace('.', '0')
        bit_string = bit_string.replace('N', '0')

        # Suppressing data from Ts and Gs
        if inc_TG:
            bit_string = list(bit_string)
            j = 0
            while j < len(bit_string):
                if seq[j] == 'T' or seq[j] == 'G':
                    bit_string[j] = '0'
                j += 1
            bit_string = ''.join(bit_string)

        bit_string = np.array(list(bit_string))
        bit_string[thresh_pos] = '0'
        bit_string = ''.join(bit_string)

        bit_strings[i] = bit_string

    X = EM_Class.BV_Object(bit_strings, mut_popavg, n_discard, ref_file,
                           ref, seq, indices)
    return X
