#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 nylander <johan.nylander@nbis.se>
# Distributed under terms of the MIT license.
# Last modified: tor dec 08, 2022  06:24

"""
Description:

    Parse BAM file and fasta file (reference) and compare the matching
    sequences according to the matching scheme (below).

    Output is a count of matches for a specific ref seq in BAM.

Usage:

    parse.py -b BAMFILE -f FASTAFILE [-g [GENELENGTH]] [-m [MISMATCHES]]
             [-t [THREADS]] [-o [OUTPUT]] [-v] [-V] [-h]

Examples:

    ./parse.py -b my.sorted.bam -f ref.fas -o counts.tsv
    ./parse.py -b my.sorted.bam -f ref.fas -g 9 -m 2 -t 4 -v -o counts.tsv

Options:

    -h, --help              show this help message and exit
    -b, --bam BAM           bam file
    -f, --fasta FASTA       fasta (references) file
    -g, --genelength [GENELENGTH] length of leading/trailing gene segment
    -m, --mismatches [MISMATCHES] number of accepted mismatches in core sequence
    -t, --threads [THREADS] number of threads
    -o, --output [OUTPUT] Output file (default: standard output)
    -v, --verbose         increase output verbosity
    -V, --version         show program's version number and exit

Prerequisites:

    tab2fasta list.tab > ref.fas
    bwa index ref.fas
    bwa mem -t 8 ref.fas query.fastq.gz | \\
        samtools sort -@8 -o my.sorted.bam -
    samtools index -@8 my.sorted.bam

Matching scheme:

    GTTCGCCATCTTCAGACTACTAGTTTGATTCGACCCACCGATACCGTTCGTATAATGTATGCTATACGAACGGTATCTCCAGGTCGAATCAAACTAGTTCCCTGAAGATGTCGATG  (Reference sequence)
    GTTCGCCAT (First gene, 9 bp)                                                                               ATGTCGATG  (Second gene, 9 bp)
             CTTCAG  (Non-biological seq, 6 bp)                                                          CTGAAG           (Non-biological seq, 6 bp)
                   ACT (First index, 3 bp. Either ACT or GGA)                                         TCC                 (Second index, 3 bp. Either AGT or TCC (wich are RC of first))
                      ACTAGTTTGATTCGACCCACCGATACCGTTCGTATAATGTATGCTATACGAACGGTATCTCCAGGTCGAATCAAACTAGT                    (Core matching sequence. Allow m mismatches (0-3?))



"""

import sys
import re
import argparse
import pysam
from pyfaidx import Fasta


# Defaults
__version__ = '0.2'
gene_length_default = 9
mismatches_default = 0
threads_default = 4
index_pos_one = 15
index_pos_two = 18
ref_core_start = 18
ref_core_end = -18
allowed_first_indices = ['ACT', 'GGA']
allowed_second_indices = ['AGT', 'TCC'] # revcomp of first


def hammingDist(seq1, seq2):
    """
    Calculate the Hamming distance in terms of number of changes
    Args: seq1 and seq2 (str)
    Return: distance (int)
    TODO: do we need upper()?
    """
    i = 0
    count = 0
    while (i < len(seq1)):
        if (seq1[i] != seq2[i]):
            count += 1
        i += 1
    return count


def getGene(alignedSegmentObject, genelength, first):
    """
    Get the "gene" sequence of length genelength in the first
    or last part of the sequence. No mismatches allowed.
    Args: alignedSegmentObject, genelength, first (boolean)
    Return: gene sequence (str)
    TODO: rewrite to get ref seq from ref seq and not from aligned pairs!
    """
    aso = alignedSegmentObject
    sp = aso.get_aligned_pairs(matches_only = True)
    if (sp[0][1] != aso.reference_start):
        return(0)

    if (sp[-1][1] != aso.reference_end - 1):
        return(0)

    if (first):
        g = genelength
        ref_gene_start = aso.reference_start
        if (sp[g][1] == g):
            ref_gene_end = g
            ref_gene_seq = aso.get_reference_sequence()[ref_gene_start:ref_gene_end]
        else:
            #if (debug):
            #    print('Warning: getGene: first ref_gene_seq not of length 9 ', file = sys.stderr)
            return(0)
        query_gene_start = sp[ref_gene_start][0]
        query_gene_end = sp[ref_gene_end][0]
        if ((query_gene_end - query_gene_start) == g):
            query_gene_seq = aso.query_sequence[query_gene_start:query_gene_end]
            #query_gene_seq = aso.get_forward_sequence()[query_gene_start:query_gene_end]
        else:
            #if (debug):
            #    print('Warning: getGene: first query_gene_seq not of length 9', file = sys.stderr)
            return(0)
    else:
        g = -genelength
        ref_gene_start = sp.index(sp[:][g])
        ref_gene_end = aso.reference_end - 1
        if ((ref_gene_start - ref_gene_end - 1) == g):
            ref_gene_seq = aso.get_reference_sequence()[ref_gene_start:ref_gene_end + 1]
        else:
            #if (debug):
            #    print('Warning: getGene: last ref_gene_seq not of length 9', file = sys.stderr)
            return(0)
        query_gene_start = sp[ref_gene_start][0]
        query_gene_end = sp[ref_gene_end][0]
        if ((query_gene_start - query_gene_end - 1) == g):
            query_gene_seq = aso.query_sequence[query_gene_start:query_gene_end + 1]
        else:
            #if (debug):
            #    print('Warning: getGene: last query_gene_seq not of length 9', file = sys.stderr)
            return(0)

    if (query_gene_seq == ref_gene_seq):
        return(query_gene_seq)
    else:
        return(0)


def checkIndex(alignedSegmentObject, i1, i2, allowed_indices, first):
    """
    Locate first and last indices.
    Allowed triplets: ACT or GGA (and RC: AGT or TCC)
    Indexes are located at position 15:18 and -18:-15
    Args:
    Return: True or False (or should we return the actual seq?)
    TODO:
    """
    aso = alignedSegmentObject
    sp = aso.get_aligned_pairs(matches_only = True)

    if (first):
        ref_index_start = sp.index(sp[:][i1])
        ref_index_end = sp.index(sp[:][i2])
        ref_index_seq = aso.get_reference_sequence()[ref_index_start:ref_index_end]
        if (ref_index_seq not in allowed_indices):
            #if (debug):
            #    print('Warning: checkIndex: first ref_index not allowed', ref_index_seq, file = sys.stderr)
            return(False)
        query_index_start = sp[ref_index_start][0]
        query_index_end = sp[ref_index_end][0]
        query_index_seq = aso.query_sequence[query_index_start:query_index_end]
        if (query_index_seq not in allowed_indices):
            #if (debug):
            #    print('Warning: checkIndex: first query_index not allowed', query_index_seq, file = sys.stderr)
            return(False)
    else:
        ref_index_start = sp.index(sp[:][-i2])
        ref_index_end = sp.index(sp[:][-i1])
        ref_index_seq = aso.get_reference_sequence()[ref_index_start:ref_index_end]
        if (ref_index_seq not in allowed_indices):
            #if (debug):
            #    print('Warning: checkIndex: last ref_index not allowed', ref_index_seq, file = sys.stderr)
            return(False)
        query_index_start = sp[ref_index_start][0]
        query_index_end = sp[ref_index_end][0]
        query_index_seq = aso.query_sequence[query_index_start:query_index_end]
        if (query_index_seq not in allowed_indices):
            #if (debug):
            #    print('Warning: checkIndex: last query_index not allowed', query_index_seq, file = sys.stderr)
            return(False)
    return(True)


def getHammingDistanceBetweenCores(alignedSegmentObject, ref_core_start, ref_core_end):
    """
    Locate the core sequence on the reference and the query.
    Core on ref is located between positions 18 and -18
    Args: alignedSegmentObject, ref_core_start, ref_core_end
    Return: Hamming distance
    """
    aso = alignedSegmentObject
    sp = aso.get_aligned_pairs(matches_only = True)
    ref_start = sp.index(sp[ref_core_start])
    ref_end = sp.index(sp[ref_core_end])
    ref_seq = aso.get_reference_sequence()[ref_start:ref_end]
    query_start = sp[ref_start][0]
    query_end = sp[ref_end][0]
    query_seq = aso.query_sequence[query_start:query_end]
    hamming_distance = hammingDist(ref_seq, query_seq)
    return(hamming_distance)


def doParse(args):
    """
    Parse bam and reference fasta files and count occurences
    of ref seq in bam.
    """
    if (args.mismatches):
        allowed_hd = args.mismatches
    else:
        allowed_hd = mismatches_default

    count_dict = {}
    ref = Fasta(args.fasta, sequence_always_upper = True)
    g = args.genelength
    for k in ref.keys():
       count_dict[k] = 0

    bf = pysam.AlignmentFile(args.bam, 'rb', threads = args.threads)
    for r in bf.fetch():
        if (args.verbose):
            print('\nQuery:', r.query_name, file = sys.stderr)
            print('Ref  :', r.reference_name, file = sys.stderr)
        if (r.query_alignment_length < r.reference_length):
            continue
        first_gene = getGene(r, args.genelength, first = True)
        if (first_gene):
            if (args.verbose):
                print('     : found first gene:', first_gene, file = sys.stderr)
            last_gene = getGene(r, args.genelength, first = False)
            if (last_gene):
                if (args.verbose):
                    print('     : found last gene: ', last_gene, file = sys.stderr)
                if (checkIndex(r, index_pos_one, index_pos_two, allowed_first_indices, first = True)):
                    if (args.verbose):
                        print('     : first query_index_seq: OK', file = sys.stderr)
                    if (checkIndex(r, index_pos_one, index_pos_two, allowed_second_indices, first = False)):
                        if (args.verbose):
                            print('     : last query_index_seq: OK', file = sys.stderr)
                        hd = getHammingDistanceBetweenCores(r, ref_core_start, ref_core_end)
                        if (args.verbose):
                            print('     : Hamming distance:', hd, file = sys.stderr)
                        if (hd <= allowed_hd):
                            count_dict[r.reference_name] += 1
                            if (args.verbose):
                                print('     : add count for:', r.reference_name, file = sys.stderr)

    if (args.output):
        with open(args.output, 'w') as f:
            for k in sorted(count_dict, key = count_dict.get, reverse = True):
                f.write(str(k) + '\t' + ref[k][:].seq + '\t' + str(count_dict[k]) + '\n')
            f.close()
    else:
        for k in sorted(count_dict, key = count_dict.get, reverse = True):
            print(k, '\t', ref[k][:].seq, count_dict[k])

    if (args.verbose):
        print('\n')
        print('Count', '\t', 'Reference', file = sys.stderr)
        for k in sorted(count_dict, key = count_dict.get, reverse = True):
            print(count_dict[k], '\t', k, file = sys.stderr)

    if (args.verbose):
        print('\nEnd of script', file = sys.stderr)


def main():
    parser = argparse.ArgumentParser(
            prog = 'parse',
            description = 'Count references from fasta file in sam/bam file',
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-b', '--bam',
            required = True,
            type = argparse.FileType('r'),
            help = 'BAM file')
    parser.add_argument('-f', '--fasta',
            required = True,
            type = str,
            help = 'fasta (references) file')
    parser.add_argument('-g', '--genelength',
            type = int, nargs = '?', default = gene_length_default,
            help = 'length of leading/trailing gene segment')
    parser.add_argument('-m', '--mismatches',
            type = int, nargs = '?', default = mismatches_default,
            help = 'number of accepted mismatches in core sequence')
    parser.add_argument('-t', '--threads',
            type = int, nargs = '?', default = threads_default,
            help = 'number of threads')
    parser.add_argument('-o', '--output',
            type = str, nargs = '?',
            help = 'output file (default: standard output)')
    parser.add_argument('-v', '--verbose',
            action = 'store_true',
            help = 'increase output verbosity')
    parser.add_argument('-V', '--version',
            action = 'version',
            version = '%(prog)s version ' + __version__)
    #parser.add_argument('--debug',
    #        action = 'store_true', 
    #        help = 'print debug messages to stderr')
    parser.set_defaults(func = doParse)
    args = parser.parse_args()
    args.func(args)

if ( __name__ == "__main__"):
    main()

