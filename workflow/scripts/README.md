# Scripts folder

- Last modified: fre dec 09, 2022  03:21
- Sign: nylander

## Description

Scripts necessary for the counts procedure.

1. `parse.py`: for parsing the result (bam file) from search. Requires python
   v.3, and python modules pysam and pyfaidx.

## Details for the `parse.py` script

### Description

Parse BAM file and fasta file (reference) and compare the matching sequences
according to the matching schema (see below).

Output is a count of matches for a specific ref seq.

### Usage

    $ parse.py -b BAMFILE -f FASTAFILE [-g [GENELENGTH]] [-m [MISMATCHES]]
        [-t [THREADS]] [-o [OUTPUT]] [-v] [-V] [-h]

### Examples

Allow no mismatches in "core" region:

    $ ./parse.py -b sorted.bam -f ref.fas -o counts.tsv

Allow two mismatches in the core region

    $ ./parse.py -b sorted.bam -f ref.fas -g 9 -m 2 -t 4 -v -o counts.tsv

### Options

    -h,--help              show help message and exit
    -b,--bam BAM           bam file
    -f,--fasta FASTA       fasta (references) file
    -g,--genelength [GENELENGTH] length of leading/trailing gene segment
    -m,--mismatches [MISMATCHES] number of accepted mismatches in core sequence
    -t,--threads [THREADS] number of threads
    -o,--output [OUTPUT] Output file (default: standard output)
    -v,--verbose         increase output verbosity
    -V,--version         show program's version number and exit

### Prerequisites

    $ tab2fasta list.tab > ref.fas
    $ bwa index ref.fas
    $ bwa mem -t 8 ref.fas query.fas | \
        samtools sort -@8 -o my.sorted.bam -
    $ samtools index -@8 my.sorted.bam

### Matching scheme

    GTTCGCCATCTTCAGACTACTAGTTTGATTCGACCCACCGATACCGTTCGTATAATGTATGCTATACGAACGGTATCTCCAGGTCGAATCAAACTAGTTCCCTGAAGATGTCGATG  (Reference sequence)
    GTTCGCCAT (First gene, 9 bp)                                                                               ATGTCGATG  (Second gene, 9 bp)
             CTTCAG  (Non-biological seq, 6 bp)                                                          CTGAAG           (Non-biological seq, 6 bp)
                   ACT (First index, 3 bp. Either ACT or GGA)                                         TCC                 (Second index, 3 bp. Either AGT or TCC (wich are RC of first))
                      ACTAGTTTGATTCGACCCACCGATACCGTTCGTATAATGTATGCTATACGAACGGTATCTCCAGGTCGAATCAAACTAGT                    (Core matching sequence. Allow m mismatches (0-3?))


