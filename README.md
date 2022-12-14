# fastp-cleaning-6439

- Last modified: mån dec 19, 2022  04:28
- Sign: nylander

## Description

Workflow using snakemake for filtering fastq files using
[fastp](https://github.com/OpenGene/fastp), followed by mapping with [bwa
mem](https://github.com/lh3/bwa) and counting of occurrences of reference
sequences.

### Input

1. gzipped paired-end Illumina `.fastq.gz` files.
2. `list.tab` file with reference sequences.

### Output

- filtered `.fq.gz` files placed in a new output directory
- filtering reports (before/after)
- merged pair-end reads `.fq.gz` files (optional)
- fasta files from filtered fastq (optional)
- mapping results (bam files)
- counts of reference sequences, one for each filtered fastq file (tsv files,
  in folder output/counts)
- workflow report in html format (in folder output)

### Filtering steps

1. Filter on [quality](https://github.com/OpenGene/fastp#quality-filter)
2. Filter on [length](https://github.com/OpenGene/fastp#length-filter)
3. Per read cutting per quality using a [sliding window from front to
   tail](https://github.com/OpenGene/fastp#per-read-cutting-by-quality-score)
4. Adapters are automatically [detected and
   trimmed](https://github.com/OpenGene/fastp#adapters)
5. PCR-deduplicaiton (optional, by editing the
   [config.yaml](config/config.yaml) file)
6. Merging of paired-end read pairs (optional, by editing the
   [config.yaml](config/config.yaml) file)
7. Convert the filtered fastq files to fasta (optional, by editing the
   [config.yaml](config/config.yaml) file)

Please see the [fastp-manual](https://github.com/OpenGene/fastp/wiki) for
details on the filtering procedures.

**Note:** for the current version of the workflow, a minimum sequence
length of 116 is applied. This would correspond to the length of the 
sequences in the example file `list.tab` that was provided. In order to
change this length, edit the file `config/config.yaml`.

### Searching/Mapping steps

Mapping is done with [bwa mem](https://github.com/lh3/bwa) using the fastq
files as queries, and the sequences in `list.tab` as the "genome".

### Parsing the search results

Searching for reference sequences in the fastq files are done according to a
specific schema using a script,
[`parse_bam.py`](workflow/scripts/parse_bam.py), which parses the bam files.
Please see the [README file](workflow/scripts/README.md) in the
workflow/scripts directory for details. Workflow defaults for the script
`parse_bam.py` are specified in the file [config.yaml](config/config.yaml).

---

## How to run locally with conda

1. Install
    - [`snakemake`](https://snakemake.readthedocs.io/en/stable/#)
    - [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Clone the repository: `git clone
https://github.com/nylander/fastp-cleaning-6439.git`
3. Put input data (gzip-compressed paired-end illumina fastq files) in folder
`fastp-cleaning-6439/input`
4. Put the file `list.tab` in `fastp-cleaning-6439/input`
5. Review the `fastp-cleaning-6439/config/config.yaml` and make sure input file name
endings (currently `_R1_001.fastq.gz`), matches your input files, select the
steps used by the pipeline, and change options for software used if needed.
6. Make sure your current working directory is `fastp-cleaning`
7. Test run `snakemake --use-conda -n`
8. Run with `snakemake --use-conda --cores N` (substitute N with the number of
cores you wish to use)

## How to run on UPPMAX

Note: On rackham we are loading conda and snakemake (v.5) as modules. This can
be done manually or using a script
([rackham/scripts/init.sh](rackham/scripts/init.sh)). For convenience, we also
start the run with the [Makefile](Makefile).

1. Log in to [UPPMAX](https://uppmax.uu.se/) (rackham.uppmax.uu.se)
2. Clone the repository in a project folder: `git clone
   https://github.com/nylander/fastp-cleaning-6439.git`
3. Edit the file `fastp-cleaning-6439/rackham/rackham.yaml` to add your CPU-project
   account number. For example: `sed -i -e 's/snic1234-56-789/snic2022-01-001/'
   rackham/rackham.yaml`.
4. Add input files (use symbolic links to save space) to
   `fastp-cleaning-6439/input/`. Don't forget the `list.tab`.
5. Review the `fastp-cleaning-6439/config/config.yaml` and make sure input file-name
   endings (currently `_R1_001.fastq.gz`), matches your input files, select the
   steps used by the pipeline, and change options for software used if needed.
6. Start a screen session: `screen -S fast-cleaning`
7. Load modules: `source rackham/scripts/init.sh`
8. Test run: `make slurm-test`
9. Run: `make slurm-run`
10. Detach from the screen session (Ctrl+A, Ctrl+D).
11. Monitor progress by, e.g., `jobinfo -u $USER`
12. To re-attach to screen session, use `screen -R fast-cleaning`

## Acknowledgements

The pipeline is an extension of
[fastp-cleaning](https://github.com/nylander/fastp-cleaning), heavily
influenced by [stag-mwc](https://github.com/ctmrbio/stag-mwc).  Thanks to
[Marcel Martin](https://github.com/marcelm) and [John
Sundh](https://github.com/johnne) for important feedback on the bam parsing
(MM) and the workflow (JS).


## License and copyright

Copyright (c) 2021, 2022 Johan Nylander

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

