citations.add(publications["pyfaidx"])

if config["run_fastp_merge"]:
    readpairs = ["R1", "R2", "R12"]
else:
    readpairs = ["R1", "R2"]

tsv_files = expand(str(COUNTSDIR/"{sample}_{readpair}.tsv"),
    sample = SAMPLES,
    readpair = readpairs)
all_outputs.extend(tsv_files)

parse_bam_config = config["parse_bam"]
extra = parse_bam_config["extra"]
g = parse_bam_config["g"]
m = parse_bam_config["m"]

rule parse_bam:
    input:
        fas = INPUTDIR/"ref.fas",
        bam = BAMDIR/"{sample}_{readpair}.bam",
        bai = BAMDIR/"{sample}_{readpair}.bam.bai"
    output:
        tsv = COUNTSDIR/"{sample}_{readpair}.tsv"
    log:
        LOGDIR/"parse/{sample}_{readpair}.log"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["parse_bam"]["n"] if "parse_bam" in cluster_config else parse_bam_config["n"]
    params:
        extra = extra
    shell:
        """
        python workflow/scripts/parse_bam.py \
            -b {input.bam} \
            -f {input.fas} \
            -g {g} \
            -m {m} \
            -t {threads} \
            -o {output.tsv} \
            &> {log}
        """

