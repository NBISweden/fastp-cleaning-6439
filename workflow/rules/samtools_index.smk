citations.add(publications["samtools"])

if config["run_fastp_merge"]:
    readpairs = ["R1", "R2", "R12"]
else:
    readpairs = ["R1", "R2"]

bai_files = expand(str(BAMDIR/"{sample}_{readpair}.bam.bai"),
    sample = SAMPLES,
    readpair = readpairs)
all_outputs.extend(bai_files)

samtools_index_config = config["samtools_index"]
extra = samtools_index_config["extra"]

rule samtools_index:
    input:
        BAMDIR/"{sample}_{readpair}.bam"
    output:
        BAMDIR/"{sample}_{readpair}.bam.bai"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["samtools_index"]["n"] if "samtools_index" in cluster_config else samtools_index_config["n"]
    shell:
        """
        samtools index -@ {threads} {input}
        """
