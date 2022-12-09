citations.add(publications["samtools"])

if config["fastp_merge"]:
    readpairs = ["R1", "R2", "R12"]
else:
    readpairs = ["R1", "R2"]

bai_files = expand(str(BAMDIR/"{sample}_{readpair}.bam.bai"),
    sample = SAMPLES,
    readpair = readpairs)
all_outputs.extend(bai_files)

samtools_config = config["samtools"]
extra = samtools_config["extra"]

rule samtools_index:
    input:
        BAMDIR/"{sample}_{readpair}.bam"
    output:
        BAMDIR/"{sample}_{readpair}.bam.bai"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["samtools"]["n"] if "samtools" in cluster_config else samtools_config["n"]
    shell:
        """
        samtools index -@ {threads} {input}
        """
