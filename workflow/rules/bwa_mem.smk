citations.add(publications["bwa"])

if config["fastp_merge"]:
    readpairs=["R1", "R2", "R12"]
else:
    readpairs=["R1", "R2"]

bam_files = expand(str(BAMDIR/"{sample}_{readpair}.bam"),
    sample=SAMPLES,
    readpair=readpairs)
all_outputs.extend(bam_files)

bwa_config = config["bwa"]
extra = bwa_config["extra"]

rule bwa_mem:
    input:
        ref = INPUTDIR/"ref.fas",
        idx = multiext("input/ref.fas", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        fqgz = FASTQDIR/"{sample}_{readpair}.fq.gz"
    output:
        bam = BAMDIR/"{sample}_{readpair}.bam"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["bwa"]["n"] if "bwa" in cluster_config else bwa_config["n"]
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.fqgz} | \
              samtools sort -@{threads} -o {output.bam} -
        """
