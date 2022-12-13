citations.add(publications["bwa"])

if config["run_fastp_merge"]:
    readpairs=["R1", "R2", "R12"]
else:
    readpairs=["R1", "R2"]

bam_files = expand(str(BAMDIR/"{sample}_{readpair}.bam"),
    sample=SAMPLES,
    readpair=readpairs)
all_outputs.extend(bam_files)

bwa_mem_config = config["bwa_mem"]

rule bwa_mem:
    input:
        ref = INPUTDIR/"ref.fas",
        idx = multiext("input/ref.fas", ".amb", ".ann", ".bwt", ".pac", ".sa"),
        fqgz = FASTQDIR/"{sample}_{readpair}.fq.gz"
    output:
        bam = BAMDIR/"{sample}_{readpair}.bam"
    log:
         LOGDIR/"bwa_mem/{sample}_{readpair}.log"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["bwa_mem"]["n"] if "bwa_mem" in cluster_config else bwa_mem_config["n"]
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.fqgz} | \
              samtools sort -@{threads} -o {output.bam} - &> {log}
        """
