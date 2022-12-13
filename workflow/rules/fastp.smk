citations.add(publications["fastp"])

trimmed_qc = expand(str(FASTQDIR/"{sample}_{readpair}.fq.gz"),
    sample = SAMPLES,
    readpair = ["R1", "R2"])
all_outputs.extend(trimmed_qc)

fastp_config = config["fastp"]
if config["run_fastp_dedup"]:
    extra = fastp_config["extra"] + " --length_required " + str(fastp_config["length_required"]) + " --dedup --dup_calc_accuracy " + str(fastp_config["dup_calc_accuracy"])
else:
    extra = fastp_config["extra"] + " --length_required " + str(fastp_config["length_required"])

rule fastp:
    input:
        read1 = INPUTDIR/config["fn_pattern"].format(sample="{sample}", readpair="R1"),
        read2 = INPUTDIR/config["fn_pattern"].format(sample="{sample}", readpair="R2")
    output:
        read1 = FASTQDIR/"{sample}_R1.fq.gz",
        read2 = FASTQDIR/"{sample}_R2.fq.gz",
        html = FASTPDIR/"{sample}.fastp.html"
    log:
        LOGDIR/"fastp/{sample}.fastp.log"
    shadow:
        "shallow"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["fastp"]["n"] if "fastp" in cluster_config else fastp_config["n"]
    params:
        extra = extra
    shell:
        """
        fastp \
            --in1 {input.read1} \
            --in2 {input.read2} \
            --out1 {output.read1} \
            --out2 {output.read2} \
            --html {output.html} \
            --thread {threads} \
            {params.extra} \
            &> {log}
        """

