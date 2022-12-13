bwa_index_config = config["bwa_index"]

rule bwa_index:
    input:
        INPUTDIR/"ref.fas"
    output:
        idx=multiext("input/ref.fas", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        LOGDIR/"bwa_index/bwa_index.log"
    params:
        algorithm="bwtsw"
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["bwa_index"]["n"] if "bwa_index" in cluster_config else bwa_index_config["n"]
    wrapper:
        "v1.20.0/bio/bwa/index"
