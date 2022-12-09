rule bwa_index:
    input:
        INPUTDIR/"ref.fas",
    output:
        idx=multiext("input/ref.fas", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index/ref.log",
    params:
        algorithm="bwtsw",
    conda:
        "../envs/fastp-cleaning.yaml"
    wrapper:
        "v1.20.0/bio/bwa/index"
