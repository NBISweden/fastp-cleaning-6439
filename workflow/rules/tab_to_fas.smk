citations.add(publications["seqkit"])

ref_fas_file = INPUTDIR/"ref.fas"
all_outputs.append(str(ref_fas_file))

rule tab_to_fasta:
    input:
        INPUTDIR/"list.tab"
    output:
        INPUTDIR/"ref.fas",
    conda:
        "../envs/fastp-cleaning.yaml"
    shell:
        "seqkit tab2fx {input} > {output}"

