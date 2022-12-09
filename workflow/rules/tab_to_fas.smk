citations.add(publications["seqkit"])

ref_fas_file = INPUTDIR/"ref.fas"
all_outputs.append(str(ref_fas_file))

tab_to_fasta_config = config["tab_to_fasta"]

rule tab_to_fasta:
    input:
        INPUTDIR/"list.tab"
    output:
        INPUTDIR/"ref.fas",
    conda:
        "../envs/fastp-cleaning.yaml"
    threads:
        cluster_config["tab_to_fasta"]["n"] if "tab_to_fasta" in cluster_config else tab_to_fasta_config["n"]
    shell:
        "seqkit tab2fx --threads {threads} {input} > {output}"

