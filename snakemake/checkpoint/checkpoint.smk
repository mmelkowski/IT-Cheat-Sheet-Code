# The rule will produced a folder with multuiple output.
# since the output is not already known the checkpoint is used instead of rule the re evaluate the dag
checkpoint seqkit_split:
    input:
        fasta="input.fasta",
    output:
        outdir=directory("results/input.fasta.split"),
    log:
        "logs/seqkit/split/input.log",
    params:
        command="split",
        extra="-i --by-id-prefix ''",
    threads: 1
    wrapper:
        "v3.0.0/bio/seqkit"

rule seqkit_seq_upper:
    input:
        fasta="results/input.fasta.split/{i}.fasta"
    output:
        fasta="results/input.fasta.split/{i}.upper.fasta"
    log:
        "logs/seqkit/seq/{sample}_{i}_seqkit_seq_upper.log"
    params:
        command="seq",
        extra="--upper-case",
    threads: 1
    wrapper:
        "v3.0.0/bio/seqkit"

# function for merging
def collect_file(wildcards):
    nodes = checkpoints.seqkit_split.get(**wildcards).output.outdir
    files = expand("results/input.fasta.split/{i}.upper.fasta",
                   i=glob_wildcards(os.path.join(nodes, '{i}.upper.fasta')).i
                   )
    return files

# Concatenate the files for one results without the new wildcard introduced
rule merge_files:
    input:
        collect_file
    output:
        "results/seq.fasta"
    log:
        "logs/merge_files.log"
    threads: 1
    shell:
        "(cat {input} > {output}) > {log} 2>&1"