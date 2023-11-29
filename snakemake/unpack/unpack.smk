# We can use an input function to return a dict.
# The keyword unpack will then transform the key to input of the rule.

def get_dict_file(wildcards):
    my_dict = {
            "fastp" : f"results/{wildcards.sample}/input_fastp.json",
            "fasta" : f"results/{wildcards.sample}/input.fasta"
    }

rule read_dict:
    input:
        unpack(get_dict_file)
    output:
        "{sample}_analysis.html"
    log:
        log_file = "logs/read_dict/{sample}.log"
    script:
        "../scripts/analysis.py"

# Is equivalent to:
"""
rule read_dict:
    input:
        fastp = f"results/{wildcards.sample}/input_fastp.json",
        fasta = f"results/{wildcards.sample}/input.fasta"
    output:
        "{sample}_analysis.html"
    log:
        log_file = "logs/read_dict/{sample}.log"
    script:
        "../scripts/analysis.py"
"""