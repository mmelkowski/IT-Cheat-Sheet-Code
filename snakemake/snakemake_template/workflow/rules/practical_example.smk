# this rule use the snakemake wrapper repository online
# the version is set by the config value for all rule update

rule fastp_pe:
    input:
        sample = get_raw_illumina_data
    output:
        trimmed = temp(["results/{sample}/illumina/preprocessing/{sample}_R1_trimmed.fastq.gz",
            "results/{sample}/illumina/preprocessing/{sample}_R2_trimmed.fastq.gz"]),
        html = "results/{sample}/illumina/preprocessing/{sample}_fastp.html",
        json = "results/{sample}/illumina/preprocessing/{sample}_fastp.json"
    log:
        "logs/preprocessing/fastp/{sample}.log"
    params:
        extra = "--qualified_quality_phred 30"
    threads: 4
    group: "illumina"
    wrapper:
        f"{config['wrapper_version']}/bio/fastp"

