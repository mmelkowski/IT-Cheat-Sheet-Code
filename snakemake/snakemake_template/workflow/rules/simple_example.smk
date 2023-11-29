rule fastp_paired:
    input:
        fq1=get_fastq_f1,
        fq2=get_fastq_f2,
    output:
        trimmed_fq1="{outdir}/{sample_id}/{sample_id}_out.R1.fastq.gz",
        trimmed_fq2="{outdir}/{sample_id}/{sample_id}_out.R2.fastq.gz",
        html="{outdir}/{sample_id}/{sample_id}_fastp.html",
        json="{outdir}/{sample_id}/{sample_id}_fastp.json",
    log:
        '{outdir}/{sample_id}/logs/fastp_paired_log__{sample_id}.log'
    params:
        extra="-l 60"
    threads:
        thread_medium
    wrapper:
        "file:snakemake-wrappers-extended/bio/fastp_pe"


rule copy_ref:
    input:
        ".tests/data/{ref}.fasta"
    output:
        "{outdir}/ref/{ref}.fasta"
    log:
        "{outdir}/ref/logs/copy__{ref}.out"
    shell:
        "cp {input} {output}"


rule bwa_index:
    input:
        "{outdir}/ref/{ref}.fasta"
    output:
        "{outdir}/ref/{ref}.fasta.amb",
        "{outdir}/ref/{ref}.fasta.ann",
        "{outdir}/ref/{ref}.fasta.bwt",
        "{outdir}/ref/{ref}.fasta.pac",
        "{outdir}/ref/{ref}.fasta.sa"
    log:
        '{outdir}/ref/logs/bwa_index_log__{ref}.out'
    wrapper:
        "file:snakemake-wrappers-extended/bio/bwa_index"


rule bwa_mem_paired:
    input:
        ref_index=rules.bwa_index.output,
        ref=rules.bwa_index.input,
        fq1="{outdir}/{sample_id}/{sample_id}_out.R1.fastq.gz",
        fq2="{outdir}/{sample_id}/{sample_id}_out.R2.fastq.gz",
    output:
        "{outdir}/{sample_id}/{sample_id}__{ref}.sam"
    log:
        "{outdir}/{sample_id}/logs/bwa_mem_paired_err__{sample_id}__{ref}.log"
    threads:
        1
    params:
        ""
    wrapper:
        "file:snakemake-wrappers-extended/bio/bwa_mem_pe"
