import os
import pandas as pd
from snakemake.utils import validate
from pathlib import Path


configfile: "config/config.yaml"

# Threads value determination from workflow.cores
thread_high = workflow.cores
thread_medium = max((workflow.cores * 0.5)//1, 1)
thread_low = max((workflow.cores * 0.25)//1, 1)
thread_one = 1

# Feed your units.tsv file to snakemake with pandas and set an index : <TEMPLATE.CHANGE>
samples = (
    pd.read_csv(
        config["samples"],
        sep="\t",
        dtype={"sample": str, "platform": str},
        comment="#",
    )
    .set_index("sample", drop=False)
    .sort_index()
)

# Feed your units.tsv file to snakemake with pandas and set an index : <TEMPLATE.CHANGE>
units = (
    pd.read_csv(
        config["units"],
        sep="\t",
        dtype={"sample": str},
        comment="#",
    )
    .set_index(["sample"], drop=False)
    .sort_index()
)

# Validate config and units files following a predefined schema
validate(config, schema="../schemas/config.schema.yaml")
validate(units, schema="../schemas/units.schema.yaml")

# safeguard user provided outdir in case a trailing "/" is used
outdir = os.path.join(config["outdir"], "").rstrip("/")
input_dir = os.path.join(config["input_dir"], "").rstrip("/")

# <TEMPLATE.CHANGE> change sample_id with your index (can be a list)
# This will generate your sample list and define your wildcards at the top-most level (target rule input/All_input) or any level that wildcards need to be explicit
def get_unit():
    return units[["sample_id", "ref"]].itertuples()

# This is a basic function to get input files from units.tsv in a programmatic way, function can be even more complex and versatile
def get_fastq_f1(wildcards):
    return os.path.join(input_dir, units.loc[wildcards.sample_id, "fq1"])

def get_fastq_f2(wildcards):
    return os.path.join(input_dir, units.loc[wildcards.sample_id, "fq2"])

def get_ref(wildcards):
    return os.path.join(input_dir, units.loc[wildcards.sample_id, "ref"])

# Other example of basic/input function
def get_raw_illumina_data(wildcards):
    if  pd.notnull(illumina_units.loc[wildcards.sample_id, 'fq2']):
        # paired end local sample
        return [f"{config['local_data']}/{Path(illumina_units.loc[wildcards.sample_id, 'fq1']).name}", 
                f"{config['local_data']}/{Path(illumina_units.loc[wildcards.sample_id, 'fq2']).name}"]
    else:
        # single end local sample
        return [f"{config['local_data']}/{Path(illumina_units.loc[wildcards.sample_id, 'fq1']).name}"]

def get_fastp_output(wildcards):
    if pd.notnull(illumina_units.loc[wildcards.sample_id, 'fq2']):
        return [f"results/{wildcards.sample_id}/illumina/preprocessing/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                f"results/{wildcards.sample_id}/illumina/preprocessing/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    else:
        return [f"results/{wildcards.sample_id}/illumina/preprocessing/{wildcards.sample_id}.fastq.gz"]

# Input function can be called like regular python function to use already written output
def get_chain_input(wildcards):
    return ",".join(get_fastp_output(wildcards))

# input function for a rule with unpack
def get_report_file(wildcards):
    report_dict = {}
    ILL = {
            "fastp" : f"results/{wildcards.sample_id}/illumina/preprocessing/{wildcards.sample_id}_fastp.json"
    }

    ONT = {
            "NanoStats" : f"results/{wildcards.sample_id}/ont/nanoplot/{wildcards.sample}_NanoStats.txt",
    }

    if "ILL" in units:
        report_dict.update(ILL)

    elif "ONT" in units:
        report_dict.update(ONT)

    return report_dict

