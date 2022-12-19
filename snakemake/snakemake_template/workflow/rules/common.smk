import os
import pandas as pd
from snakemake.utils import validate


configfile: "config/config.yaml"


# Feed your units.tsv file to snakemake with pandas and set an index : <TEMPLATE.CHANGE>
units = pd.read_csv(config["units"], dtype=str, sep="\t").set_index(
    ["sample_id"], drop=False
)

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

# Validate config and units files following a predefined schema
validate(config, schema="../schemas/config.schema.yaml")
validate(units, schema="../schemas/units.schema.yaml")
