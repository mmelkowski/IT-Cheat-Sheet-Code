> TO DELETE: This template is a quick folder for a conda based pipeline.
> This is based on the template made by Skander Haitiri at PathoQuest based itself on the snakemake guideline.

# Snakemake workflow: `<name>`

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.14.0-brightgreen.svg)](https://snakemake.github.io)

This repository aims to be a template to use for new snakemake pipelines, ideally we collectively keep it up to date with snakemake's latest recommandations.
It is based on snakemake's [template](https://github.com/snakemake-workflows/snakemake-workflow-template) and [distribution and reproducibility guidelines](https://snakemake.readthedocs.io/en/v7.14.2/snakefiles/deployment.html#distribution-and-reproducibility).

A Snakemake workflow for `<description>`

## Template usage

>## TO DO BEFORE USE
>
>* replace all `<description>`, `<name>` , `<tag>` placehoders
>
---

## Workflow building

This repository comes with wrappers submodules
to add a rule using a local version of a wrapper:

```python
rule metaquast:
    input:
        contigs=[".tests/data/contigs1.fasta", ".tests/data/contigs2.fasta"],
        ref=".tests/data/genome.fasta",
    output:
        directory("results/metaquast_out"),
    log:
        "results/logs/metaquast.log",
    params:
        extra="--min-contig 5 --min-identity 95.0",
    wrapper:
        "file:snakemake-wrappers-extended/bio/metaquast"
        # use file: to use local version or https:// for remote
```

## Workflow usage

### Local mode

* Create environement

```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

* Run workflow

```bash
snakemake --profile config/profile/test
```

### Docker mode

* Create environment

```bash
docker build -t <name> .
```

* From docker image

```bash
docker run --rm --mount type=bind,src=$PWD/config/,dst=/pipeline/config --mount type=bind,src=$PWD/results/,dst=/pipeline/results <name> --profile config/profiles/test
```

* Or from inside a docker container

```bash
docker run --rm -it --mount type=bind,src=$PWD/config/,dst=/pipeline/config --mount type=bind,src=$PWD/results/,dst=/pipeline/results --entrypoint '/bin/bash' <name>

micromamba run -n snakemake snakemake --profile config/profiles/test
```
