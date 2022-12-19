[![Snakemake](https://img.shields.io/badge/snakemake-≥7.14.0-brightgreen.svg)](https://snakemake.github.io)

# The Wrapper Repository

The Wrapper Repository is a collection of reusable wrappers extending [snakemake's wrapper repository](https://github.com/snakemake/snakemake-wrappers) with additional in-house wrappers.

# Contributing

To extend this repository with a new wrapper we need to provide:

* `environment.yaml` to install dependencies using conda
* `meta.yaml` a metadata file contining info about input/output, special parameters and anything we deem necessary as well as a reference to dependencies
* `wrapper.py` the entrypoint that can deal with arbitrary `input` `output`
* `test` a folder containing a Snakefile and some dummy data that will be used for testing
The resulting folder structure for metaquast rule is as follows:

```bash
bio
└── metaquast
    ├── environment.yaml
    ├── meta.yaml
    ├── test
    │   ├── Snakefile
    │   ├── contigs1.fasta
    │   ├── contigs2.fasta
    │   └── genome.fasta
    └── wrapper.py

2 directories, 7 files
```

## Unit tests

This repository comes with a predefined `test.py` with some boilerplate functions.

If all the files mentioned above are present, we can directly add this snippet:

```python
@skip_if_not_modified
def test_metaquast():
    run(
        "bio/metaquast",
        ["snakemake", "--cores", "1", "metaquast_out", "--use-conda", "-F"]
    )
```

In this example `metaquast_out` stands for the target of that rule, in this case a directory.

If you're using conda and have pytest installed you can run:

```bash
pytest test.py
```

Alternatively, you can run:

```bash
python -m pytest test.py
```

## Why wrappers?

Wrappers are the best way to control the `functional` part of a rule without taking away the freedom to adapt it to different use cases

### Example

```python
rule metaquast:
    input:
        contigs=["contigs1.fasta", "contigs2.fasta"],
        ref="genome.fasta",
    output:
        directory("metaquast_out"),
    log:
        "logs/metaquast.log",
    params:
        extra="--min-contig 5 --min-identity 95.0",
    wrapper:
        "main/bio/metaquast"
```

In this rule the functional part is the `wrapper` directive, by calling it we ensure the command and software versions are always the same, while being able to use `input functions`,`wildcards` and `params`...etc.

We can therefore use the samme wrapper for :

* `workflow_A` that will assess the quality of 2 assemblies.
  * `contigs1.fasta`
  * `contigs2.fasta`
* `workflow_B` that will assess the quality depending on reference type using an input function.

```python
rule metaquast:
    input:
        unpack(metaquast_input)
    output:
        "metaquast_out_{reference}"
    .
    .
    wrapper:
        "main/bio/metaquast"
```

with `metaquast_input` :

```python
def metaquast_input(wildcards):
    if wildcards.reference == "human":
        return {
            "contigs": ["contigs1.fasta", "contigs2.fasta"],
            "ref": "genome.fasta",
        }
    return {
        "contigs": ["contigs1.fasta", "contigs2.fasta"],
        "ref": config["standard_reference"],
    }

```

NOTE: in this kind of setup we can use all snakemake features e.g:

```python
    return {
        "contigs" : rules.some_rule.output["contigs"],
        "ref": config["standard_reference"],
    }
```

All the while knowing that the underlying command and software versions used are the same.

Visit [the snakemake wrapper documentation](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html) for more information about contributing.
