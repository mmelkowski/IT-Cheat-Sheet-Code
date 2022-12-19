__author__ = "Skander Hatira"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "metaquast  -r {snakemake.input.ref} -o {snakemake.output} --threads {snakemake.threads} {snakemake.params.extra} {snakemake.input.contigs} {log}"
)
