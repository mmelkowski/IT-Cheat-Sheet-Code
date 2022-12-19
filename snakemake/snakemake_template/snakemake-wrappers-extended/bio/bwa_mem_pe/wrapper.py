__author__ = "Mickael MELKOWSKI"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    'bwa mem'
    ' {snakemake.params}'
    ' -t {snakemake.threads}'
    ' {snakemake.input.ref}'
    ' {snakemake.input.fq1}'
    ' {snakemake.input.fq2}'
    ' > {snakemake.output}'
    ' 2> {snakemake.log}'
)
