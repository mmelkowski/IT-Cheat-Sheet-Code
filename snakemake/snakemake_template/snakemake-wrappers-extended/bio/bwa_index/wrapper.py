__author__ = "Mickael MELKOWSKI"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    'bwa index'
    ' {snakemake.input}'
    ' &> {snakemake.log}'
)
