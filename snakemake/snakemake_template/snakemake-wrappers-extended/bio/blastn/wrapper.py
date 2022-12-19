__author__ = "Skander Hatira"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
database = snakemake.params.get("database", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    " blastn -query {snakemake.input} -db {database} {log} -out {snakemake.output}"
)
