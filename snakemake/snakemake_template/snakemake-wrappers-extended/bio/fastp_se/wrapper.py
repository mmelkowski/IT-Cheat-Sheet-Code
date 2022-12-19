__author__ = "Skander Hatira"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "fastp  -i  {snakemake.input.single}  -o {snakemake.output.trimmed}  -w {snakemake.threads}  --html {snakemake.output.html} --json {snakemake.output.json} -D {snakemake.params.extra} {log}"
)
