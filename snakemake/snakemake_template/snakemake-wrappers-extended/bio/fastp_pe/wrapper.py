__author__ = "Skander Hatira"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell(
    "fastp  -i  {snakemake.input.fq1} -I {snakemake.input.fq2} -o {snakemake.output.trimmed_fq1} -O {snakemake.output.trimmed_fq2} -w {snakemake.threads}  --html {snakemake.output.html} --json {snakemake.output.json} -D {snakemake.params.extra} {log}"
)
