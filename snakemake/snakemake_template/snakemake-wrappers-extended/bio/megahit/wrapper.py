__author__ = "Skander Hatira"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
prefix = snakemake.params.get("prefix", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
minlen = snakemake.params.get("minlen", "")
inp = (
    f"-r {snakemake.input.single}"
    if len(snakemake.input) == 1
    else f"-1 {snakemake.input.fq1} -2 {snakemake.input.fq2}",
)

shell(
    "rm -r {snakemake.output.out_dir}; megahit  {inp} -t {snakemake.threads}  --min-contig-len {minlen} -o {snakemake.output.out_dir} --out-prefix {prefix}   --continue  {log}"
)
