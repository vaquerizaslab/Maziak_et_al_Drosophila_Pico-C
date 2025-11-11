"""Script to run 'bwa-mem2 mem'."""

import logging
from snakemake.shell import shell

logging.basicConfig(
    format='%(asctime)s %(message)s',
    level=logging.INFO,
    datefmt='%Y/%m/%d %H:%M:%S'
)

threads = snakemake.threads
index = '-x ' + snakemake.params.index
output = '-S ' + snakemake.output[0]

n = len(snakemake.input)
assert (n == 1 or n == 2), \
    "input must have 1 (single-end) or 2 (paired-end) files."
reads = ' '.join(f'-{i} {{}}' for i in range(1, n+1)).format(*snakemake.input)

alignment_options = "--maxins 2000 -N 1 "

output_options = "--no-unal --no-mixed --no-discordant "

command = "bowtie2 --threads {threads} --time " \
    "-U {reads} {index} {output} " \
    "{alignment_options}" \
    "{output_options}"

logging.info("Running: " + command.format(**locals()))
shell(command)
