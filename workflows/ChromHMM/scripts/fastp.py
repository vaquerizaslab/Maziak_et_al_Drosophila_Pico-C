"""Script to run 'fastp'."""

import logging
from pathlib import Path
from snakemake.shell import shell


logging.basicConfig(
    format='%(asctime)s %(message)s',
    level=logging.INFO,
    datefmt='%Y/%m/%d %H:%M:%S'
)

threads = snakemake.threads

n = len(snakemake.input)
assert (n == 1 or n == 2), \
    "input must have 1 (single-end) or 2 (paired-end) files."
if n == 1:
    reads = "--in1 {}".format(*snakemake.input)
else:
    reads = "--in1 {} --in2 {}".format(*snakemake.input)

if snakemake.output != []:
    if n == 1:
        trimmed = "--out1 {}".format(*snakemake.output)
    else:
        trimmed = "--out1 {} --out2 {}".format(*snakemake.output)
else:
    trimmed = ""

extras = "--correction --detect_adapter_for_pe " if n == 2 else ""

trimming_options = "-3 --cut_tail_window_size {} --cut_tail_mean_quality {}".format(
    snakemake.params.size, snakemake.params.quality
)

dedup_options = "--dedup --dup_calc_accuracy 4 "

html = "--html {}.html".format(snakemake.params.report_base)
json = "--json {}.json".format(snakemake.params.report_base)

Path(snakemake.params.report_base).parent.mkdir(parents=True)

command = "fastp --version && " \
    "fastp --thread {threads} " \
    "{reads} {trimmed} " \
    "{trimming_options} {dedup_options} " \
    "{json} {html} " \
    "{extras}"

logging.info("Running: " + command.format(**locals()))
shell(command)
