"""Script to download FASTQ files from ENA."""

import logging
from pathlib import Path
import pandas as pd
from snakemake.shell import shell


logging.basicConfig(
    format='%(asctime)s %(message)s',
    level=logging.INFO,
    datefmt='%Y/%m/%d %H:%M:%S'
)

metadata = pd.read_csv(snakemake.config['metafile'], sep='\t')
run_ids = metadata['run_accession'].tolist()
md5sums = [x.split(';') for x in metadata['fastq_md5'].tolist()]

checksum = {}
if len(md5sums[0]) == 2:
    for run_id, md5 in zip(run_ids, md5sums):
        checksum[f"{run_id}_1"] = md5[0]
        checksum[f"{run_id}_2"] = md5[1]
else:
    for run_id, md5 in zip(run_ids, md5sums):
        checksum[run_id] = md5[0]

row = metadata[metadata['run_accession'] == snakemake.wildcards.id].squeeze()
if ';' in row.get('fastq_ftp'):
    assert 'run' in snakemake.wildcards.keys(), \
        "output in single-end format, but metadata is for paired-end data."
    download_link = row.get('fastq_ftp').split(';')[int(snakemake.wildcards.run) - 1]
    sample_name = f"{snakemake.wildcards.id}_{snakemake.wildcards.run}"
else:
    assert 'run' not in snakemake.wildcards.keys(), \
        "output in paired-end format, but metadata is for single-end data."
    download_link = row.get('fastq_ftp')
    sample_name = f"{snakemake.wildcards.id}"

dl_command = f"wget --retry-connrefused -O {snakemake.output[0]} ftp://{download_link}"
md5_command = f"echo '{checksum[sample_name]}  {snakemake.output[0]}' | md5sum -c"

logging.info("Running: " + dl_command)
shell(dl_command)
logging.info("Running: " + md5_command)
shell(md5_command)
shell(f"chmod 444 {snakemake.output[0]}")
