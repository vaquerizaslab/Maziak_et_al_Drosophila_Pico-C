#pairtools.smk

#-------------------#
# Alignment 
#-------------------#

rule alignment:
    input:
        R1 = "fastq/{id}_1.fastq.gz",
        R2 = "fastq/{id}_2.fastq.gz"
    output:
        temp("outputs/microc_bams/{id}.sam")
    threads: 8
    resources:
        time= "28:15:00",
        mem= "10G"
    shell:
        "bwa mem -t {threads} -T0 -5SP "
        "genome/BWAdm6Index_ucsc/genome "
        "{input.R1} {input.R2} "
        "-o {output}" 


#-------------------#
# Parsing
#-------------------#

rule parse: 
    input:
        rules.alignment.output
    output:
        temp("pairtools/{id}.parsed.pairsam")
    threads: 4
    resources:
        time= "48:15:00",
        mem= "50G"
    shell:
        "pairtools parse2 "
        "--min-mapq 3 "
        "--report-position outer "
        "--max-inter-align-gap 30 "
        "--chroms-path chrom.sizes.ucsc "
        "{input} > {output}"

rule split: 
    input:
        rules.parse.output
    output:
        temp("pairtools/split/{id}.parsed.pairs")
    resources:
        time= "48:15:00",
        mem= "50G"
    shell:
        "pairtools split "
        "--output-pairs {output} "
        "{input}"

rule sort:
    input:
        rules.split.output
    output:
        "pairtools/pairs/{id}.sorted.pairs"
    resources:
        time= "48:15:00",
        mem= "40G"
    threads: 3
    shell:
        "pairtools sort --nproc {threads} "
        " --tmpdir=./ "
        "{input} > {output}"

## For QC - dedup in paper was done with FAN-C ##
rule pt_dedup: 
    input: 
        rules.sort.output
    output: 
        pairs = "pairtools/pairs/{id}.dedup.pairs",
        stats = "pairtools/qc/pairtools_stats/{id}.stats"
    resources:
        time= "48:15:00",
        mem= "40G"    
    shell: 
        "pairtools dedup "
        "--max-mismatch 1 "
        "{input} "
        "-o {output.pairs} "
        "--output-stats {output.stats}"
