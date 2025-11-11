#rules/pt2fanc.smk

#-------------------#
# Set-up fragments 
#-------------------#

rule bin_genome:
    output:
        "genome/{genome}_{fragments}_fragments.bed"
    resources:
        time= "00:15:00",
        mem= "1G"
    params:
        chrs = lambda wildcards: ",".join(chromosomes[wildcards.genome]),
        genome_fa = touch(lambda wildcards: config["genome_fasta"][wildcards.genome])
    shell:
        "fanc fragments -c {params.chrs} {params.genome_fa} {wildcards.fragments} {output}" 


#-------------------#
# Make .hic
#-------------------#

rule pt2fanc:
    input:
        ancient(rules.sort.output)
    output:
        "fanc/pairs/{id}.fanc.pairs" 
    resources:
        time= "08:15:00",
        mem= "30G"
    shell:
        "fanc pairs {input} {output} "
        "-g dm_100_fragments.bed"

rule dedup:
    input:
        ancient(rules.pt2fanc.output)
    output:
        "fanc/stats/{id}.txt"
    threads: 2
    resources:
        time= "08:15:00",
        mem= "50G"
    shell:
        "fanc --version && "
        "fanc pairs "
        "--filter-pcr-duplicates 1 "
        "--statistics {output} "
        "{input}"


#-------------------#
# Make .hic 
#-------------------#

rule make_picoc:
    input:
        pairs = ancient(rules.pt2fanc.output),
        stats = ancient(rules.dedup.output)
    output:
        "/fanc/hic/{id}.hic" 
    resources:
        time= "08:15:00",
        mem= "40G"
    shell:
        "fanc hic -tmp {input.pairs} {output}"

## Used for PCA ##
rule bin_hic_ice:
    input: 
        ancient(rules.make_picoc.output)
    output:
        "fanc/hic/{id}_{res}.hic"
    resources:
        time= "08:15:00",
        mem= "40G"
    threads: 2
    shell:
        "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto --normalise --norm-method ice "
        "{input} {output}"

#-------------------#
# Merge
#-------------------#

rule merge_hic:
    input:
        ancient(lambda wc: expand(
            rules.make_picoc.output,
            id = stage2sample[wc.stage]))
    output:
        "{stage}/{stage}.hic"
    resources:
        time= "08:15:00",
        mem= "50G"
    threads: 4
    shell:
        "fanc hic -t {threads} {input} "
        "{output}"

rule bin_hic_merge:
    input: 
        ancient(rules.merge_hic.output)
    output:
        "{stage}/hic_ice/{stage}_{res}.hic"
    resources:
        time= "08:15:00",
        mem= "70G"
    threads: 2
    shell:
        "fanc hic -t {threads} -tmp "
        "-b {wildcards.res} --low-coverage-auto --normalise --norm-method ice "
        "{input} {output}"