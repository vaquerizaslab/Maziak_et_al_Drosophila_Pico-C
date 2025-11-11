#extract_lowmap.smk

# NC9 had the lowest coverage - whi it was used for downsampling.

rule copy_NC9:
    input:
        NC9 = ancient("NC9/hic_ice/NC9_{res}.hic")
    output:
        "fanc/downsampled/NC9_{res}.hic"
    shell:
        "cp "
        "{input.NC9} "
        "/home/nmaziak/mnt/storage/vaquerizas/nmaziak/20240219_microc_analysis/downsampled/"

# ###################
# Downsample
# ###################

rule downsample_new:
    input:
        stage = ancient("{stage}/hic_ice/{stage}_{res}.hic"),
        NC9 = ancient(rules.copy_NC9.output)
    output:
        "fanc/downsampled/{stage}_{res}.hic"
    resources:
        time= "08:15:00",
        mem= "10G",
    retries: 5
    shell:
        "fanc downsample "
        "{input.stage} "
        "{input.NC9} "
        "{output}"

rule bin_hic_downsample:
    input: 
        ancient(rules.downsample_new.output)
    output:
        "fanc/downsampled_ice/{stage}_{res}.hic"
    resources:
        time= "08:15:00",
        mem= "70G"
    threads: 2
    shell:
        "fanc hic -t {threads} -tmp "
        "--low-coverage-auto --normalise --norm-method ice "
        "{input} {output}"

rule copy_NC9_downsample_ice:
    input:
        NC9 = ancient("NC9/hic_ice/NC9_{res}.hic")
    output:
        "fanc/downsampled_ice/NC9_{res}.hic"
    shell:
        "cp "
        "{input.NC9} "
        "fanc/downsampled/"

#----------------------------------# 
# Extract regions of low mappability
#----------------------------------#

rule export_marginals_microc:
    input:
        "fanc/downsampled/hic/{stage}_500bp.hic"
    output:
        temp("fanc/downsampled/marginals/preprocessed/{stage}_500bp_marginals.bed")
    shell: 
        "python scripts/export_marginals.py {input} {output}"

rule extract_marginals:
    input:
        ancient(rules.export_marginals_microc.output)
    output:
        temp("fanc/downsampled/marginals/processed/{stage}_500bp_extracted.bed")
    shell: 
        "awk '(NR>0) && ($5 == 0.0 ) ' "
        "{input} "
        "> {output}"

rule cat_marginals:
    input:
        ancient(expand("fanc/downsampled/marginals/processed/{stage}_500bp_extracted.bed", 
        stage = ["NC9", "NC10", "NC11", "NC12", "NC13", "NC14"]))
    output:
        temp("fanc/downsampled/marginals/final/all_stages_marginals_unmerged.bed")
    shell:
        "cat {input} > "
        "{output}"

rule sort_marginals:
    input: 
        ancient(rules.cat_marginals.output)
    output:
        temp("fanc/downsampled/marginals/final/all_stages_marginals_unmerged.sorted.bed")
    shell:
        "sort -k1,1 -k2,2n {input} > {output} "

rule merge_all:
    input: 
        ancient(rules.sort_marginals.output)
    output:
        temp("fanc/downsampled/marginals/final/all_stages_marginals.bed")
    shell:
        "bedtools merge -d 2100 "
        "-i {input} "
        "> {output}"

rule flank_all:
    input:
        bed = ancient(rules.merge_all.output), 
        genome = "genome/chrom.sizes.ucsc"
    output: 
        "outputs/marginals/final/dm6_marginals_2kb.bed"
    shell: 
        "bedtools slop -i {input.bed} -g {input.genome} -b 2000 "
        "> {output}"




