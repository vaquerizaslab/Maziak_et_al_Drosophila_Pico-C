#boundaries.smk
# ###################
# Insulation at 500 bp
# ###################

rule downsample_insulation_score_500bp:
    input: 
        ancient("fanc/downsampled_ice/{stage}_500bp.hic")
    output:
        "outputs/boundaries/downsampled/insulation_boundaries/{stage}_500bp.insulation"
    resources:
        time= "08:15:00",
        mem= "20G"
    shell:
        "fanc insulation {input} {output} "
        "-w 12kb 16kb 20kb 25kb 30kb"


rule downsample_insulation_score_bed_500bp:
    input: 
        ancient(rules.downsample_insulation_score_500bp.output)
    output:
        "outputs/boundaries/downsampled/insulation_boundaries/{stage}_500bp.insulation_12kb.bed"
    resources:
        time= "08:15:00",
        mem= "10G"
    shell:
        "fanc insulation -o bed {input}"


rule downsample_insulation_score_bw_500bp:
    input:
        done = ancient(rules.downsample_insulation_score_bed_500bp.output), #to make sure it runs one at a time on each sample or else fanc complains
        corr = ancient(rules.downsample_insulation_score_500bp.output)
    output:
        "outputs/boundaries/downsampled/insulation_boundaries/{stage}_500bp.insulation_12kb.bigwig"
    resources:
        time= "08:15:00",
        mem= "10G"
    shell:
        "fanc insulation -o bigwig {input.corr}"


# ###################
# Boundaries at 500 bp
# ###################

rule downsample_boundaries_500bp:
    input: 
        ancient(rules.downsample_insulation_score_500bp.output)
    output:
        temp("outputs/boundaries/downsampled/insulation_boundaries/{stage}_500bp.insulation_boundaries_{window}.bed")
    resources:
        time= "08:15:00",
        mem= "10G"
    shell:
        "fanc boundaries {input} {output} "
        "-w {wildcards.window}"


rule filter_boundaries:
    input:
        ancient(rules.downsample_boundaries_500bp.output)
    output:
        temp("outputs/boundaries/downsampled/insulation_boundaries/filt/{stage}_{filt}_{window}_500bp_boundaries.bed")
    shell:
        "awk '(NR>0) && ($5 > {wildcards.filt} ) ' "
        "{input} "
        "> {output} "


##----------------------------------------------------------##
 ##      remove marginal boundaries  && make into ucsc    ##
##----------------------------------------------------------##

rule remove_marginals_boundaries:
    input:
        boundaries = ancient(rules.filter_boundaries.output),
        marginals = "outputs/marginals/final/dm6_marginals_2kb.bed"
    output:
        temp("outputs/boundaries/downsampled/insulation_boundaries/filt/{stage}_{filt}_{window}_500bp_boundaries.nomarginals.bed")
    shell:
        "bedtools intersect -v "
        "-a {input.boundaries} "
        "-b {input.marginals} "
        "> {output}"


# ###################
# Remove centromeric calls 
# ###################

rule remove_centromeres_boundaries: 
    input: 
        boundaries = ancient("outputs/boundaries/downsampled/insulation_boundaries/filt/{stage}_0.45_16kb_500bp_boundaries.nomarginals.bed"),
        cento = "genome/centromeres_dm6_Zenketal2021.bed"
    output: 
        "outputs/boundaries/final/{stage}_0.45_16kb_500bp_boundaries.nocento.bed"
    shell: 
        "bedtools subtract -A "
        "-a {input.boundaries} "
        "-b {input.cento} "
        " > {output}"