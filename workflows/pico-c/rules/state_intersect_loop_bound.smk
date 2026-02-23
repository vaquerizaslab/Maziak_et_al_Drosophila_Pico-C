# state_intersect_loop_bound.smk

##-------------------------------------------##
    #Process Loops#
##-------------------------------------------##

rule make_long_subset: 
    input: 
        "mustache/final/0.01/{stage}/{stage}_250_1kb_2kb_4kb_16kb.clean.bed",
    output: 
        "loops_bed/{stage}_subset_long.bed"
    shell: 
        "python3 scripts/make_bed.py {input} | "
        "sort -k1,1 -k2,2n "
        "> {output}"

##-------------------------------------------##
    #make local regions slop#
##-------------------------------------------##

rule slop_loops_subset:
    input:
        bed = rules.make_long_subset.output,
        cento_marg = "cento_marg_dm6.bed"
    output: 
        "loops_bed/slop/{stage}_subset_{slop}.bed"
    shell: 
        "bedtools slop "
        "-b {wildcards.slop} "
        "-g chrom.sizes.ucsc "
        "-i {input.bed} | "
        "sort -k1,1 -k2,2n | "
        "bedtools subtract -a stdin -b {input.cento_marg} "
        "> {output}"

rule slop_boundaries_subset:
    input:
        bed = "downsampled/insulation_boundaries/filt/no_cento/{stage}_0.45_16kb_500bp_boundaries.nocento.bed",
        cento_marg = "cento_marg_dm6.bed"
    output: 
        "filt/no_cento/{stage}_0.45_16kb_500bp_boundaries.nomarginals_{slop}.bed" 
    shell: 
        "bedtools slop "
        "-b {wildcards.slop} "
        "-g chrom.sizes.ucsc "
        "-i {input.bed} | "
        "sort -k1,1 -k2,2n | "
        #"bedtools merge | "
        "bedtools subtract -a stdin -b {input.cento_marg} "
        "> {output}"

##-------------------------------------------##
    #simplify bedfile for states#
##-------------------------------------------##

rule add_state:
    input:
        "ChromHMM/beds/{id}_20_altcolstate.bed"
    output:
        "ChromHMM/beds/simple/{id}.bed"
    shell:
        "awk -v OFS='\t' 'NR > 0 {{$4=\"state\"$4; print}}'  {input} | "
        "awk '{{print $1, $2, $3, $4}}' | "
        "awk -F'[[:blank:]]' -v OFS=\"\t\" '{{$1=$1; print}}' "
        "> {output}"


##-------------------------------------------##
    #overlap slop state enrichment#
##-------------------------------------------##

rule intersect_states_loops_subset_slop:
    input:
        loops = rules.slop_loops_subset.output,
        states = "ChromHMM/beds/simple/{id}.bed",
    output:
        "state_intersect/loops/slop/{stage}_subset_long.{id}.{slop}.bed"
    shell: 
        "bedtools intersect -wo "
        "-a {input.loops} "
        "-b {input.states} "
        "> {output} "


rule intersect_states_boundaries_subset_slop:
    input:
        loops = rules.slop_boundaries_subset.output,
        states = "ChromHMM/beds/simple/{id}.bed",
    output:
        "state_intersect/boundaries/slop/{stage}_0.45_16kb_500bp_boundaries.nomarginals.{id}.{slop}.bed"
    shell: 
        "bedtools intersect -wo "
        "-a {input.loops} "
        "-b {input.states} "
        "> {output} "


##-------------------------------------------##
    #overlap state enrichment#
##-------------------------------------------##

rule intersect_states_boundaries:
    input:
        boundaries = "downsampled/insulation_boundaries/filt/no_cento/{stage}_0.45_16kb_500bp_boundaries.nocento.bed",
        states = "ChromHMM/beds/simple/{id}.bed"
    output:
        "state_intersect/boundaries/{stage}_0.45_16kb_500bp_boundaries.nomarginals.{id}.bed"
    shell: 
        "bedtools intersect -wo "
        "-a {input.boundaries} "
        "-b {input.states} "
        "> {output} "


rule intersect_states_loops_subset:
    input:
        loops = "loops_bed/{stage}_subset_long.bed",
        states = "ChromHMM/beds/simple/{id}.bed"
    output:
        "state_intersect/loops/{stage}_subset_long.{id}.bed"
    shell: 
        "bedtools intersect -wo "
        "-a {input.loops} "
        "-b {input.states} "
        "> {output} "


