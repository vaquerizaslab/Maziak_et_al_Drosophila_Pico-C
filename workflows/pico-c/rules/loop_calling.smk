#loop_calling.smk"


# ###################
# Loop Calling
# ###################

#Remove centromeric regions

rule subset_hic:
    input:
        "{stage}/{stage}.hic"
    output:
        "outputs/subset/hic/{stage}_subset.hic"
    resources:
        time= "08:15:00",
        mem= "10G"
    shell:
        "fanc hic "
        "--subset chr2L:1-22000000,chr2R:6000000-25286936,chr3L:1-23000000,chr3R:4100000-32079331,chrX:1-21500000 "
        "{input} "
        "{output}"

rule bin_hic_subset:
    input: 
        ancient(rules.subset_hic.output)
    output:
        "outputs/subset/hic/binned/{stage}_125bp.hic"
    resources:
        time= "08:15:00",
        mem= "50G"
    threads: 2
    shell:
        "fanc hic -t {threads} -tmp "
        "-b 125bp --low-coverage-auto --ice-correct "
        "{input} {output}"

rule make_cool_subset:
    input:
        ancient(rules.bin_hic_subset.output)
    output:
        temp("outputs/subset/cool/{stage}_125bp.cool")
    resources:
        time= "88:15:00",
        mem= "40G"
    threads: 2
    shell:
        "fanc to-cooler "
        "-t {threads} "
        "--uncorrected "
        "--no-multi "
        "{input} {output}" 

rule cooler_zoomify_subset:
    input:
        rules.make_cool_subset.output
    output: 
        "outputs/cool/subset/{stage}_125bp_subset_madmaxiter0.mcool"
    resources:
        time= "08:15:00",
        mem= "20G"
    threads: 8 
    shell:
        "cooler zoomify "
        "--nproc {threads} "
        "--balance "
        "--balance-args "
        "'--mad-max 0 --max-iters 1000' "
        "{input} "
        "-o {output}"



# ###################
# Loop Calling
# ###################

rule mustache_lowres: 
    input:
        rules.cooler_zoomify_subset.output
    output:
        "mustache/lowres/{stage}_16kb_0.8_loops.tsv"
    threads: 4
    resources:
        time= "88:15:00",
        mem= "10G"
    shell:
        "mustache --file {input} "
        "--sparsityThreshold 0.7 "
        "--resolution 16kb "
        "--chromosome chr2L chr2R chr3R chr3L chrX "
        "--pThreshold 0.01 " 
        "--sigmaZero 0.8 "
        "-p {threads} "
        "--outfile {output}"


# #####
# FDR
# #####


rule cutoff_treshold_lowres:
   input:
       rules.mustache_lowres.output
   output:
       "mustache/lowres/filt/0.01/{stage}_16kb_0.8_loops.tsv" 
   threads: 2
   shell:
        "awk -F'\t' -v OFS='\t' 'NR == 1 || $7 <= 0.01' {input} > {output}" 

rule count_dist_lo: 
    input:
        rules.cutoff_treshold_lowres.output 
    output:
        "mustache/lowres/filt2/0.01/{stage}_16kb_0.8_loops_size.tsv"
    threads: 1
    shell:
        "awk -F'\t' -v OFS='\t' 'NR == 1 {{ print $0, \"NinthColumn\"; next }} {{ print $0, $5 - $3 }}' "
        "{input} "
        "> {output}"


# ###################
# scale detect
# ###################

rule cutoff_detect_lowres: ##need to put dist
   input:
       rules.count_dist_lo.output
   output:
       "mustache/lowres/{stage}_16kb_0.8_loops_size5.tsv" 
   threads: 2
   shell:
        "tail -n +2 {input} | "
        "awk -F'\t' -v OFS='\t' '($8 <= 5 && $9 >= 1750000)' | "
        "sortBed -i stdin | "
        "cut -f1-6 -d $'\t' "
        "> {output}" 
    

# ###################
# Loop Calling
# ###################

rule mustache_hires: 
    input:
        rules.cooler_zoomify_subset.output
    output:
        "mustache/hires/{stage}_{res}_1.6_loops.tsv"
    threads: 4
    resources:
        time= "88:15:00",
        mem= "10G"
    shell:
        "mustache --file {input} "
        "--sparsityThreshold 0.7 "
        "--resolution {wildcards.res} "
        "--chromosome chr2L chr2R chr3R chr3L chrX "
        "--pThreshold 0.1 "
        "--sigmaZero 1.6 "
        "-p {threads} "
        "--outfile {output}"

# #####
# FDR
# #####

rule cutoff_treshold_hires:
   input:
       "mustache/hires/{stage}_{res}_1.6_loops.tsv"
   output:
       "mustache/hires/filt/{cutoff}/{stage}_{res}_1.6_loops.tsv" 
   threads: 2
   shell:
        "awk -F'\t' -v OFS='\t' 'NR == 1 || $7 <= {wildcards.cutoff}' {input} > {output}" 

rule count_dist_hi: 
    input:
        rules.cutoff_treshold_hires.output 
    output:
        "mustache/hires/filt2/{cutoff}/{stage}_{res}_1.6_loops_size.tsv"
    threads: 1
    shell:
        "awk -F'\t' -v OFS='\t' 'NR == 1 {{ print $0, \"NinthColumn\"; next }} {{ print $0, $5 - $3 }}' "
        "{input} "
        "> {output}"

# ###################
# scale detect
# ###################

rule cutoff_detect_and_dist_hires_250:
   input:
       "mustache/hires/filt2/{cutoff}/{stage}_250_1.6_loops_size.tsv"
   output:
       "mustache/final/{cutoff}/{stage}/{stage}_250_1.6_loops.tsv" 
   threads: 2
   shell:
        "tail -n +2 {input} | "
        "awk -F'\t' -v OFS='\t' '($9 >= 5000)' | "
        "sortBed -i stdin | "
        "cut -f1-6 -d $'\t' "
        "> {output}" 

rule cutoff_detect_and_dist_hires_1k: 
   input:
       "mustache/hires/filt2/{cutoff}/{stage}_1kb_1.6_loops_size.tsv"
   output:
       "mustache/final/{cutoff}/{stage}/{stage}_1kb_1.6_loops_size_5.tsv" 
   threads: 2
   shell:
        "tail -n +2 {input} | "
        "awk -F'\t' -v OFS='\t' '($8 <= 5 && $9 >= 50000)' | "
        "sortBed -i stdin | "
        "cut -f1-6 -d $'\t' "
        "> {output}" 

rule cutoff_detect_and_dist_hires_2k: 
   input:
       "mustache/hires/filt2/{cutoff}/{stage}_2kb_1.6_loops_size.tsv"
   output:
       "mustache/final/{cutoff}/{stage}/{stage}_2kb_1.6_loops_size_5.tsv" 
   threads: 2
   shell:
        "tail -n +2 {input} | "
        "awk -F'\t' -v OFS='\t' '($8 <= 5 && $9 >= 175000)' | "
        "sortBed -i stdin | "
        "cut -f1-6 -d $'\t' "
        "> {output}" 


rule cutoff_detect_and_dist_hires_4k: 
   input:
       "mustache/hires/filt2/{cutoff}/{stage}_4kb_1.6_loops_size.tsv"
   output:
       "mustache/final/{cutoff}/{stage}/{stage}_4kb_1.6_loops_size_5.tsv" 
   threads: 2
   shell:
        "tail -n +2 {input} | "
        "awk -F'\t' -v OFS='\t' '($8 <= 5 && $9 >= 600000)' | "
        "sortBed -i stdin | "
        "cut -f1-6 -d $'\t' "
        "> {output}" 


# ###################
# combine
# ###################

rule combine_1:
    input:
        a = "mustache/final/{cutoff}/{stage}/{stage}_250_1.6_loops.tsv",
        b =  "mustache/final/{cutoff}/{stage}/{stage}_1kb_1.6_loops_size_5.tsv" 
    output:
        q = temp("mustache/final/{cutoff}/{stage}/{stage}_both.bed"),
        r = temp("mustache/final/{cutoff}/{stage}/{stage}_250.bed"),
        s = temp("mustache/final/{cutoff}/{stage}/{stage}_1kb.bed"),
        final = "mustache/final/{cutoff}/{stage}/{stage}_250_1kb.combo.bed"
    shell:
        "pairToPair -slop 1000 -a {input.a} -b {input.b} > {output.q} && "
        "pairToPair -slop 1000 -a {input.a} -b {input.b} -type notboth > {output.r} && "
        "pairToPair -slop 1000 -a {input.b} -b {input.a} -type notboth > {output.s} && "
        "cat {output.q} {output.r} {output.s} | "
        "cut -f1-6 -d $'\t' | "
        "sortBed -i stdin "
        "> {output.final}"


rule combine_2:
    input:
        a = rules.combine_1.output.final,
        b =  "mustache/final/{cutoff}/{stage}/{stage}_2kb_1.6_loops_size_5.tsv" 
    output:
        q = temp("mustache/final/{cutoff}/{stage}/{stage}_both250_1kb.bed"),
        r = temp("mustache/final/{cutoff}/{stage}/{stage}_250_1kb.bed"),
        s = temp("mustache/final/{cutoff}/{stage}/{stage}_2kb.bed"),
        final = "mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb.combo.bed"
    shell:
        "pairToPair -slop 1000 -a {input.a} -b {input.b} > {output.q} && "
        "pairToPair -slop 1000 -a {input.a} -b {input.b} -type notboth > {output.r} && "
        "pairToPair -slop 1000 -a {input.b} -b {input.a} -type notboth > {output.s} && "
        "cat {output.q} {output.r} {output.s} | "
        "cut -f1-6 -d $'\t' | "
        "sortBed -i stdin "
        "> {output.final}"

rule combine_3:
    input:
        a = rules.combine_2.output.final,
        b =  "mustache/final/{cutoff}/{stage}/{stage}_4kb_1.6_loops_size_5.tsv" 
    output:
        q = temp("mustache/final/{cutoff}/{stage}/{stage}_both250_1kb_2kb.bed"),
        r = temp("mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb.bed"),
        s = temp("mustache/final/{cutoff}/{stage}/{stage}_4kb.bed"),
        final = "mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb_4kb.combo.bed"
    shell:
        "pairToPair -slop 1000 -a {input.a} -b {input.b} > {output.q} && "
        "pairToPair -slop 1000 -a {input.a} -b {input.b} -type notboth > {output.r} && "
        "pairToPair -slop 1000 -a {input.b} -b {input.a} -type notboth > {output.s} && "
        "cat {output.q} {output.r} {output.s} | "
        "cut -f1-6 -d $'\t' | "
        "sortBed -i stdin "
        "> {output.final}"


rule combine_4:
    input:
        a = rules.combine_3.output.final,
        b =  rules.cutoff_detect_lowres.output 
    output:
        q = temp("mustache/final/{cutoff}/{stage}/{stage}_both250_1kb_2kb_4kb.bed"),
        r = temp("mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb_4kb.bed"),
        s = temp("mustache/final/{cutoff}/{stage}/{stage}_16kb.bed"),
        final = "mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb_4kb_16kb.combo.bed"
    shell:
        "pairToPair -slop 1000 -a {input.a} -b {input.b} > {output.q} && "
        "pairToPair -slop 1000 -a {input.a} -b {input.b} -type notboth > {output.r} && "
        "pairToPair -slop 1000 -a {input.b} -b {input.a} -type notboth > {output.s} && "
        "cat {output.q} {output.r} {output.s} | "
        "cut -f1-6 -d $'\t' | "
        "sortBed -i stdin "
        "> {output.final}"

#Remove loops where anchors overlap regions of low-mappability
rule remove_marginal_loops:
    input: 
        loops = "mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb_4kb_16kb.combo.bed",
        marginals = "outputs/marginals/final/dm6_marginals_2kb.bed"
    output:
        "mustache/final/{cutoff}/{stage}/{stage}_250_1kb_2kb_4kb_16kb.clean.bed"
    shell:
        "pairToBed -a {input.loops} -b {input.marginals} -type neither "
        "> {output}" 

# ###################
# clodius
# ###################

#For HiGlass visualization
rule clodius_aggregate_loops: 
    input: 
        rules.remove_marginal_loops.output
    output:
        "mustache/final/clodius/{cutoff}/clean/{stage}_250_1kb_2kb_4kb_16kb.clean.multires" 
    shell:
        "clodius aggregate bedpe "
        "--assembly dm6 "
        "--has-header "
        "--chr1-col 1 --from1-col 2 --to1-col 3 "
        "--chr2-col 4 --from2-col 5 --to2-col 6 "
        "--output-file {output} "
        "{input} " 