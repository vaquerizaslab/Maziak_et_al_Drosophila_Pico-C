#coolpup_aggregates.smk

#used for paper#
##-------------------------------------------------------------------------##
 ##  make expected files for faster normalization when running coolpup.py ## 
##-------------------------------------------------------------------------##

rule calculate_expected: 
    input:
        "mcool/{stage}_{res}.mcool"
    output:
        "outputs/expected/{stage}_expected_{res}.tsv"
    threads: 4
    resources:
        time= "02:00:00",
        mem= "10G"
    shell: 
        "cooltools expected-cis "
        "--nproc {threads} "
        "--ignore-diags 2 "
        "-o {output} "
        "{input}::resolutions/{wildcards.res}"


##-------------------------------------------------------------------------##
##                  coolpup.py for boundaries                              ## 
##-------------------------------------------------------------------------##


rule pileup_boundaries_all:
    input:
        mcool = "mcool/{stage}_{res}.mcool",
        boundaries = "outputs/boundaries/final/{stage}_0.45_16kb_500bp_boundaries.nocento.bed",
        expected = rules.calculate_expected.output
    output:
        "outputs/coolpup/boundaries/clpy/{stage}_{res}_boundaries_{flank}.all.clpy"
    resources:
        time= "08:00:00",
        mem= "30G"
    threads: 7
    shell:
        "coolpup.py {input.mcool}::resolutions/{wildcards.res} "
        "{input.boundaries} "
        "--features_format bed "
        "--local "
        "--expected {input.expected} "
        "--ignore_diags 2 "
        "--nproc {threads} "
        "--flank {wildcards.flank} "
        "--outname {output}"


rule plot_pileup_boundaries_all:
    input:
        rules.pileup_boundaries_all.output
    output: 
        "outputs/coolpup/boundaries/plots/{stage}_{res}_boundaries_{flank}.all.pdf"
    resources:
        time= "02:00:00",
        mem= "10G"
    shell: 
        "plotpup.py "
        "--font_scale 3.5 "
        "--vmin 0.5 "
        "--vmax 1.5 "
        "--dpi 1000 "
        "--height 3 "
        "--input_pups {input} "
        "--output {output}"


##-------------------------------------------------------------------------##
##                 coolpup.py for Loops                                    ## 
##-------------------------------------------------------------------------##


rule pileup_loops_bydist:
    input:
        mcool = "mcool/{stage}_{res}.mcool",
        loops = "mustache/final/0.01/{stage}/{stage}_250_1kb_2kb_4kb_16kb.clean.bed",
        expected = rules.calculate_expected.output
    output:
        "outputs/coolpup/loops/clpy/{stage}_{res}_loops_{flank}.bydist.clpy"
    resources:
        time= "28:00:00",
        mem= "50G"  
    threads: 7
    shell:
        "coolpup.py {input.mcool}::resolutions/{wildcards.res} "
        "{input.loops} "
        "--features_format bedpe "
        "--expected {input.expected} "
        "--mindist 5000 "
        "--by_distance 0 100000 500000 1000000 2000000 "
        "--ignore_diags 2 "
        "--nproc {threads} "
        "--flank {wildcards.flank} "
        "--outname {output}"


rule plot_pileup_loops_bydist:
    input:
        rules.pileup_loops_bydist.output
    output: 
        pdf = "outputs/coolpup/loops/plots/{stage}_{res}_loops_{flank}.bydist.pdf",
        bed = "outputs/coolpup/loops/plots/{stage}_{res}_loops_{flank}.bydist.bedpe"
    resources:
        time= "02:00:00",
        mem= "10G"
    shell: 
        "plotpup.py "
        "--vmin 0.25 "
        "--vmax 25 "
        "--dpi 600 "
        "--font_scale 2 "
        "--height 2 "
        "--center 1 "
        "--no_score "
        "--input_pups {input} "
        "--output {output.pdf} "
        "--stripe_sort center_pixel "



##-------------------------------------------------------------------------##
##                 coolpup.py for Compartments                             ## 
##-------------------------------------------------------------------------##


rule calculate_expected_comp: 
    input:
        "mcool/{stage}_{res}.mcool"
    output:
        "outputs/coolpup/compartments/expected_cis/{stage}_expected_{res}.tsv"
    threads: 4
    resources:
        time= "01:00:00",
        mem= "10G"
    shell: 
        "cooltools expected-cis "
        "--nproc {threads} "
        "--ignore-diags 0 "
        "-o {output} "
        "{input}::resolutions/{wildcards.res}"

#Remove any additional compartments with regions of low mappability
rule remove_marginals_comp: 
    input:
        comp = "outputs/compartments_CALDER/ICE/NC14_3kb/sub_compartments/all_sub_compartments_3kb_sub_{comp}.bed" , 
        margs = "outputs/marginals/final/dm6_marginals_2kb.bed"
    output: 
        temp("outputs/compartments_CALDER/beds/NC14_sub_compartments_3kb_{comp}.bed")
    shell: 
        "bedtools intersect -v "
        "-a {input.comp} "
        "-b {input.margs} "
        "-f 0.1 " 
        "> {output}" 

#For removing out of bound regions - this has been an issue with calder 
rule remove_outofbound:
    input: 
        comp = rules.remove_marginals_comp.output,
        chrom = "chrom.sizes.bed"
    output:
        "outputs/compartments_CALDER/coolpup_beds/NC14_sub_compartments_3kb_{comp}.bed"
    shell: 
        "bedtools intersect -a {input.comp} -b {input.chrom} "
        "> {output}"

rule pileups_comp:
    """
    Create pile-ups of A compartments 
    """
    input:
        mcool = "mcool/{stage}_{res}.mcool",
        bed = "outputs/compartments_CALDER/coolpup_beds/NC14_sub_compartments_3kb_{comp}.bed",
        expected = "outputs/coolpup/compartments/expected_cis/{stage}_expected_{res}.tsv"
    output:
        clpy = "outputs/coolpup/compartments/clpy/{stage}_comp{comp}_r{res}.clpy"
    resources:
        mem="150G",
        time="00-04:00:00" 
    threads: 3
    shell:
        "coolpup.py "
        "{input.mcool}::resolutions/{wildcards.res} "
        "{input.bed} "
        "--expected {input.expected} "
        "--ignore_diags 0 "
        "--rescale --rescale_flank 1 --rescale_size 99 "
        "--nproc {threads} "
        "--outname {output.clpy}"

rule plot_comp:
    """
    Plot pile-ups 
    """
    input:
        clpy = "outputs/coolpup/compartments/clpy/{stage}_comp{comp}_r{res}.clpy"
    output: 
        pdf = "outputs/coolpup/compartments/plots/comp_{stage}/{stage}_comp{comp}_r{res}.pdf"
    resources:
        mem="3G",
        time="00-01:00:00"
    shell: 
        "plotpup.py "
        "--vmax 1.75 "
        "--vmin 0.75 "
        "--not_symmetric "
        "--height 3 "
        "--dpi 400 "
        "--no_score "
        "--scale log "
        "--input_pups {input.clpy} "
        "--output {output.pdf}"