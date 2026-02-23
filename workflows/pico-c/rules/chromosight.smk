# chromosight.smk

#This follows https://chromosight.readthedocs.io/en/latest/notebooks/g1_metaphase_yeast_example.html
#subsampled scoring will slightly vary run-to-run for samples being subsampled, but it will consistently give the same trend

# ---------------# 
##  Boundaries  ##
# ---------------# 


rule chromosight_boundaries:
    input: 
        mcool = "mcool/{stage}_125bp.mcool",
        bound = "boundary_clusters/NC14_all_boundaries.bed"
    output: 
        "chromosight/scores/nc14_boundaries/{stage}_{res}_sub187837762.pdf"
    resources:
        time= "20:00:00",
        mem= "30G"
    threads: 10
    shell: 
        "chromosight quantify "
        "--pattern borders "
        "--subsample 187837762 "
        #"--perc-undetected=40 "
        #"--perc-zero 40 "
        "--n-mads=0 "
        "-t {threads} "
        "{input.bound} "
        "{input.mcool}::resolutions/{wildcards.res} "
        "chromosight/scores/nc14_boundaries/{wildcards.stage}_{wildcards.res}_sub187837762"


rule chromosight_boundaries_Jabba:
    input: 
        mcool = "mcool/{stage}_125bp.mcool",
        bound = "boundary_clusters/NC14_all_boundaries.bed"
    output: 
        "chromosight/scores/nc14_boundaries/{stage}_{res}_sub310191338.pdf"
    resources:
        time= "20:00:00",
        mem= "30G"
    threads: 10
    shell: 
        "chromosight quantify "
        "--pattern borders "
        "--subsample 310191338 "
        #"--perc-undetected=40 "
        #"--perc-zero 40 "
        "--n-mads=0 "
        "-t {threads} "
        "{input.bound} "
        "{input.mcool}::resolutions/{wildcards.res} "
        "chromosight/scores/nc14_boundaries/{wildcards.stage}_{wildcards.res}_sub310191338"

# ---------------# 
##  Loops  ##
# ---------------# 

rule chromosight_loops_nc14:
    input: 
        mcool = "mcool/{stage}_125bp.mcool",
        loops = "loop_anchor_clustering/loops/NC14_250_1kb_2kb_4kb_16kb.clean.chromosight.bedpe"
    output: 
        "chromosight/scores/nc14_loops/{stage}_{res}_sub187837762.tsv"
    resources:
        time= "20:00:00",
        mem= "30G"
    threads: 10
    shell: 
        "chromosight quantify "
        "--subsample 187837762 "
        "--perc-undetected=40 "
        "--perc-zero 40 "
        "--n-mads=0 "
        "-t {threads} "
        "{input.loops} "
        "{input.mcool}::resolutions/{wildcards.res} "
        "chromosight/scores/nc14_loops/{wildcards.stage}_{wildcards.res}_sub187837762"


rule chromosight_loops_Jabba:
    input: 
        mcool = "mcool/{stage}_125bp.mcool",
        loops = "loop_anchor_clustering/loops/NC14_250_1kb_2kb_4kb_16kb.clean.chromosight.bedpe"
    output: 
        "chromosight/scores/nc14_loops_dkd/{stage}_{res}_sub310191338.tsv"
    resources:
        time= "20:00:00",
        mem= "30G"
    threads: 10
    shell: 
        "chromosight quantify "
        "--subsample 310191338 "
        "--perc-undetected=40 "
        "--perc-zero 40 "
        "--n-mads=0 "
        "-t {threads} "
        "{input.loops} "
        "{input.mcool}::resolutions/{wildcards.res} "
        "chromosight/scores/nc14_loops_dkd/{wildcards.stage}_{wildcards.res}_sub310191338"


# ------------------# 
##  Combine Loops  ##
# ------------------# 

rule subset_chromo_loops_nc14: 
    input: 
        expand("chromosight/scores/nc14_loops/{{stage}}_{res}_sub187837762.tsv", 
            res = ["1000", "2000", "4000", "8000"])
    output:
        "chromosight/scores/nc14_loops/{stage}_all_loops.tsv"
    shell: 
        "python scripts/merge_chromo_scores.py {input} > {output}"


rule subset_chromo_loops_nc14_Jabba: 
    input: 
        expand("chromosight/scores/nc14_loops_dkd/{{stage}}_{res}_sub310191338.tsv", 
            res = ["1000", "2000", "4000", "8000"])
    output:
        "chromosight/scores/nc14_loops_dkd/{stage}_all_loops.tsv"
    shell: 
        "python scripts/merge_chromo_scores.py {input} "
        "> {output}"

