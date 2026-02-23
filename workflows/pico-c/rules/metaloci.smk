#metaloci.smk 


##----------------------##
##      Bin genome      ##
##----------------------##

rule makewindows_genome: 
    input: 
        "genome/dm6.genome"
    output:
        temp("genome/dm6.genome.bins") 
    shell: 
        "bedtools makewindows "
        "-g {input} "
        "-w 200000 "
        "-s 100000 "
        "> {output}" 

rule remove_lowmapping: 
    input:
        bins = rules.makewindows_genome.output,
        margs = "genome/cento_marg_dm6.bed"
    output:
        "genome/dm6.genome.bed" 
    shell: 
        "bedtools intersect -v "
        "-a {input.bins} "
        "-b {input.margs} "
        "-f 0.05 "
        "> {output}" 

rule make_regionsfile: 
    input: 
        "genome/dm6.genome.bed"
    output: 
        "genome/regionsfile.txt"
    shell: 
        "awk '{{print $1\":\"$2\"-\"$3}}' "
        "{input} > {output}"



##----------------------##
##       Metaloci       ##
##----------------------##


## Pol II ## 
rule metaloci_polII_nc12:
    input: 
        hic = "/home/nmaziak/mnt/storage/vaquerizas/nmaziak/20240219_microc_analysis/NC12/hic_ice/NC12_500bp.hic",
        signal = "../relevant_files/metaloci_bw/NC12_polII_merge.bw"
    output:
        directory("metaloci_scores/polII/NC12/NC12_polII_{region}")
    shell:
        "python3 scripts/metaloci.py "
        "-o {output}  -c 1.2 "
        "{input.hic} "
        "{input.signal} "
        "{wildcards.region}"


rule metaloci_polII_nc14:
    input: 
        hic = "/home/nmaziak/mnt/storage/vaquerizas/nmaziak/20240219_microc_analysis/NC14/hic_ice/NC14_500bp.hic",
        signal = "../relevant_files/metaloci_bw/MBT_nc14em_PRJNA266120_polIIser5_merge.bw"
    output:
        directory("metaloci_scores/polII/NC14/NC14_polII_{region}")
    shell:
        "python3 scripts/metaloci.py "
        "-o {output}  -c 1.2 "
        "{input.hic} "
        "{input.signal} "
        "{wildcards.region}"


## other marks ##
rule metaloci_marks_nc12:
    input: 
        hic = "/home/nmaziak/mnt/storage/vaquerizas/nmaziak/20240219_microc_analysis/NC12/hic_ice/NC12_500bp.hic",
        signal = "../relevant_files/metaloci_bw/preMBT_{factor}_merge.bw"
    output:
        directory("metaloci_scores/{factor}/NC12/NC12_{factor}_{region}")
    shell:
        "python3 scripts/metaloci.py "
        "-o {output}  -c 1.2 "
        "{input.hic} "
        "{input.signal} "
        "{wildcards.region}"

rule metaloci_marks_nc14:
    input: 
        hic = "/home/nmaziak/mnt/storage/vaquerizas/nmaziak/20240219_microc_analysis/NC14/hic_ice/NC14_500bp.hic",
        signal = "../relevant_files/metaloci_bw/MBT_{factor}_merge.bw"
    output:
        directory("metaloci_scores/{factor}/NC14/NC14_{factor}_{region}")
    shell:
        "python3 scripts/metaloci.py "
        "-o {output}  -c 1.2 "
        "{input.hic} "
        "{input.signal} "
        "{wildcards.region}"



