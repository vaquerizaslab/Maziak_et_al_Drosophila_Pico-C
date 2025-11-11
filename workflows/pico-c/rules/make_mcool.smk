#-------------------#
# Make mcool files
#-------------------#

rule make_cool:
    input:
        "{stage}/hic_ice/{stage}_{res}.hic"
    output:
        temp("cool/{stage}_{res}.cool")
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

rule cooler_zoomify:
    input:
        rules.make_cool.output
    output: 
        "mcool/{stage}_{res}.mcool"
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