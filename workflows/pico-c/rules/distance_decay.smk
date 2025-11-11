#distance_decay.smk

# ###################
# Distance decay
# ###################


rule calc_distance_decay:
    input: 
        ancient("{stage}/hic_ice/{stage}_{res}.hic")
    output:
        "outputs/distance_decay/{stage}/{stage}_{res}_expected_values_all.txt",
        "outputs/distance_decay/{stage}/{stage}_{res}_expected_values_per_chrom.txt"
    resources:
        time= "08:15:00",
        mem= "50G"
    shell: 
        "python scripts/calc_distance_decay.py {input} "
        "outputs/distance_decay/{wildcards.stage}/{wildcards.stage}_{wildcards.res}"