
#-------------------#
# Plot
#-------------------#

rule plot_merge:
    input:
        ancient("{stage}/hic_ice/{stage}_{res}.hic")
    output:
        "qc/region_plots/{stage}_{res}.png"
    resources:
        time= "08:15:00",
        mem= "50G"
    shell:
        "fancplot -o {output} "
        "chr2L:1mb-7mb "
        "-p triangular "
        "-vmin 0.0001 -vmax 0.1 "
        "-l {input}"


#-------------------#
# Find resolution
#-------------------#

rule resolution_merge:
    input:
        ancient(expand("{{stage}}/hic_ice/{{stage}}_{res}.hic", 
            stage = stages, 
            res = test_res))
    output:
        stats = "{stage}/hic_ice/{stage}.resolution.tsv",
        plot = "{stage}/hic_ice/{stage}.resolution.pdf"
    resources:
        time= "08:15:00",
        mem= "50G"
    threads: 2
    shell:
        "python scripts/calculate_resolution.py -t {threads} "
        "-p {output.plot} "
        "{wildcards.stage}/hic_ice/ "
        "{output.stats}"


#-------------------#
# Run FANC PCA
#-------------------#

rule pca: 
    input: 
        ancient(expand("fanc/hic/{id}_{res}.hic", 
            id = IDs.sample_id, 
            res = ["100kb"]))
    output:
        "fanc/pca/fanc_pca_{res}_nodiag_ice_s100k.pca"
    shell: 
        "fanc pca "
        "--min-distance {wildcards.res} --max-distance 1Mb "
        "--no-zeros "
        "--sample-size 100000 "
        "{input} "
        "{output}"
    
