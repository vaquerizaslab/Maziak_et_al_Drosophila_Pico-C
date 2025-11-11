#compartments.smk

#------------------------------#
# Generate inputs for `CALDER` #
#------------------------------#

# Get sparse Hi-C matrices with `fanc dump` for each individual chromosome
rule fanc_dump:
    input:
        hic = "{stage}/hic_ice/{stage}_{bin}.hic"
    output:
        matrix = "outputs/compartments_CALDER/sparse_matrix/{stage}_{bin}/{stage}_{bin}_{chrom}_sparseMatrixReg.txt"
    resources:
        mem="10G",
        time="00-01:00:00" # format: "DD-HH:MM:SS"
    threads: 1
    shell:
        "fanc dump -tmp "
        "--subset {wildcards.chrom}--{wildcards.chrom} "
        "{input} "
        "{output.matrix}"

# Slim sparse matrices to pos_x pos_y contact_value

rule slim_sparse_matrix:
    """
    Slim sparse matrix
    """
    input:
        matrix = rules.fanc_dump.output
    output:
        slim_matrix = "outputs/compartments_CALDER/sparse_matrix/{stage}/{stage}_{bin}_{chrom}_sparseMatrixReg_slim.txt"
    threads: 1
    shell:
        "cat {input} | cut -f2,5,7 > {output}"



#----------------------------------#
# Run `CALDER` for each chromosome #
#----------------------------------#

## Run R script from terminal through Snakemake

rule run_CALDER:
    """
    Run CALDER for each chromosome
    """
    input:
        matrix = expand("outputs/compartments_CALDER/sparse_matrix/{stage}/{stage}_{{bin}}_{chrom}_sparseMatrixReg_slim.txt", 
            stage = stages, 
            chrom= ["chr2L", "chr2R", "chr3L", "chr3R", "chrX"]),
        feature_track = "marks_bw/MBT_H3K36me3_merge.bw"
    output:
        comp = "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments.bed"
    params:
        infolder = "outputs/compartments_CALDER/sparse_matrix/{stage}",
        outfolder = "outputs/compartments_CALDER/ICE/{stage}_{bin}"
    resources:
        mem="20G",
        time="00-06:00:00" # format: "DD-HH:MM:SS"
    threads: 2
    shell:
        """
        module load singularityce/3.11.3
        singularity exec --bind /path/ /path/bioconductor_docker_RELEASE_3_17.sif Rscript --vanilla scripts/sc_compCALDER_v2.R -f {params.infolder} -t {input.feature_track} -b {wildcards.bin} -o {params.outfolder} 
        """


rule comaprtment_sort:
    input: 
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments.bed"
    output:
        temp("outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_sorted.bed")
    shell: 
        "python3 scripts/sort_bed_compartments.py {input} {output}"

     
rule comaprtment_nocento: 
    input:
        compartments = "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_sorted.bed",
        cento = "genome/centromeres_dm6_Zenketal2021.bed"
    output:
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_nocento.bed"
    shell: 
        "bedtools subtract "
        "-a {input.compartments} "
        "-b {input.cento} "
        "> {output}"

#--------------------#
# Parse compartments #
#--------------------#

rule parse_sub: 
    input: 
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_nocento.bed"
    output: 
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_sub.bed"
    shell: 
        """
        awk 'BEGIN {{OFS="\\t"}} 
            {{
                split($4, name_parts, ".");
                subcomp = name_parts[1] "." name_parts[2];
                
                if (NR == 1) {{
                    chrom = $1;
                    start = $2;
                    end = $3;
                    compartment = subcomp;
                }} else {{
                    if ($1 == chrom && $2 == (end + 1) && subcomp == compartment) {{
                        end = $3;  # Merge with previous region
                    }} else {{
                        print chrom, start, end, compartment;
                        chrom = $1;
                        start = $2;
                        end = $3;
                        compartment = subcomp;
                    }}
                }}
            }} 
            END {{ print chrom, start, end, compartment; }}' {input} > {output}
        """

rule parse_AB:
    input: 
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_nocento.bed"
    output:
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_ab.bed"
    shell: 
        """
        awk 'BEGIN {{OFS="\\t"}} 
            {{
                subcomp = substr($4, 1, 1);  # Extract only the first character (A or B)
                
                if (NR == 1) {{
                    chrom = $1;
                    start = $2;
                    end = $3;
                    compartment = subcomp;
                }} else {{
                    if ($1 == chrom && $2 == (end + 1) && subcomp == compartment) {{
                        end = $3;  # Merge with previous region
                    }} else {{
                        print chrom, start, end, compartment;
                        chrom = $1;
                        start = $2;
                        end = $3;
                        compartment = subcomp;
                    }}
                }}
            }} 
            END {{ print chrom, start, end, compartment; }}' {input} > {output}
        """


rule split_sub: 
    input:
        rules.parse_sub.output
    output:
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_sub_{comp}.bed"
    shell:
        "python3 scripts/compartments_split_bed.py {input}"


rule split_AB: 
    input:
        rules.parse_AB.output
    output:
        "outputs/compartments_CALDER/ICE/{stage}_{bin}/sub_compartments/all_sub_compartments_{bin}_ab_{comp}.bed"
    shell: 
        "python3 scripts/compartments_split_bed.py {input}"