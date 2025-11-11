# 3D Genome Reorganization Foreshadows Zygotic Genome Activation in *Drosophila*

This repository contains the computational analyses associated with the study  
**“3D Genome Reorganization Foreshadows Zygotic Genome Activation in *Drosophila***”.

The code provided here complements the computational methods described in the manuscript. Please note that while we aim for easy reproducibility, this repository contains the original code that was developed over time and used for the analyses, and may not represent a fully self-contained workflow for reproducing all results in the paper.

Processed and raw FASTQ datasets generated in this study and used in the workflows are available under the accession *[E-MTAB-14477](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14477)* at the ArrayExpress database.

---

## Overview

The repository contains the four main computational workflows described in the paper, along with relevant scripts and supporting files.

1. **Pico-C Data Processing**  
   - Preprocessing, alignment, and quality control of Pico-C datasets.  
   - Generation and normalization of contact maps for downstream analysis.
   - Boundary detection, loop calling, and compartment identification from Pico-C data. Source data and outputs of this analysis are additionally available in the *Supplementary Data* accompanying the manuscript.
   - This will cover a bulk of the paper's data
   > Note: knockdown Pico-C's and α-Amanitin and their respective controls were processed in the same manner. 

2. **RNA-seq Analysis (α-Amanitin and Control Embryos) and ATAC-seq Processing (Control vs. Zelda–GAF Double Knockdowns)**  
   - Initial alignment and processing steps. 
   - Read mapping, peak calling, and normalization.  
   - Comparative chromatin accessibility analysis between genotypes.  

3. **ChromHMM Analysis**  
   - Example Snakemake workflows used to process publicly available datasets.

In addition, relevant R Markdown (`.Rmd`) files used for figure-specific analyses are provided in the `Rmd` folder.

---
