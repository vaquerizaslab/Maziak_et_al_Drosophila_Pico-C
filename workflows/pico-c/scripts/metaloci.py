#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import math
from pathlib import Path

import numpy as np
import pandas as pd
from splot.esda import plot_local_autocorrelation
from splot.esda import plot_moran

import metaloci
from metaloci.parsers import get_parser

import warnings
warnings.filterwarnings('ignore')

# Variables to set for function
plotsize = 10      # Size of plots
logcutoff = 1.     # Min cut-off for log10(Hi-C interactions)
minsig    = 1.     # Min signal data per bin (pseudo-counts)
influence = 1.5    # Influence of consecutive particles
bfact     = 2      # Buffer space for Voronoi
mli_per   = 10_000  # Permutations to assess MLI p-value
signipval = 0.05   # Significance for selecting a MLI p-value


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    # Function arguments as variables
    gene: str = args.gene_name
    region: str = args.coordinates
    cutoff: float = args.cutoff
    smoothing: float = args.smoothing
    dataset: str = args.dataset
    midp: int = args.midpoint
    sigdata: str = "WhateverYouWant"
    hicfile: str = args.hic
    sigfile: str = args.signal
    highlightfile: str = args.regions_of_interest
    outdir: str = args.output_directory.rstrip('/')
    cmap: str = args.colourmap
    vmin = args.min_value
    vmax = args.max_value

    # Make sure the output directory exists
    Path(outdir).mkdir(parents=True, exist_ok=True)

    # Get matrix
    corrected_matrix = metaloci.HiC.load(hicfile, region, log=True)
    print(f"Bin size {hicfile}: {corrected_matrix.resolution:,}bp for region {corrected_matrix.region}")

    metaloci.plot_matrix(
        corrected_matrix,
        title="Corrected Hi-C data",
    )
    plt.savefig(f"{outdir}/{gene}_{dataset}_hic.pdf")

    # Get only the significant pairs
    sig_matrix = corrected_matrix.significant_interactions(cutoff=cutoff)

    metaloci.plot_matrix(
        sig_matrix,
        title="Signficant interactions"
    )
    plt.savefig(f"{outdir}/{gene}_{dataset}_sighic.pdf")

    # Modify the matrix and transorm to restraints
    restraints = sig_matrix.to_Restraints(smoothing=smoothing)
    print(f"Sparse matrix contains {restraints.counts:,} restraints")

    metaloci.plot_matrix(
        restraints,
        title="Restraints"
    )
    plt.savefig(f"{outdir}/{gene}_{dataset}_restraints.pdf")

    # Get the KK layout
    print("Calculating Kamada-Kawai layout...")
    kkgraph = restraints.to_KamadaKawaiGraph()

    # Save the KK layout as JSON
    kkgraph.save(f"{outdir}/{gene}_{dataset}_graph.json")

    # Collect the highlights
    if highlightfile is not None:
        highlights = metaloci.bed_highlights(highlightfile, kkgraph.region, kkgraph.resolution)
    else:
        highlights = None

    # See the Kamada-Kawai layout
    print("Plotting Kamada-Kawai...")
    metaloci.draw_kk_layout(kkgraph, cmap="coolwarm", highlights=highlights)
    plt.savefig(f"{outdir}/{gene}_{dataset}_kk.pdf")

    # Collect signal data
    print(f"Getting signal data for {sigdata}...")
    if sigfile.endswith(".bw") or sigfile.endswith(".bigwig"):
        signal = metaloci.bigwig_signal(sigfile, kkgraph.region, kkgraph.resolution, minsig)
        signal = np.log10(np.nan_to_num(signal, nan=np.nanmin(signal)))
        show_colourbar = True
    elif sigfile.endswith(".bed"):
        signal = metaloci.bed_signal(sigfile, kkgraph.region, kkgraph.resolution)
        cmap = metaloci.bed_colourmap(sigfile)
        show_colourbar = False
    print(f"Signal stats - min: {np.min(signal)}, max: {np.max(signal)}, length: {len(signal)}")

    # See the map with signal colors
    print("Mapping into Kamada-Kawai signal data...")
    _, ax = metaloci.draw_kk_layout(
        kkgraph,
        signal=signal,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        highlights=highlights,
        colourbar=show_colourbar
    )
    if sigfile.endswith(".bed"):
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=idx+1, markerfacecolor=cmap((idx + 0.5) / cmap.N), markersize=8)
            for idx
            in range(cmap.N)
            if idx + 1 in signal
        ]
        ax.legend(handles=legend_elements, loc="center right", prop={"size": 8}, frameon=False)

    plt.savefig(f"{outdir}/{gene}_{dataset}_kk_signal.pdf")
    plt.close()

    # From points to Voronoi polygons
    print("Voronoing...")
    vgraph = kkgraph.to_VoronoiGraph()

    _, ax = metaloci.draw_gaudi_plot(
        vgraph,
        signal=signal,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        highlights=highlights,
        options={"legend": show_colourbar, "legend_kwds": {"shrink": 0.8}}
    )
    if sigfile.endswith(".bed"):
        # voronoi_data.plot(ax=axs[1], column="v", cmap=cmap, edgecolor="lightgrey", scheme="User_Defined", classification_kwds={"bins": np.arange(1.5, 21.5)})
        ax.legend(handles=legend_elements, loc="center right", prop={"size": 8}, frameon=False)

    plt.savefig(f"{outdir}/{gene}_{dataset}_gaudi.pdf")
    plt.close()

    if sigfile.endswith(".bed"):
        exit()

    # Calculate Global Moran's I
    gmoran = vgraph.global_moran(signal)
    if (math.isnan(gmoran.I)):
        print(">>> ERROR... MATH NAN... signal data contains nans")
        exit()
    print("Global Moran Index: {:4.2f} with p-val: {:6.4f} ".format(gmoran.I, gmoran.p_sim))

    # Plot Global Moran's I
    plot_moran(gmoran, zstandard=False, figsize=(10, 4))
    plt.savefig(f"{outdir}/{gene}_{dataset}_gmoran.pdf".format(gene, dataset))
    plt.close()

    # Calculate Local Moran's I
    lmoran = vgraph.local_moran(signal, permutations=mli_per)
    sig_points = len(lmoran.p_sim[(lmoran.p_sim < signipval) & (lmoran.q == 1)])
    print(f"There are a total of {sig_points} significant points in Local Moran Index")

    # Plot Local Moran's I
    print("Final LMI plot...")
    plot_local_autocorrelation(lmoran, vgraph.shapes, "signal", p=signipval, cmap="coolwarm", scheme="EqualInterval")
    plt.savefig(f"{outdir}/{gene}_{dataset}_lmoran.pdf")
    plt.close()

    # Create a BED file that contains some Moran's information for each bin
    beddata = metaloci.create_table(
        vgraph,
        local_i=lmoran.Is,
        pval=lmoran.p_sim,
        quadrant=lmoran.q,
        global_i=gmoran.I
    )
    beddata.to_csv(f"{outdir}/{gene}_{dataset}_moran.bed", sep="\t", index=False)

    # Print significant metaloci
    mask = (lmoran.p_sim <= signipval) & (lmoran.q <= 2)
    hhpoints = vgraph.shapes[mask].idx

    if highlights is not None:
        for highlight in highlights:
            nmetaloci = sum(1 for bin in highlight.bins if bin in hhpoints.values)
            print(
                f">>> Highlighted region for {highlight.name} has {nmetaloci} significant metaloci out of "
                f"{len(highlight.bins)} bins it covers."
            )
