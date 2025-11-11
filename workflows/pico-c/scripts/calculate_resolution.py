#!/usr/bin/env python
import argparse
import glob
from multiprocessing import Pool
import os
import sys

import fanc
from fanc.tools.general import human_format
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def resolution_parser():
    parser = argparse.ArgumentParser(
        prog='fanc resolution',
        description='Determine the resolution of a Hi-C experiment'
    )

    parser.add_argument(
        'input',
        help='Directory containing binned Hi-C matrices.'
    )

    parser.add_argument(
        'output',
        nargs='?',
        help='Output file for results. If not provided, will write to stdout'
    )

    parser.add_argument(
        '-i', '--interactions-threshold',
        type=int,
        default=1000,
        help='Threshold of interactions per bin. Default: %(default)d'
    )

    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=1,
        help='Number of threads used for parallel processing of files. Default: %(default)d'
    )

    parser.add_argument(
        '-p', '--plot',
        help='Path for saving resolution plot'
    )

    return parser


def _bins_passing_threshold(hic_file: str, threshold: int) -> float:
    """Percentage of bins passing threshold."""
    hic = fanc.load(hic_file, mode='r')
    marginals = hic.marginals(masked=False, norm=False)
    return hic.bin_size, sum(marginals > threshold) / len(marginals)

def resolution(argv, **kwargs):
    parser = resolution_parser()
    args = parser.parse_args(argv[1:])

    input_dir = args.input
    output_file = os.path.expanduser(args.output) if args.output is not None else None
    plot_file = os.path.expanduser(args.plot) if args.plot is not None else None
    threshold = args.interactions_threshold
    threads = args.threads

    input_files = sorted(
        glob.iglob(f"{input_dir}/*.hic"),
        key=lambda file: fanc.load(file, mode='r').bin_size
    )

    with Pool(threads) as p:
        results = p.starmap(
            _bins_passing_threshold,
            ((file, threshold) for file in input_files)
        )

    with open(output_file, 'w') if output_file is not None else sys.stdout as f:
        for result, file in sorted(zip(results, input_files), reverse=True):
            f.write(f"{os.path.basename(file)}\t{result[0]}\t{result[1]}\n")

    if plot_file is not None:
        import matplotlib
        matplotlib.use('agg')

        bin_size, bins_passing = list(zip(*results))
        bins_passing = [x * 100 for x in bins_passing]

        fig, ax = plt.subplots()
        ax.plot(bin_size, bins_passing)
        ax.axhline(80, linestyle='--', color='black')
        ax.set_xscale('log')
        ax.set_xticks(bin_size)
        ax.invert_xaxis()
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: f"{human_format(x)}b")
        )
        ax.set_xlabel('Bin size')
        ax.set_ylabel(f"Bins with > {threshold} contacts %")
        fig.savefig(plot_file)
        plt.close(fig)


if __name__ == '__main__':
    resolution(sys.argv)

