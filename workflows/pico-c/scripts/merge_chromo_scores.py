#!/usr/bin/env python3

import sys
import pandas as pd
import os

def extract_resolution(filename):
    """Extract resolution as the second underscore-separated token from the basename."""
    base = os.path.basename(filename)
    parts = base.split("_")
    try:
        return int(parts[1])
    except (IndexError, ValueError):
        return None

def main():
    if len(sys.argv) < 2:
        print("Usage: python merge_chromo_scores.py file1.tsv file2.tsv ... > output.tsv", file=sys.stderr)
        sys.exit(1)

    dfs = []
    for file in sys.argv[1:]:
        res = extract_resolution(file)
        if res is None:
            print(f"Could not extract resolution from filename: {file}", file=sys.stderr)
            continue

        try:
            df = pd.read_csv(file, sep="\t")
        except Exception as e:
            print(f"Failed to read {file}: {e}", file=sys.stderr)
            continue

        df["resolution"] = res
        # Create a unique loop ID based on genomic coordinates (not bin1/bin2)
        df["loop_id"] = (
            df["chrom1"] + ":" + df["start1"].astype(str) + "-" + df["end1"].astype(str) + "_" +
            df["chrom2"] + ":" + df["start2"].astype(str) + "-" + df["end2"].astype(str)
        )
        dfs.append(df)

    if not dfs:
        print("No valid input files processed.", file=sys.stderr)
        sys.exit(1)

    all_data = pd.concat(dfs, ignore_index=True)

    # Convert score to numeric, invalid entries become NaN
    all_data["score"] = pd.to_numeric(all_data["score"], errors="coerce")

    # Drop rows with missing scores (NaN)
    all_data = all_data.dropna(subset=["score"])

    # Sort by loop_id and resolution (lowest res = highest priority)
    all_data = all_data.sort_values(by=["loop_id", "resolution"])

    # Keep the highest-resolution entry with a valid score for each loop
    best_loops = all_data.drop_duplicates(subset="loop_id", keep="first")

    # Drop helper column
    best_loops = best_loops.drop(columns=["loop_id"])

    # Output to stdout
    best_loops.to_csv(sys.stdout, sep="\t", index=False)

if __name__ == "__main__":
    main()
