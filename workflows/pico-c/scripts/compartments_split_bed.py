import pandas as pd
import sys
import os

# Check for command-line arguments
if len(sys.argv) != 2:
    print("Usage: python split_bed.py <input_file>")
    sys.exit(1)

input_file = sys.argv[1]
input_dir = os.path.dirname(input_file)
input_filename = os.path.basename(input_file)

# Extract bin size from the filename
bin_size = input_filename.split("_")[3]  # Assumes the bin size is the 4th part of the filename

# Load the BED file
df = pd.read_csv(input_file, sep="\t", header=None)

# Determine the unique compartment types (A/B or A.1, A.2, etc.)
df[3] = df[3].apply(lambda x: "".join(x.split(".")[:2]))
unique_compartments = df[3].unique()

# Generate filename without extension
base_filename, ext = os.path.splitext(input_filename)

# Split and save separate BED files
for compartment in unique_compartments:
    subset = df[df[3] == compartment].copy()
    subset["compartment"] = compartment  # Add the compartment name to the output
    output_file = os.path.join(input_dir, f"{base_filename}_{compartment}{ext}")
    subset.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Saved: {output_file}")

print("Splitting complete.")
