# Bulk-RNA-seq-project

#!/usr/bin/env python
# coding: utf-8

import os
import glob
import pandas as pd
import time

# path to folder with featureCounts outputs
path = "/home/tiasha/bulk_rnaseq_work/quants"

files = glob.glob(os.path.join(path, "*_featurecounts.txt"))
files = sorted(files)
print("Files found:", files)

all_counts = []

for file in files:
    start_time = time.time()
    df = pd.read_csv(file, sep="\t", comment="#")
    # the last column contains the counts produced by featureCounts
    sample_col = df.columns[-1]
    sample_name = os.path.basename(file).replace("_featurecounts.txt", "")
    df = df[["Geneid", sample_col]]
    df.rename(columns={sample_col: sample_name}, inplace=True)
    all_counts.append(df)
    elapsed = (time.time() - start_time) / 60
    print(f"Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min")

# merge all on Geneid
counts_matrix = all_counts[0]
for df in all_counts[1:]:
    counts_matrix = counts_matrix.merge(df, on="Geneid", how="outer")

output_file = os.path.join(path, "GSE106305_counts_matrix_3011.csv")
counts_matrix.to_csv(output_file, index=False)

print("\nAll files processed!")
print("Merged matrix shape:", counts_matrix.shape)
print("Saved to:", output_file)
