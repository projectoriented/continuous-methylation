#!/usr/bin/env python3
"""
Script is tuned for Snakemake and purpose is for merging all percent methylated intermediary files.
"""

import pandas as pd

# Get the DSS columns from one file
if snakemake.wildcards.analysis_type == "two_group":
    target_columns = ["chr", "start", "end", "length", "nCG", "meanMethy1", "meanMethy2", "diff.Methy",
                      "areaStat"]
elif snakemake.wildcards.analysis_type == "model_based":
    target_columns = ["chr", "start", "end", "length", "nCG", "areaStat"]

dss_df = pd.read_table(snakemake.input.all_samples[0], header=0, usecols=target_columns)

# Get the sample columns
samples_df = pd.concat(
    [pd.read_table(x, header=0, usecols=lambda c: c not in target_columns) for x in snakemake.input.all_samples],
    axis=1
)

pd.concat(
    [dss_df, samples_df], axis=1
).to_csv(snakemake.output.chrom_summary_table, sep='\t', header=True, index=False)