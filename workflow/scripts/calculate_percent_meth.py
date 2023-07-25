#!/usr/bin/env python3
"""
Script is tuned for Snakemake and purpose is for calculating percent methylated within a region.
"""
import pandas as pd

df = pd.read_table(snakemake.input.sample_bed, header=None, names=["chrom", "start", "n_valid", "n_mod"])
dss_df = pd.read_table(snakemake.input.dss_out, header=0)


# # Calculate the percentage methylated
def get_percnt_methylated(row):
    n_valid_sum = df.loc[df.start.between(row.start, row.end), "n_valid"].sum()
    n_mod_sum = df.loc[df.start.between(row.start, row.end), "n_mod"].sum()

    percnt_methylated = n_mod_sum / n_valid_sum

    return float('{:.2f}'.format(percnt_methylated))


try:
    sample_name = snakemake.wildcards.other_sample
except AttributeError:
    sample_name = snakemake.wildcards.sample

try:
    column_name = f'{sample_name}_{snakemake.wildcards.hap}'
except AttributeError:
    column_name = f'{sample_name}'


dss_df[column_name] = dss_df.apply(lambda row: get_percnt_methylated(row=row), axis=1)
dss_df.fillna(0.0, inplace=True)

dss_df.to_csv(snakemake.output.sample_summary_table, sep='\t', header=True, index=False)
