#!/bin/python3

import os
import glob
import numpy as np
import pandas as pd

###
# Script for gathering significant results from rMATS output
# Significant results = IncLevelDifference and FDR above and below their respective, user-set thresholds
###

PSI = float(snakemake.params.psi)
padj = float(snakemake.params.padj)
EVENTS = glob.glob(snakemake.params.rmats_output + "*MATS.JCEC.txt")
output_folder = snakemake.params.script_output_folder

print("Getting significant results from rMATS output.")
for event in EVENTS:
    file_name = os.path.basename(event).replace(".txt", "")
    # Load the file using Pandas
    df = pd.read_csv(event, sep="\t")  # Ensure it's tab-delimited
    
    # Convert the comma-separated values in the four read count columns to lists of floats and compute their averages
    def get_average(col_values):
        values = list(map(float, str(col_values).split(",")))  # Convert comma-separated string to float list
        return np.mean(values) if values else 0  # Compute the mean, handle empty lists safely
    df["avg_IJC_SAMPLE_1"] = df["IJC_SAMPLE_1"].apply(get_average)
    df["avg_SJC_SAMPLE_1"] = df["SJC_SAMPLE_1"].apply(get_average)
    df["avg_IJC_SAMPLE_2"] = df["IJC_SAMPLE_2"].apply(get_average)
    df["avg_SJC_SAMPLE_2"] = df["SJC_SAMPLE_2"].apply(get_average)
    # Compute the overall average across all four columns
    df["overall_avg"] = df[["avg_IJC_SAMPLE_1", "avg_SJC_SAMPLE_1", "avg_IJC_SAMPLE_2", "avg_SJC_SAMPLE_2"]].mean(axis=1)
    # Filtering conditions
    condition_psi_padj = df["IncLevelDifference"].abs().ge(PSI) & df["FDR"].le(padj)
    condition_overall_avg_10 = df["overall_avg"].ge(10)
    # Apply both conditions
    filtered_df = df[condition_psi_padj & condition_overall_avg_10]
    # Drop extra columns before saving
    filtered_df = filtered_df.drop(columns=["avg_IJC_SAMPLE_1", "avg_SJC_SAMPLE_1", "avg_IJC_SAMPLE_2", "avg_SJC_SAMPLE_2", "overall_avg"])
    # Save results
    output_path = os.path.join(snakemake.params.script_output_folder, f"{file_name}.top.txt")
    filtered_df.to_csv(output_path, sep="\t", index=False)
    # Display that file has been dealt with
    print(f"Processed: {file_name}")

print("Getting all significant DAS results' ids in DAS.txt")
# Generator to create set of all significant ids from previously created .top files
list_uniq_DAS = {
    line.split()[1] 
    for filename in os.listdir(output_folder) 
    if os.path.isfile(os.path.join(output_folder, filename))
    for i, line in enumerate(open(os.path.join(output_folder, filename), "r"))
    if i > 0 and len(line.split()) > 1  # Skip header and ensure at least two elements
}
with open(snakemake.output[0], "w") as DAS:
    DAS.write("\n".join(list_uniq_DAS) + "\n")
