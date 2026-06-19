# This script compiles the summary of somatic mutability rates for short STRs based on the output from short_somatic_perIndividual.R. It calculates the number of somatic reads, the expected number of reads, and the rate of somatic mutations for each allele length and repeat sequence. The results are filtered to include only expansions (jump == 1) with more than 500 individuals and are saved to a CSV file.
# Usage: python compile_summary.py <input_file_list>
#   <input_file_list>: a text file with one summary file path per line (default: file_list.txt)

import pandas as pd
import numpy as np
import re
import sys

# Read the input file list
input_list = sys.argv[1] if len(sys.argv) > 1 else "file_list.txt"
with open(input_list) as f:
    summary_files = [line.strip() for line in f if line.strip()]

tab_sum = pd.concat(
    [pd.read_csv(fp, sep=r"\s+", engine="python", encoding="utf-8-sig") for fp in summary_files],
    ignore_index=True
)
tab_sum.columns = tab_sum.columns.str.strip()

# Define the standard error of the mean function
def sem(x):
    return np.std(x, ddof=1) / np.sqrt(len(x))

repLen = 3

# Number of individuals per allele:
print(tab_sum.columns)
number_per_allele = tab_sum[["V1", "A1len", "A1midconsensus", "A2len", "A2midconsensus"]].drop_duplicates()
allele_cts = pd.DataFrame({
    'A': pd.concat([number_per_allele['A1len'], number_per_allele['A2len']], ignore_index=True),
    'Aseq': pd.concat([number_per_allele['A1midconsensus'].astype(str), 
                       number_per_allele['A2midconsensus'].astype(str)], ignore_index=True)
}).groupby(['A', 'Aseq']).size().reset_index(name='n')

# Somatic read counts: get the number of reads that appear likely a somatic expansion (+1,+2) 
# or contraction (-1,-2) from each allele length / repeat sequence
mask_somatic = (tab_sum['jump'].notna() & 
                tab_sum['somatic_neg2'].notna() & 
                tab_sum['somatic_neg1'].notna() & 
                tab_sum['somatic_pos1'].notna() & 
                tab_sum['somatic_pos2'].notna())

somatic_data = tab_sum[mask_somatic].copy()
somatic_data['relevantSOM'] = somatic_data.apply(
    lambda row: (row['somatic_neg2'] if row['jump'] == -2 
                 else row['somatic_neg1'] if row['jump'] == -1 
                 else row['somatic_pos1'] if row['jump'] == 1 
                 else row['somatic_pos2']), axis=1
)
somatic_data['relevantALLELE'] = somatic_data.apply(
    lambda row: str(row['A1midconsensus']) if row['fromLEN'] == row['A1len'] 
                else str(row['A2midconsensus']), axis=1
)

somatic_cts = somatic_data.groupby(['jump', 'fromLEN', 'relevantALLELE'])['relevantSOM'].sum().reset_index()
somatic_cts.rename(columns={'relevantSOM': 'nSomatic_reads'}, inplace=True)

# Number of reads expected: number of trustworthy reads (from reads originating from main alleles) 
# that passed likely-somatic filters
mask_denominator = (tab_sum['jump'].isna() & 
                    tab_sum['somatic_neg2'].notna() & 
                    tab_sum['somatic_neg1'].notna() & 
                    tab_sum['somatic_pos1'].notna() & 
                    tab_sum['somatic_pos2'].notna())

denominator_data = tab_sum[mask_denominator].copy()
denom_grouped = denominator_data.groupby(['V1', 'V3']).agg({
    'somatic_neg2': 'sum',
    'somatic_neg1': 'sum',
    'somatic_pos1': 'sum',
    'somatic_pos2': 'sum'
}).reset_index()

# Expand to have one row per jump value
denominator_expanded = []
for _, row in denom_grouped.iterrows():
    for jump, col in [(-2, 'somatic_neg2'), (-1, 'somatic_neg1'), (1, 'somatic_pos1'), (2, 'somatic_pos2')]:
        denominator_expanded.append({
            'V1': row['V1'],
            'V3': row['V3'],
            'jump': jump,
            'nREADs': row[col]
        })

denominator_df = pd.DataFrame(denominator_expanded)

# Group by V3 and jump
denominator_cts = denominator_df.groupby(['V3', 'jump']).agg({
    'nREADs': ['count', 'mean', sem]
}).reset_index()
denominator_cts.columns = ['V3', 'jump', 'nDENOMINATOR', 'DENOMINATOR', 'SEM_DENOMINATOR']
denominator_cts['fromLEN'] = denominator_cts['V3'] - denominator_cts['jump'] * repLen

# Merge dataframes
df_complete = somatic_cts.merge(
    denominator_cts[['jump', 'fromLEN', 'nDENOMINATOR', 'DENOMINATOR', 'SEM_DENOMINATOR']], 
    on=['jump', 'fromLEN']
)
df_complete = df_complete.merge(
    allele_cts, 
    left_on=['fromLEN', 'relevantALLELE'], 
    right_on=['A', 'Aseq']
)

# Calculate rates
rates = df_complete.groupby(['jump', 'fromLEN', 'relevantALLELE', 'nDENOMINATOR', 'n']).apply(
    lambda x: (x['nSomatic_reads'].sum() / x['n'].iloc[0]) / x['DENOMINATOR'].iloc[0]
).reset_index()
rates.rename(columns={rates.columns[-1]: 'rate'}, inplace=True)

# Filter and format results
results = rates[rates['jump'] == 1].copy()
results = results.sort_values('fromLEN')
results = results[results['n'] > 500]
results['Alen'] = results['fromLEN'] / 3
results['relevantA'] = results['relevantALLELE'].apply(lambda x: re.sub(r'AGC', '.', x))
results = results[['fromLEN', 'Alen', 'relevantA', 'rate']]

# Reset index for clean output
results = results.reset_index(drop=True)
results.index = results.index + 1  # Start index at 1 like R

print(results)

# Optionally save to CSV
results.to_csv('output.csv')
