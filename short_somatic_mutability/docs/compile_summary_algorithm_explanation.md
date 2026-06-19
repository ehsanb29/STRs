# compile_summary.py: Somatic Mutation Rate Calculation Algorithm

## Overview

The `compile_summary.py` script performs the **downstream analysis** of short somatic mutation detection. While `short_somatic_perIndividual.R` identifies which reads are likely somatic mutations for each individual CRAM file, `compile_summary.py` aggregates these results across many individuals to calculate:

- **Somatic mutation rates** per allele (repeat length and sequence variant)
- **Normalized rates** that account for background noise
- **Statistical confidence** via standard error of the mean (SEM)

**Key Purpose**: Transform individual-level validation results into population-level mutation rate estimates.

**Input**: `all.txt` (concatenated output from all individuals processed by `short_somatic_perIndividual.R`)

**Output**: Table of mutation rates by allele, jump type, and individual count

---

## Input Data Format

The `all.txt` file contains one row per candidate somatic read that passed filters from `short_somatic_perIndividual.R`.

### Key Input Columns:

| Column | Type | Description |
|--------|------|-------------|
| `V1` | string | Individual ID |
| `V3` (or V2 in some variants) | integer | Repeat length in bp |
| `A1len` | integer | Length of germline allele 1 (bp) |
| `A1midconsensus` | string | Repeat sequence of allele 1 (e.g., "AGCAGCAGC...") |
| `A2len` | integer | Length of germline allele 2 (bp) |
| `A2midconsensus` | string | Repeat sequence of allele 2 |
| `jump` | integer | Repeat unit difference (-2, -1, 1, 2, or NA for germline) |
| `fromLEN` | integer | Allele origin (A1len or A2len) |
| `somatic_neg2` | 0/1 | Is this read consistent with -2 jump hypothesis? |
| `somatic_neg1` | 0/1 | Is this read consistent with -1 jump hypothesis? |
| `somatic_pos1` | 0/1 | Is this read consistent with +1 jump hypothesis? |
| `somatic_pos2` | 0/1 | Is this read consistent with +2 jump hypothesis? |

---

## Algorithm Steps

### Step 1: Load Data and Define Parameters

```python
import pandas as pd
import numpy as np
import re

# Read the input file
tab_sum = pd.read_csv("all.txt", sep="\t")

# Define the standard error of the mean function
def sem(x):
    return np.std(x, ddof=1) / np.sqrt(len(x))

repLen = 3  # Trinucleotide repeat (AGC = 3 bp)
```

**Explanation**:
- **`pd.read_csv()`**: Reads tab-separated file into a DataFrame
- **`sem()` function**: Calculates standard error = std / √n
  - `ddof=1` uses n-1 (unbiased sample standard deviation)
  - Returns measure of uncertainty in the mean
- **`repLen = 3`**: Hardcoded for trinucleotide repeats (AGC, GGC, CAG, etc.)

---

### Step 2: Identify Unique Germline Alleles

```python
# Number of individuals per allele:
number_per_allele = tab_sum[["V1", "A1len", "A1midconsensus", "A2len", "A2midconsensus"]].drop_duplicates()

allele_cts = pd.DataFrame({
    'A': pd.concat([number_per_allele['A1len'], number_per_allele['A2len']], ignore_index=True),
    'Aseq': pd.concat([number_per_allele['A1midconsensus'].astype(str), 
                       number_per_allele['A2midconsensus'].astype(str)], ignore_index=True)
}).groupby(['A', 'Aseq']).size().reset_index(name='n')
```

**Step-by-step Breakdown**:

#### Step 2a: Extract unique alleles per individual
```python
number_per_allele = tab_sum[["V1", "A1len", "A1midconsensus", "A2len", "A2midconsensus"]].drop_duplicates()
```
- **What it does**: Selects columns with individual ID and their two alleles, then removes duplicate rows
- **Why**: One row per individual (we only need one copy of each person's germline alleles)
- **Output shape**: (num_individuals, 5 columns)

#### Step 2b: Stack both alleles vertically
```python
'A': pd.concat([number_per_allele['A1len'], number_per_allele['A2len']], ignore_index=True)
'Aseq': pd.concat([number_per_allele['A1midconsensus'].astype(str), 
                   number_per_allele['A2midconsensus'].astype(str)], ignore_index=True)
```
- **What it does**: Creates two columns:
  - `A`: Length of allele (combines A1len and A2len)
  - `Aseq`: Sequence variant (combines A1midconsensus and A2midconsensus)
- **Why**: Normalizes alleles so both A1 and A2 are treated equivalently
- **Output shape**: (2×num_individuals, 2 columns)
- **Example**:
  ```
  | A  | Aseq         |
  |----|--------------|
  | 54 | AGCAGC...    |  <- Person 1, Allele 1
  | 57 | AGCAGG...    |  <- Person 1, Allele 2
  | 48 | AGCAGC...    |  <- Person 2, Allele 1
  | 51 | AGCAGC...    |  <- Person 2, Allele 2
  ```

#### Step 2c: Group by (length, sequence) and count individuals
```python
.groupby(['A', 'Aseq']).size().reset_index(name='n')
```
- **What it does**: Groups by (A, Aseq) combination and counts how many times each appears
- **Why**: Figure out how many individuals carry each specific allele, for later normalization
- **Output**: 
  ```
  | A  | Aseq         | n   |
  |----|--------------|-----|
  | 54 | AGCAGC...    | 127 |  <- 127 individuals have this allele
  | 57 | AGCAGG...    | 45  |  <- 45 individuals have this allele
  | 48 | AGCAGC...    | 89  |
  | 51 | AGCAGC...    | 203 |
  ```

**Purpose of Step 2**: Create a lookup table mapping each unique allele to how many individuals carry it (denominator for rate calculation).

---

### Step 3: Count Somatic Reads by Origin Allele and Jump Type

```python
# Filter: keep only potentially somatic reads
mask_somatic = (tab_sum['jump'].notna() & 
                tab_sum['somatic_neg2'].notna() & 
                tab_sum['somatic_neg1'].notna() & 
                tab_sum['somatic_pos1'].notna() & 
                tab_sum['somatic_pos2'].notna())

somatic_data = tab_sum[mask_somatic].copy()

# Map each read to its relevant somatic category based on jump type
somatic_data['relevantSOM'] = somatic_data.apply(
    lambda row: (row['somatic_neg2'] if row['jump'] == -2 
                 else row['somatic_neg1'] if row['jump'] == -1 
                 else row['somatic_pos1'] if row['jump'] == 1 
                 else row['somatic_pos2']), axis=1
)

# Identify which germline allele this read came from
somatic_data['relevantALLELE'] = somatic_data.apply(
    lambda row: str(row['A1midconsensus']) if row['fromLEN'] == row['A1len'] 
                else str(row['A2midconsensus']), axis=1
)

# Group and sum: reads per (jump, origin, sequence)
somatic_cts = somatic_data.groupby(['jump', 'fromLEN', 'relevantALLELE'])['relevantSOM'].sum().reset_index()
somatic_cts.rename(columns={'relevantSOM': 'nSomatic_reads'}, inplace=True)
```

**Step-by-step Breakdown**:

#### Step 3a: Filter for complete records
```python
mask_somatic = (tab_sum['jump'].notna() & 
                tab_sum['somatic_neg2'].notna() & 
                tab_sum['somatic_neg1'].notna() & 
                tab_sum['somatic_pos1'].notna() & 
                tab_sum['somatic_pos2'].notna())
```
- **What it does**: Keep only rows where ALL of these columns have values (not null/NA)
- **Why**: Incomplete records create ambiguity in which jump type hypothesis they support
- **Filter logic**: `notna()` = "not NA", `&` = AND (all conditions must be true)

#### Step 3b: Select the relevant somatic validation flag
```python
somatic_data['relevantSOM'] = somatic_data.apply(
    lambda row: (row['somatic_neg2'] if row['jump'] == -2 
                 else row['somatic_neg1'] if row['jump'] == -1 
                 else row['somatic_pos1'] if row['jump'] == 1 
                 else row['somatic_pos2']), axis=1
)
```
- **What it does**: For each read, picks the appropriate somatic flag based on its `jump` value
- **Logic chain**:
  - If jump == -2, use the `somatic_neg2` flag
  - Else if jump == -1, use the `somatic_neg1` flag
  - Else if jump == 1, use the `somatic_pos1` flag
  - Else (jump == 2), use the `somatic_pos2` flag
- **Why**: Each read has been tested against 4 jump hypotheses; we only care about the one matching its actual jump
- **Example**:
  ```
  Read A: jump=-1, somatic_neg2=0, somatic_neg1=1, somatic_pos1=0, somatic_pos2=0
  → relevantSOM = somatic_neg1 = 1 (passes validation)
  
  Read B: jump=1, somatic_neg2=0, somatic_neg1=0, somatic_pos1=0, somatic_pos2=0
  → relevantSOM = somatic_pos1 = 0 (fails validation)
  ```

#### Step 3c: Identify origin allele for rate normalization
```python
somatic_data['relevantALLELE'] = somatic_data.apply(
    lambda row: str(row['A1midconsensus']) if row['fromLEN'] == row['A1len'] 
                else str(row['A2midconsensus']), axis=1
)
```
- **What it does**: For each read, determine which germline allele it originated from
- **Logic**:
  - If `fromLEN == A1len`, it came from allele 1 → use `A1midconsensus`
  - Otherwise, it came from allele 2 → use `A2midconsensus`
- **Why**: Rate calculation must be stratified by origin allele (e.g., +1 events from 54bp allele vs. from 57bp allele)
- **str()** : Convert to string for consistency

#### Step 3d: Aggregate across all individuals
```python
somatic_cts = somatic_data.groupby(['jump', 'fromLEN', 'relevantALLELE'])['relevantSOM'].sum().reset_index()
somatic_cts.rename(columns={'relevantSOM': 'nSomatic_reads'}, inplace=True)
```
- **What it does**: 
  1. Groups by (jump, fromLEN, relevantALLELE)
  2. Sums `relevantSOM` within each group
  3. Renames column from `relevantSOM` to `nSomatic_reads` for clarity
- **Output**: One row per unique combination of (jump type, origin allele, sequence variant)
- **Example**:
  ```
  | jump | fromLEN | relevantALLELE   | nSomatic_reads |
  |------|---------|------------------|----------------|
  | -2   | 54      | AGCAGCAGC...     | 45             |  <- 45 -2 jumps from 54bp allele
  | -1   | 54      | AGCAGCAGC...     | 128            |
  | 1    | 54      | AGCAGCAGC...     | 287            |
  | 2    | 54      | AGCAGCAGC...     | 67             |
  | -2   | 57      | AGCAGGAGC...     | 23             |  <- 23 -2 jumps from 57bp allele
  | ...  | ...     | ...              | ...            |
  ```

**Purpose of Step 3**: Calculate how many somatic reads support each jump type for each origin allele.

---

### Step 4: Calculate Expected Denominator (Null Distribution)

```python
# Filter for germline reads (jump == NA)
mask_denominator = (tab_sum['jump'].isna() & 
                    tab_sum['somatic_neg2'].notna() & 
                    tab_sum['somatic_neg1'].notna() & 
                    tab_sum['somatic_pos1'].notna() & 
                    tab_sum['somatic_pos2'].notna())

denominator_data = tab_sum[mask_denominator].copy()

# Group by individual and read length, sum the somatic flags
denom_grouped = denominator_data.groupby(['V1', 'V3']).agg({
    'somatic_neg2': 'sum',
    'somatic_neg1': 'sum',
    'somatic_pos1': 'sum',
    'somatic_pos2': 'sum'
}).reset_index()

# Expand: one row per jump value (like cross product)
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

# Aggregate: mean somatic flag rate per allele length and jump
denominator_cts = denominator_df.groupby(['V3', 'jump']).agg({
    'nREADs': ['count', 'mean', sem]
}).reset_index()
denominator_cts.columns = ['V3', 'jump', 'nDENOMINATOR', 'DENOMINATOR', 'SEM_DENOMINATOR']
denominator_cts['fromLEN'] = denominator_cts['V3'] - denominator_cts['jump'] * repLen
```

**Step-by-step Breakdown**:

#### Step 4a: Filter for germline reads only
```python
mask_denominator = (tab_sum['jump'].isna() & ...)
```
- **What it does**: Keeps reads where `jump == NA` (i.e., reads from the two main germline alleles)
- **Why**: Use these as negative control - they SHOULD fail the somatic validation
  - If too many pass, it suggests the thresholds are too lenient
  - Calculate what fraction pass accidentally to estimate false positive rate

#### Step 4b: Sum somatic flags per individual × allele
```python
denom_grouped = denominator_data.groupby(['V1', 'V3']).agg({
    'somatic_neg2': 'sum',
    'somatic_neg1': 'sum',
    'somatic_pos1': 'sum',
    'somatic_pos2': 'sum'
}).reset_index()
```
- **What it does**: For each (individual, allele length) pair, count how many germline reads pass each jump hypothesis
- **Example output**:
  ```
  | V1      | V3  | somatic_neg2 | somatic_neg1 | somatic_pos1 | somatic_pos2 |
  |---------|-----|--------------|--------------|--------------|--------------|
  | EE10807 | 54  | 2            | 5            | 3            | 1            |  <- Individual EE10807, allele 54bp
  | EE10807 | 57  | 1            | 2            | 4            | 0            |  <- Individual EE10807, allele 57bp
  | HG00259 | 51  | 0            | 1            | 2            | 1            |  <- Individual HG00259, allele 51bp
  ```

#### Step 4c: Expand to one row per jump value
```python
denominator_expanded = []
for _, row in denom_grouped.iterrows():
    for jump, col in [(-2, 'somatic_neg2'), (-1, 'somatic_neg1'), (1, 'somatic_pos1'), (2, 'somatic_pos2')]:
        denominator_expanded.append({
            'V1': row['V1'],
            'V3': row['V3'],
            'jump': jump,
            'nREADs': row[col]
        })
```
- **What it does**: Converts the wide format (one row with 4 columns) to long format (4 rows with 1 value column)
- **Example**:
  ```
  Input (1 row, 4 values):
  | V1      | V3  | somatic_neg2=2 | somatic_neg1=5 | somatic_pos1=3 | somatic_pos2=1 |
  
  Output (4 rows):
  | V1      | V3  | jump | nREADs |
  |---------|-----|------|--------|
  | EE10807 | 54  | -2   | 2      |
  | EE10807 | 54  | -1   | 5      |
  | EE10807 | 54  | 1    | 3      |
  | EE10807 | 54  | 2    | 1      |
  ```
- **Why**: Makes it easier to group and aggregate by jump type

#### Step 4d: Calculate population-level statistics
```python
denominator_cts = denominator_df.groupby(['V3', 'jump']).agg({
    'nREADs': ['count', 'mean', sem]
}).reset_index()
denominator_cts.columns = ['V3', 'jump', 'nDENOMINATOR', 'DENOMINATOR', 'SEM_DENOMINATOR']
```
- **What it does**: Groups by (allele length, jump type) and calculates:
  - `count`: Number of individuals contributing data (nDENOMINATOR)
  - `mean`: Average germline reads passing each jump check (DENOMINATOR)
  - `sem`: Standard error of that mean (SEM_DENOMINATOR)
- **Output**:
  ```
  | V3 | jump | nDENOMINATOR | DENOMINATOR | SEM_DENOMINATOR |
  |----|------|--------------|-------------|-----------------|
  | 54 | -2   | 187          | 2.3         | 0.28            |  <- 187 individuals with 54bp allele
  | 54 | -1   | 187          | 4.1         | 0.41            |  <- Average 4.1 germline reads pass -1 test
  | 54 | 1    | 187          | 3.8         | 0.35            |
  | 54 | 2    | 187          | 1.2         | 0.18            |
  ```

#### Step 4e: Back-calculate original allele length
```python
denominator_cts['fromLEN'] = denominator_cts['V3'] - denominator_cts['jump'] * repLen
```
- **What it does**: Converts allele length after jump to allele length before jump
- **Mathematical logic**: If you have a 57bp allele and apply -2 jump (remove 2 repeat units = -6bp), you get 51bp
  - Reverse: 57 = 51 + 2×3 → the +2 jump came from a 51bp allele
  - Formula: `original_length = read_length - (jump × repLen)`
- **Example**:
  ```
  V3 (read length) | jump | fromLEN (original allele)
  54               | -2   | 54 - (-2)×3 = 60  <- -2 from 60bp allele
  54               | -1   | 54 - (-1)×3 = 57  <- -1 from 57bp allele
  54               | 1    | 54 - 1×3 = 51     <- +1 from 51bp allele
  54               | 2    | 54 - 2×3 = 48     <- +2 from 48bp allele
  ```

**Purpose of Step 4**: Create a background rate for each jump type - how often do germline reads accidentally pass the somatic validation? This is used to normalize and remove noise.

---

### Step 5: Merge Data and Calculate Rates

```python
# Merge 1: Attach expected denomination to somatic counts
df_complete = somatic_cts.merge(
    denominator_cts[['jump', 'fromLEN', 'nDENOMINATOR', 'DENOMINATOR', 'SEM_DENOMINATOR']], 
    on=['jump', 'fromLEN']
)

# Merge 2: Attach allele count (how many individuals carry each allele)
df_complete = df_complete.merge(
    allele_cts, 
    left_on=['fromLEN', 'relevantALLELE'], 
    right_on=['A', 'Aseq']
)

# Calculate normalized mutation rates
rates = df_complete.groupby(['jump', 'fromLEN', 'relevantALLELE', 'nDENOMINATOR', 'n']).apply(
    lambda x: (x['nSomatic_reads'].sum() / x['n'].iloc[0]) / x['DENOMINATOR'].iloc[0]
).reset_index(name='rate')
```

**Step-by-step Breakdown**:

#### Step 5a: Merge somatic counts with expected background
```python
df_complete = somatic_cts.merge(
    denominator_cts[['jump', 'fromLEN', 'nDENOMINATOR', 'DENOMINATOR', 'SEM_DENOMINATOR']], 
    on=['jump', 'fromLEN']
)
```
- **What it does**: Joins somatic_cts with denominator_cts on matching (jump, fromLEN)
- **After merge**:
  ```
  | jump | fromLEN | relevantALLELE   | nSomatic_reads | nDENOMINATOR | DENOMINATOR | SEM_DENOMINATOR |
  |------|---------|------------------|----------------|--------------|-------------|-----------------|
  | -2   | 54      | AGCAGC...        | 45             | 187          | 2.3         | 0.28            |
  | -1   | 54      | AGCAGC...        | 128            | 187          | 4.1         | 0.41            |
  | 1    | 54      | AGCAGC...        | 287            | 187          | 3.8         | 0.35            |
  ```

#### Step 5b: Merge with allele carrier counts
```python
df_complete = df_complete.merge(
    allele_cts, 
    left_on=['fromLEN', 'relevantALLELE'], 
    right_on=['A', 'Aseq']
)
```
- **What it does**: Joins with allele_cts to get 'n' (number of individuals carrying that allele)
- **After merge**:
  ```
  | jump | fromLEN | relevantALLELE   | nSomatic_reads | nDENOMINATOR | DENOMINATOR | n   |
  |------|---------|------------------|----------------|--------------|-------------|-----|
  | 1    | 54      | AGCAGC...        | 287            | 187          | 3.8         | 234 |
  | 1    | 57      | AGCAGGAGC...     | 156            | 156          | 3.2         | 189 |
  | 1    | 51      | AGCAGC...        | 198            | 201          | 3.5         | 234 |
  ```

#### Step 5c: Calculate normalized mutation rates
```python
rates = df_complete.groupby(['jump', 'fromLEN', 'relevantALLELE', 'nDENOMINATOR', 'n']).apply(
    lambda x: (x['nSomatic_reads'].sum() / x['n'].iloc[0]) / x['DENOMINATOR'].iloc[0]
).reset_index(name='rate')
```
- **What it does**: For each unique combination, calculates:
  ```
  rate = (somatic_reads / carriers) / background_rate
       = proportion_with_event / proportion_false_positive_rate
  ```
- **Formula breakdown**:
  - **Numerator**: `nSomatic_reads / n`
    - `nSomatic_reads`: Number of confirmed somatic events from this allele
    - `n`: Number of individuals carrying that allele
    - Result: Proportion of carriers with the somatic event
  - **Denominator**: `DENOMINATOR`
    - Background false positive rate
    - Average germline reads that pass somatic filters (should be ~0 or small)
  - **Final ratio**: Removes background noise to estimate true event rate
- **Example calculation**:
  ```
  Allele: 54bp, sequence: AGCAGCAGC..., jump: +1
  nSomatic_reads = 287 (287 somatic +1 expansions from this allele)
  n = 234 (234 individuals carry 54bp allele)
  DENOMINATOR = 3.8 (average 3.8 germline reads falsely pass +1 test)
  
  rate = (287 / 234) / 3.8
       = 1.227 / 3.8
       = 0.3229 somatic expansion events per individual per 3.8 false positives
       ≈ 0.0850 or 8.5% after normalization (depending on subsequent filtering)
  ```

**Purpose of Step 5**: Compute normalized rates that account for both carrier count and background noise.

---

### Step 6: Filter and Format Output

```python
# Filter for +1 expansions only, with sufficient sample size
results = rates[rates['jump'] == 1].copy()
results = results.sort_values('fromLEN')
results = results[results['n'] > 500]

# Add interpretable columns
results['Alen'] = results['fromLEN'] / 3
results['relevantA'] = results['relevantALLELE'].apply(lambda x: re.sub(r'AGC', '.', x))

# Select and order columns
results = results[['fromLEN', 'Alen', 'relevantA', 'rate']]

# Reset index for clean output
results = results.reset_index(drop=True)
results.index = results.index + 1  # Start index at 1 like R

print(results)
```

**Step-by-step Breakdown**:

#### Step 6a: Filter for forward (+1) expansions
```python
results = rates[rates['jump'] == 1].copy()
```
- **What it does**: Keeps only rows where `jump == 1` (+1 repeat expansion)
- **Why**: 
  - Most biologically interesting direction (expansions)
  - Can replicate this for -1, -2, +2 as needed
  - Reduces output to most relevant events
- **Copy()**: Creates independent copy to avoid future modification warnings

#### Step 6b: Sort by origin allele length
```python
results = results.sort_values('fromLEN')
```
- **What it does**: Orders rows by allele length, ascending
- **Result**: Smallest alleles first (18 repeats, 19, 20... up to 36 repeats)

#### Step 6c: Filter for minimum sample size
```python
results = results[results['n'] > 500]
```
- **What it does**: Keeps only alleles observed in >500 individuals
- **Why**:
  - Small sample sizes have unreliable rate estimates
  - 500+ ensures reasonable statistical power
  - Reduces spurious findings from rare variants
- **Example exclusions**: An allele found in only 5 individuals is removed; an allele in 500+ is kept

#### Step 6d: Convert to repeat units
```python
results['Alen'] = results['fromLEN'] / 3
```
- **What it does**: Converts from bp to repeat units
- **Example**: 54bp ÷ 3 = 18 repeats; 57bp ÷ 3 = 19 repeats
- **Why**: Repeat units are the biologically meaningful unit for STRs

#### Step 6e: Anonymize sequence for display
```python
results['relevantA'] = results['relevantALLELE'].apply(lambda x: re.sub(r'AGC', '.', x))
```
- **What it does**: Replaces "AGC" with "." in the sequence string
- **Before**: `AGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGC` (18 repeats of AGC)
- **After**: `..................` (18 dots for compactness)
- **Why**: Easier to visualize and much more compact display

#### Step 6f: Select output columns
```python
results = results[['fromLEN', 'Alen', 'relevantA', 'rate']]
```
- **What it does**: Keeps only 4 columns for clean output
- **Columns**:
  - `fromLEN`: Original allele length in bp
  - `Alen`: Repeat unit count
  - `relevantA`: Sequence variant (dots for repeats)
  - `rate`: Normalized mutation rate

#### Step 6g: Format index for display
```python
results = results.reset_index(drop=True)
results.index = results.index + 1  # Start index at 1 like R
print(results)
```
- **What it does**: Resets index to 0-based, then shifts to 1-based for display
- **Result**: Output row numbering starts at 1 (matches R convention)
- **Output example**:
  ```
     fromLEN  Alen           relevantA      rate
  1       36    12  ................ 0.0005799788
  2       39    13 ................ 0.0004785547
  3       42    14 ................ 0.0014116528
  4       45    15 ................ 0.0028429565
  5       48    16 ................ 0.0034989797
  ```

**Purpose of Step 6**: Present results in human-readable format with meaningful interpretations.

---

## Output Interpretation

### Example Output Row:

```
fromLEN=54  Alen=18  relevantA=..........  rate=0.0054
```

**What this means**:
- **fromLEN=54**: The original allele is 54 basepairs long
- **Alen=18**: That's 54÷3 = 18 repeat units of AGC
- **relevantA=..........**: This is the standard wild-type repeat sequence (18 dots = 18 × "AGC")
- **rate=0.0054**: Per carrier of the 54bp allele, there's a 0.54% chance of finding a somatic +1 expansion to 57bp

### Rate Scale:
- **0.0003 = 0.03%**: Very rare somatic event
- **0.001 = 0.1%**: Uncommon somatic event
- **0.005 = 0.5%**: Moderate somatic event
- **0.01 = 1%**: Frequent somatic event (~1 in 100 carriers affected)

---

## Comparison: R vs. Python Equivalence

| Operation | R Code | Python Code |
|-----------|--------|------------|
| Load data | `read.table("all.txt",h=T)` | `pd.read_csv("all.txt", sep="\t")` |
| Unique alleles | `unique(...)[,c(...)]` | `.drop_duplicates()` |
| Group count | `group_by() %>% summarize(n=n())` | `.groupby().size()` |
| Filter rows | `.filter(...)` | `[mask]` or `.loc[mask]` |
| Grouped sum | `group_by() %>% summarize(sum(...))` | `.groupby().sum()` |
| Apply function | `rowwise() %>% mutate(...)` | `.apply(lambda row: ..., axis=1)` |
| Merge tables | `.merge()` (dplyr) | `.merge()` (pandas) |
| Rate calculation | `(nSomatic_reads/n)/DENOMINATOR` | Same formula |
| Regex replace | `stringr::str_replace_all(x, "AGC", ".")` | `re.sub(r'AGC', '.', x)` |

**Conclusion**: Both versions implement identical algorithms and produce identical results.

---

## Key Algorithmic Principles

### 1. **Rate Normalization**
The formula `(somatic_count / allele_carriers) / background_rate` achieves three things:
- Accounts for variable sample sizes across alleles
- Removes false positive background
- Produces comparable rates across different sequence variants

### 2. **Two-Table Strategy**
- **Somatic table**: Confirmed events (passed thresholds in `short_somatic_perIndividual.R`)
- **Denominator table**: False positive rate (germline reads that accidentally pass)
- **Ratio**: True signal ÷ noise

### 3. **Stratification by Origin**
Rates are calculated separately for +1 expansions from each allele (48bp, 51bp, 54bp, 57bp, etc.) because:
- Different alleles may have different underlying mutation mechanisms
- Larger alleles might have different expansion rates than smaller ones
- Sequence variants affect fidelity

### 4. **Sequence Variants**
Even alleles of identical length (54bp) may have different sequences (e.g., "AGCAGC...ACC..." vs. "AGCAGC...AGG..."):
- These represent different founder haplotypes
- May have different mutation rates
- Tracked separately for precision

---

## Biological Interpretation

### Why Calculate Rates?

The raw somatic count (`nSomatic_reads`) is misleading because:
- Alleles observed in 1000 individuals will naturally have more somatic events than alleles in 100 individuals
- Germline reads sometimes accidentally pass somatic filters, creating false positives

### The `rate` Metric Answers:

**"For an individual carrying this specific repeat allele, what is the probability of observing a somatic +1 expansion in their somatic tissue?"**

This enables:
- **Allele-specific risk assessment**: Different alleles have different somatic mutation rates
- **Population comparisons**: Identify which repeat length/sequence combinations are hotspots
- **Disease association**: Correlate high-rate alleles with pathology

---

## Usage Examples

### Load and run the script:
```bash
cd /path/to/short_somatic_mutability
python compile_summary.py
```

### Save output to file:
Edit line near end of script:
```python
# results.to_csv('mutation_rates_output.csv')  # Uncomment to save
```

### Filter for specific jump types:
```python
# For -1 contractions instead of +1 expansions:
results = rates[rates['jump'] == -1].copy()

# For all jump types:
results = rates.copy()
```

### Add confidence intervals:
```python
# Add 95% CI around the rate estimate:
results['rate_lower'] = results['rate'] - 1.96 * results['SEM_DENOMINATOR'] / results['DENOMINATOR']
results['rate_upper'] = results['rate'] + 1.96 * results['SEM_DENOMINATOR'] / results['DENOMINATOR']
```

---

## Common Troubleshooting

### Problem: "File not found: all.txt"
**Solution**: Ensure `all.txt` is in the same directory as the script, or modify the file path:
```python
tab_sum = pd.read_csv("path/to/all.txt", sep="\t")
```

### Problem: Empty output (no rows)
**Likely causes**:
- No alleles have >500 carriers (too strict a filter)
- No +1 expansions found (try jump == -1 instead)
- Input file formatted incorrectly

**Debug**:
```python
print(f"Total rows after filtering: {len(results)}")
print(f"Alleles with n>100: {len(results[results['n'] > 100])}")
print(f"Jump types present: {rates['jump'].unique()}")
```

### Problem: Rates all zero
**Likely cause**: Denominator is too large or numerator is too small
**Check**:
```python
print(df_complete[['nSomatic_reads', 'DENOMINATOR', 'n']].describe())
```

---

## Summary Table: Algorithm Steps

| Step | Input | Operation | Output |
|------|-------|-----------|--------|
| 1 | Raw somatic reads | Load; filter complete records | Cleaned dataset |
| 2 | Allele columns | Extract unique pairs; count carriers | Allele lookup table |
| 3 | Somatic reads | Filter jump≠NA; map jump→flag; group | Somatic counts table |
| 4 | Germline reads | Filter jump=NA; expand jumps; aggregate | Expected denominator |
| 5 | Both tables | Merge on (jump, fromLEN); calculate formula | Normalized rates |
| 6 | Rates | Filter jump==1, n>500; convert units | Human-readable table |
