# Short Somatic Mutability Algorithm - Detailed Explanation

## Overview
The `short_somatic_perIndividual.R` script detects somatic DNA repeat expansions and contractions in individual samples. It identifies reads that differ from germline alleles, validates them using quality metrics, and flags potentially true somatic mutations.

---

## Workflow: From CRAM File to somatic.sh to R Analysis

### The somatic.sh Orchestration Script

The `somatic.sh` bash script serves as an **orchestrator** that:
1. Extracts repeat region reads from a sequencing CRAM file
2. Filters and pre-processes those reads
3. Prepares the input file (`IID_{ID}.txt`)
4. Calls `short_somatic_perIndividual.R` to perform analysis

**Usage**:
```bash
bash somatic.sh <CRAM_FILE> <LOCUS_SPECIFICATION>
```

**Example**:
```bash
bash somatic.sh EE10807Y.cram "TCF4_chr18_55586154_55586228_AGGAGGAGC_AGCATGAAA_AGC_L6,9_R,"
```

### somatic.sh Execution Flow

#### Step 1: Parse Command-Line Arguments
```bash
FILE="$1"; 
ID=$( basename "$FILE" | sed 's/\.final\.cram$//' | sed 's/\.cram$//' ); 
DAT=$2;
```

**What it does**:
- **FILE**: Input CRAM file path (e.g., "EE10807Y.cram")
- **ID**: Extracts individual identifier from filename using `basename` and `sed`
  - Removes `.final.cram` or `.cram` extension
  - Example: "EE10807Y.cram" → "EE10807Y"
- **DAT**: The locus specification string (passed as argument 2)

#### Step 2: Parse Locus Specification
```bash
CHR=$(echo $DAT | cut -d'_' -f 2);           # Field 2: chromosome
START=$(echo $DAT | cut -d'_' -f 3);         # Field 3: start position
END=$(echo $DAT | cut -d'_' -f 4);           # Field 4: end position
STARTSEQ=$(echo $DAT | cut -d'_' -f 5);      # Field 5: left flank sequence
ENDSEQ=$(echo $DAT | cut -d'_' -f 6);        # Field 6: right flank sequence
```

**Example with "TCF4_chr18_55586154_55586228_AGGAGGAGC_AGCATGAAA_AGC_L6,9_R,"**:
- `CHR="chr18"`
- `START="55586154"`
- `END="55586228"`
- `STARTSEQ="AGGAGGAGC"` (9bp: 6bp flank + 3bp repeat motif)
- `ENDSEQ="AGCATGAAA"` (9bp: 3bp repeat motif + 6bp flank)

#### Step 3: Extract Reads from Target Region Using Samtools

```bash
samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -F 0x4 -F 0x100 -F 0x200 -F 0x400 -F 0x800 \
  "$FILE" ${CHR}:${START}-${END}
```

**Samtools flags explained**:
- `-T ref.fa`: Reference genome file (converts SAM cigar to actual sequences)
- `-F 0x4`: Exclude unmapped reads (flag 4)
- `-F 0x100`: Exclude secondary alignments (flag 100)
- `-F 0x200`: Exclude vendor QC failures (flag 200)
- `-F 0x400`: Exclude duplicate reads (flag 400)
- `-F 0x800`: Exclude supplementary alignments (flag 800)
- `${CHR}:${START}-${END}`: Retrieve only reads in the target region

**Output**: SAM format stream with aligned reads from chr18:55586154-55586228

SAM format columns (used in awk):
- `$2`: FLAG (bit information about read)
- `$5`: MAPQ (mapping quality)
- `$7`: RNEXT (mate reference name, "=" means same as read)
- `$10`: SEQ (sequence)
- `$11`: QUAL (quality string)

#### Step 4: Filter and Extract Repeat Motif Using AWK

```bash
awk -v s=$STARTSEQ -v e=$ENDSEQ -v iid=$ID '
BEGIN {STARTSEQ_TOADD = substr(s,7,3); ENDSEQ_TOADD=substr(e,1,3)}
$5 >= 30 && $7=="=" {
  start=split($10,a,s);end=split($10,b,e); 
  if (start == 2 && end == 2) {
    split(a[2],segMID,e);
    print iid,(and($2,16)==0?"forward":"reverse"),length(STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD), STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD,$10,$11
  }
}' > IID_${ID}.txt
```

**Step-by-step breakdown**:

##### Step 4a: Initialize repeat motif segments
```bash
BEGIN {STARTSEQ_TOADD = substr(s,7,3); ENDSEQ_TOADD=substr(e,1,3)}
```
- **STARTSEQ_TOADD**: Extract characters 7-9 from STARTSEQ (the repeat motif part)
  - Example: "AGGAGG**AGC**" (position 7-9) = "AGC"
- **ENDSEQ_TOADD**: Extract characters 1-3 from ENDSEQ (the repeat motif part)
  - Example: "**AGC**ATGAAA" (position 1-3) = "AGC"
- **Why**: These are the bases that define the repeat boundaries

##### Step 4b: Quality and pairing filter
```bash
$5 >= 30 && $7=="=" {
```
- **$5 >= 30**: MAPQ (mapping quality) must be ≥30 (high confidence mapping)
- **$7 == "="**: Read is properly paired with mate on same chromosome

##### Step 4c: Detect flanking sequences and extract middle segment
```bash
start=split($10,a,s);end=split($10,b,e);
```
- **split($10,a,s)**: Split read sequence (column 10) by STARTSEQ delimiter
  - Splits "...AGGAGGAGC[repeat motif here]AGCATGAAA..." by "AGGAGGAGC"
  - Result: `start=2` if found (2 parts: before + after), or `start=1` if not found
  - Array `a[1]` = part before flank, `a[2]` = part after flank
- **split($10,b,e)**: Split by ENDSEQ
  - Similarly detects end flank presence

```bash
if (start == 2 && end == 2) {
```
- **Validation**: Both flanks must be found exactly once (start==2 AND end==2)
- **Rejects**: Reads missing flanks (likely misaligned or extracted incorrectly)

##### Step 4d: Extract the middle repeat segment
```bash
split(a[2],segMID,e);
```
- **What it does**: Takes the part after STARTSEQ (`a[2]`) and splits by ENDSEQ
- **Result**: `segMID[1]` = the middle between the two flanks (the repeat region)
- **Example**:
  ```
  Original read: "...AGGAGGAGC[AGC AGC AGC]AGCATGAAA..."
  After first split by STARTSEQ: a[2] = "[AGC AGC AGC]AGCATGAAA..."
  After second split by ENDSEQ: segMID[1] = "[AGC AGC AGC]"
  ```

##### Step 4e: Construct output and determine strand
```bash
print iid,(and($2,16)==0?"forward":"reverse"),
      length(STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD), 
      STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD,
      $10,$11
```

**Output fields** (written to `IID_${ID}.txt`):

| Field | Source | Description |
|-------|--------|-------------|
| Col 1 | `iid` | Individual ID |
| Col 2 | `and($2,16)==0?"forward":"reverse"` | Strand (bitwise AND with flag 16) |
| Col 3 | `length(...)` | Repeat length in bp (motif + flanking repeats on each side) |
| Col 4 | `STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD` | Full repeat sequence (3bp + middle + 3bp) |
| Col 5 | `$10` | Original SAM read sequence |
| Col 6 | `$11` | Original SAM quality string |

**Strand detection**: `and($2,16)==0`
- Bitwise AND of FLAG with 16 (reverse strand flag)
- Result 0 = read on forward strand
- Result ≠0 = read on reverse complement

#### Step 5: Count Extracted Reads
```bash
nREADS=$(cat IID_${ID}.txt | wc -l)
echo "Number of reads in region: $nREADS"
```
- Counts how many valid reads passed all filters
- Prints diagnostic message

#### Step 6: Conditional Execution of R Analysis
```bash
if [ $nREADS -gt 0 ]; then
  Rscript short_somatic_perIndividual.R $ID $DAT
fi
```
- **Why conditional?**: No point running R if there are zero reads
- **Passes to R**: Individual ID and locus specification
- **R receives**: Pre-processed reads in `IID_${ID}.txt`

### Example Data Flow

**Original CRAM file excerpt** (SAM format after samtools view):
```
EE10807Y  0  chr18  55586154  30  ...  =  ...  AGGAGGAGCAGCAGCAGCAGCATGAAA  ########...
EE10807Y 16  chr18  55586155  30  ...  =  ...  AGGAGGAGCAGCAGCAGCAGCATGAAA  ########...
...
```

**After somatic.sh filtering** (`IID_EE10807Y.txt`):
```
EE10807Y	forward	18	AGCAGCAGCAGCAGCAGC	AGGAGGAGCAGCAGCAGCAGCATGAAA	########...
EE10807Y	reverse	18	AGCAGCAGCAGCAGCAGC	AGGAGGAGCAGCAGCAGCAGCATGAAA	########...
...
```

**This becomes input to R** as:
- V1 (ID): "EE10807Y"
- V2 (Strand): "forward", "reverse"
- V3 (Length): 18 (in bp, will be 6 repeats of 3bp AGC motif)
- V4 (Repeat seq): "AGCAGCAGCAGCAGCAGC"
- V5 (Full seq): "AGGAGGAGCAGCAGCAGCAGCATGAAA"
- V6 (Quality): ASCII quality string

---

## Input Data

### Command Line Arguments
```
Rscript short_somatic_perIndividual.R ID DAT
```

- **ID**: Individual identifier (e.g., "EE10807Y")
- **DAT**: Locus description string in format:
  ```
  LOCUS_CHR_STARTBP_ENDBP_STARTstr_ENDstr_REP_L[positions]_R[positions]
  ```
  Example: `TCF4_chr18_55586154_55586228_AGGAGGAGC_AGCATGAAA_AGC_L6,9_R,`

### Parsed Parameters
From the DAT string, the script extracts:
- **STARTstr**: Left flank sequence (6bp repeat flanking sequence + 3bp repeat motif)
- **ENDstr**: Right flank sequence (3bp repeat motif + 6bp repeat flanking sequence)
- **repLen**: Length of repeat motif in bp (e.g., 3 for trinucleotide repeats like AGC)
- **L_FLANK**: Positions of "key" bases in left flank (bases where errors could affect repeat count)
- **R_FLANK**: Positions of "key" bases in right flank

### Input File
**File format**: `IID_{ID}.txt` (tab-separated, no header)

| Column | Field | Description |
|--------|-------|-------------|
| V1 | Individual ID | Sample identifier |
| V2 | Strand | "forward" or "reverse" |
| V3 | Repeat Length | Length of repeat sequence in bp |
| V4 | Repeat Sequence | Actual repeat sequence extracted from read |
| V5 | Read Sequence | Full extracted read sequence |
| V6 | Quality String | PHRED quality scores (ASCII encoded) |

---

## Algorithm Steps

### Step 1: Build Consensus Sequences
**Function**: `consensus(SEQUENCES, QUALITIES, START)`

This function is called within the grouping operation to create a high-confidence representative sequence from multiple sequencing reads of the same allele.

#### Code Expression:
```r
SEQ = consensus(SEQUENCES= as.character(V5),QUALITIES= as.character(V6),START= STARTstr)
```

#### Detailed Step-by-Step Processing:

**1a. Locate the Start Position in Each Read**
```r
START_LOCS = as.data.frame(stringr::str_locate(string= SEQUENCES, START))
```
- **What it does**: Uses regex `str_locate()` to find where STARTstr (the left flanking sequence) appears in each read's sequence (V5)
- **Why?**: The start location varies between reads due to read length variation. Finding it allows proper alignment
- **Output**: START_LOCS contains 2 columns: start position and end position of the match

**1b. Split Sequences and Qualities into Character Vectors**
```r
SEQUENCES = stringr::str_split(SEQUENCES,"")  # Convert "AGCAGCAGC" → ["A","G","C","A","G","C","A","G","C"]
QUALITIES = stringr::str_split(QUALITIES,"")  # Convert quality string into individual characters
```
- **What it does**: Breaks sequences into individual bases and quality scores for position-by-position comparison
- **Why?**: Allows comparison of bases at the same genomic position across multiple reads

**1c. Create a Position Grid for Alignment**
```r
SEQdf = data.frame(SPOT=seq(-156,156,by=1), SEQ="", QUAL="", stringsAsFactors = F)
```
- **What it does**: Creates a reference grid with positions from -156 to +156 relative to the repeat start
- **Why?**: Provides a fixed coordinate system to align reads despite their different lengths. Position 0 marks the START location
- **Range explanation**: -156 to +156 encompasses the flanks (6bp + 3bp on each side) plus the repeat region itself

**1d. Stack reads at the Reference Position**
```r
for(i in 1:length(SEQUENCES)){
  ADD = data.frame(SPOT=seq(1,length(SEQUENCES[[i]])) - START_LOCS[i,2] - 2, 
                   SEQ=SEQUENCES[[i]], QUAL=QUALITIES[[i]])
  SEQdf[SEQdf$SPOT %in% ADD$SPOT,"SEQ"] = paste0(SEQdf[SEQdf$SPOT %in% ADD$SPOT,"SEQ"], ADD$SEQ)
  SEQdf[SEQdf$SPOT %in% ADD$SPOT,"QUAL"] = paste0(SEQdf[SEQdf$SPOT %in% ADD$SPOT,"QUAL"], ADD$QUAL)
}
```
- **What it does**: 
  - For each read i, calculates its positions relative to the START location
  - `seq(1, length()) - START_LOCS[i,2] - 2` centers the repeat at position 0
  - Concatenates sequences and qualities from all reads at each position
- **Why?**: Creates a multiple sequence alignment where all reads are anchored at the same start position
- **Output**: SEQdf now contains stacked sequences and quality scores at each position

**1e. Call Consensus Base using Quality-Filtered Majority Vote**
```r
for(i in 1:nrow(SEQdf)){
  opts = stringr::str_split(SEQdf$SEQ[i],"")[[1]][utf8ToInt(SEQdf$SEQ[i])-33 >= 25]
  if(length(opts) > 0){
    sortedVAL = sort(table(opts), decreasing=TRUE)
    if(length(sortedVAL)==1){
      CONSENSUS = paste0(CONSENSUS, names(sortedVAL[1]), collapse = "")
    }else{
      # Handle ties or multiple options
      CONSENSUS = paste0(CONSENSUS, names(sortedVAL[1]), collapse = "")
    }
  }else{
    CONSENSUS = paste0(CONSENSUS, "X", collapse="")  # Unknown base if no high-quality bases
  }
}
```
- **What it does**:
  - `utf8ToInt(SEQdf$SEQ[i]) - 33` converts ASCII quality characters to PHRED scores
  - Filters to only bases where PHRED ≥ 25 (≈ 99.7% accuracy)
  - Counts frequency of each base (A, C, G, T)
  - Selects the most frequent base as consensus
  - If no high-quality bases exist, inserts "X" (unknown)
- **Why filter on Q≥25?**: Ensures we only use high-confidence bases. Lower quality bases are essentially noise
- **Output**: A single consensus sequence string (e.g., "AGCAGCAGC...")

**Summary of Step 1**:
- **Input**: Multiple reads with sequences and quality scores
- **Process**: Align, stack, and vote on bases using quality-weighted majority rule
- **Output**: One high-confidence consensus sequence per allele group

**Why Consensus Building is Critical**:
Multiple sequencing reads of the same allele will have random sequencing errors distributed independently. A consensus approach:
- **Reduces error rate**: Single error appearing in 1/10 reads gets outweighed by correct base in 9/10 reads
- **Stabilizes comparisons**: Two alleles can be compared reliably without error-induced variation
- **Enables somatic detection**: Real somatic mutations appear consistently in multiple reads; errors appear randomly
- **Mathematical principle**: With n reads and error rate p, consensus error rate ≈ p^(n/2), drastically lower for n≥3

---

### Step 2: Filter for Valid Heterozygotes
Apply sequential filters to identify individuals with two distinct alleles that can serve as reliable reference points for detecting somatic variants.

#### 2a. Group by Individual and Repeat Length
**Code Expression**:
```r
tab_sum_temp = tab %>% group_by(V1, V3) %>% 
  mutate(nREADs = n(),
         midSEG = names(sort(table(as.character(V4)), decreasing=T)[1]),
         midSEGlen = stringr::str_length(names(sort(table(as.character(V4)), decreasing=T)[1])),
         SEQ = consensus(SEQUENCES= as.character(V5), QUALITIES= as.character(V6), START= STARTstr),
         QUAL = ifelse(n() == 1, as.character(V6), "Consensus"))
```

**Detailed Explanation**:
- **`group_by(V1, V3)`**: Groups reads by individual (V1) and repeat length in bp (V3)
  - **Why this grouping?**: Reads with the same length likely come from the same allele. Individuals with two different lengths are likely heterozygotes
  - **Example**: Individual "EE10807Y" with reads of length 54bp grouped separately from reads of length 57bp
  
- **`nREADs = n()`**: Count number of reads in each group
  - **What it means**: How many sequencing reads support this allele variant
  - **Later use**: Identify the most-supported allele as the "main" germline allele
  
- **`midSEG = names(sort(table(as.character(V4)), decreasing=T)[1])`**: Find the most frequent repeat sequence
  - **Step-by-step**:
    - `table(as.character(V4))`: Creates frequency table of distinct sequences (V4 field)
    - `sort(..., decreasing=T)`: Orders by frequency (highest first)
    - `[1]`: Takes the most frequent
    - `names(...)`: Extracts the sequence string
  - **Why?**: Even reads of the same length may have different sequences due to point mutations. The majority sequence is consensus for that allele
  - **Example**: 80 reads with "AGCAGCAGC", 5 reads with "AGCAGCTGC" → midSEG = "AGCAGCAGC"
  
- **`SEQ = consensus(...)`**: Builds consensus from all reads in this group
  - **Calls the Step 1 function**: Uses quality-filtered voting to generate the definitive sequence
  - **Output**: High-confidence sequence string

#### 2b. Require Multiple Alleles (Heterozygotes Only)
**Code Expression**:
```r
tab_sum_temp = tab_sum_temp %>% 
  group_by(V1) %>% 
  mutate(nALLELES = length(unique(midSEGlen))) %>% 
  filter(nALLELES >= 2)
```

**Detailed Explanation**:
- **`group_by(V1)`**: Regroup at individual level (collapse across alleles)
- **`nALLELES = length(unique(midSEGlen))`**: Count how many distinct repeat lengths exist for this individual
  - **What it means**: Number of distinct alleles
  - **Example**: Individual with reads of 54bp and 57bp has nALLELES=2
  
- **`filter(nALLELES >= 2)`**: Keep only individuals with ≥2 distinct allele lengths
  - **Why this filter?**: 
    - Homozygotes (only one allele length) can't be used as reference because we can't distinguish germline from somatic
    - Heterozygotes provide two reference alleles, making it possible to detect variants
    - **Biological motivation**: Somatic mutations are rare; if you see two different lengths, the longer one is germline, not somatic

#### 2c. Require Minimum Read Support
**Code Expression**:
```r
tab_sum_temp = tab_sum_temp %>% 
  group_by(V1, nALLELES) %>% 
  mutate(A1len = midSEGlen[order(nREADs, decreasing=T)][1], 
         A2len = midSEGlen[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1]) %>% 
  group_by(V1) %>%
  filter(nREADs[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1] >= 3 & 
         sum(nREADs[match(unique(midSEGlen), midSEGlen)] > 2) >= 2)
```

**Detailed Explanation**:
- **`A1len = midSEGlen[order(nREADs, decreasing=T)][1]`**: Identify the most-supported allele
  - `order(nREADs, decreasing=T)`: Ranks alleles by read count (highest first)
  - `[1]`: Takes the most-supported allele length
  - **Biological interpretation**: The main germline allele is the one with most reads
  
- **`A2len = ... != A1len`**: Identify the second allele
  - Takes the most-supported allele that is NOT A1len
  - **Why?**: The second-most frequent allele in the sample
  
- **`filter(...[midSEGlen != A1len][1] >= 3)`**: Require second allele to have ≥3 reads
  - **Why 3 reads minimum?**: 
    - Balances sensitivity and specificity
    - 1 read could be a random error
    - 3+ independent reads confirms real allele variant
    - Also allows consensus building to work reliably (majority voting needs votes)
  
- **`filter(sum(nREADs[match(...)] > 2) >= 2)`**: Require at least 2 alleles with >2 reads each
  - **What it checks**: Both alleles must have meaningful support (not just the second allele)

#### 2d. Require Meaningful Difference Between Alleles
**Code Expression**:
```r
filter( abs(unique(A1len) - unique(A2len)) >= 5*repLen )
```

**Detailed Explanation**:
- **`5*repLen`**: Minimum difference is 5 repeat units (in bp)
  - **Example**: For trinucleotide repeats (repLen=3), need 5×3=15bp difference
  - If A1len=54bp, A2len must be ≤39bp or ≥69bp
  
- **Why this threshold?**:
  - **Excludes noise**: Sequencing/alignment errors can cause ±1-2bp changes, but ≥5 units is unambiguous
  - **Biological validity**: Different repeat lengths of >5 units typically represent different founder haplotypes
  - **Somatic event scale**: Somatic mutations scope (±1-2 units) is much smaller than this threshold
  
- **Output**: Only individuals with genuinely distinct alleles remain

**Why Step 2 is Critical**:
This step ensures that every individual analyzed:
1. **Has two alleles** → Can distinguish germline from somatic
2. **Has sufficient read support** → Consensus sequences are reliable
3. **Has meaningful differences** → Real genetic variants, not artifacts

**Impact**: Highly stringent filtering ensures high-confidence germline baseline for detecting somatic variants

---

### Step 3: Identify Alleles and Consensus Sequences
For each valid individual, determine the two germline alleles that will serve as reference.

#### Code Expression:
```r
tab_sum_temp = tab_sum_temp %>% 
  group_by(V1, nALLELES) %>% 
  mutate(A1len = midSEGlen[order(nREADs, decreasing=T)][1], 
         A2len = midSEGlen[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1],
         A1consensus = SEQ[order(nREADs, decreasing=T)][1], 
         A2consensus = SEQ[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1],
         A1midconsensus = midSEG[order(nREADs, decreasing=T)][1], 
         A2midconsensus = midSEG[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1])
```

#### Detailed Breakdown:

**3a. Select Most-Supported Allele as A1**
```r
A1len = midSEGlen[order(nREADs, decreasing=T)][1]
A1consensus = SEQ[order(nREADs, decreasing=T)][1]
A1midconsensus = midSEG[order(nREADs, decreasing=T)][1]
```
- **What it does**: Ranks alleles by read count and takes the top one
- **Why**: The most frequent allele is the main germline allele (not a somatic variant)
- **Stored information**:
  - **A1len**: Length in basepairs (e.g., 57)
  - **A1consensus**: Full consensus sequence (e.g., "AGGAGGAGCAGCAGC...")
  - **A1midconsensus**: Just the repeat motif (e.g., "AGCAGCAGC")

**3b. Select Second-Most-Supported Allele as A2**
```r
A2len = midSEGlen[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1]
A2consensus = SEQ[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1]
A2midconsensus = midSEG[order(nREADs, decreasing=T)][midSEGlen[order(nREADs, decreasing=T)] != A1len][1]
```
- **How it works**:
  - `[... != A1len]` excludes the A1 allele
  - Takes the most-supported remaining allele
- **Biological meaning**: The second germline allele (from the other chromosome)

**Summary of Step 3**:
After this step, each individual has stored:
- Two reference allele lengths
- Two reference consensus sequences (full and motif-only)
- These serve as the "germline ground truth" for somatic detection

---

### Step 4: Identify Potentially Somatic Reads
Filter for reads that deviate from the two germline alleles, suggesting they might represent somatic mutations.

#### Code Expression:
```r
tab_sum = tab_sum_temp %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(SOM_POT = ifelse(abs(midSEGlen - A1len) <= 2*repLen | 
                           abs(midSEGlen - A2len) <= 2*repLen, 1, 0)) %>% 
  filter(SOM_POT == 1)
```

#### Detailed Explanation:

**4a. Calculate Distance to Closest Allele**
```r
abs(midSEGlen - A1len)  # Distance to allele 1 in bp
abs(midSEGlen - A2len)  # Distance to allele 2 in bp
```
- **What it means**: How many basepairs away is this read from each germline allele?
- **Example**: 
  - A1len = 54bp, A2len = 57bp
  - Read with midSEGlen = 56bp is 2bp from A2, 3bp from A1
  - In repeat units: 56-57=-1 unit from A2, 56-54=2 units from A1

**4b. Apply Distance Threshold**
```r
abs(midSEGlen - A1len) <= 2*repLen | abs(midSEGlen - A2len) <= 2*repLen
```
- **Conversion to repeat units**: `2*repLen` bp = 2 repeat units (±2)
- **Meaning**: Keep reads within ±2 repeat units of ANY germline allele
- **Example with repLen=3**:
  - A1=54bp, A2=57bp
  - Accept reads with lengths: 48-60bp (±6bp from A1) OR 51-63bp (±6bp from A2)
  - Reject: 45bp or 66bp (these are >±2 units away)

**4c. Mark and Filter**
```r
SOM_POT = ifelse(..., 1, 0)
filter(SOM_POT == 1)
```
- **SOM_POT**: "Somatic potential" flag (1=yes, possibly somatic; 0=no, too different)
- **Why filter here?**: Reads >±2 units away are almost certainly artifacts, not somatic mutations
  - Could be misalignments or unrelated sequence variants
  - Somatic mutations occur as single-step changes (±1 or ±2 units max)

**Why This Threshold (±2 repeat units)?**
- **Sensitivity**: Catches real somatic events up to ±2 units
- **Specificity**: Excludes likely errors that are >2 units away
- **Biological basis**: Slippage during somatic DNA replication typically adds/removes 1-2 repeat units
- **Buffer for noise**: Some reads may have small misalignments or technical errors within ±2 units

**Output**: Only reads with SOM_POT=1 proceed. These are "candidate somatic reads" needing further validation.

---

### Step 5: Determine Origin Allele and Calculate Jump
For each candidate somatic read, determine which germline allele it likely originated from and by how much it differs.

#### Code Expression:
```r
tab_sum = tab_sum %>% 
  rowwise() %>% 
  mutate(fromLEN = ifelse(midSEGlen == A1len | midSEGlen == A2len, NA, 
                          c(A1len, A2len)[which.min(c(abs(midSEGlen - A1len), 
                                                       abs(midSEGlen - A2len)))]),
         fromSEQ = ifelse(midSEGlen == A1len, A1consensus,
                         ifelse(midSEGlen == A2len, A2consensus,
                               c(as.character(A1consensus), as.character(A2consensus))
                               [which.min(c(abs(midSEGlen - A1len), 
                                           abs(midSEGlen - A2len)))])),
         jump = ifelse(is.na(fromLEN), NA, (midSEGlen - fromLEN)/repLen)) %>% 
  filter(jump %in% c(-2, -1, 1, 2, NA))
```

#### Detailed Breakdown:

**5a. Identify if Read is a Main Allele**
```r
ifelse(midSEGlen == A1len | midSEGlen == A2len, NA, ...)
```
- **What it does**: Checks if the read's length exactly matches A1len or A2len
- **If TRUE**: `fromLEN = NA` (this is a germline allele, not a somatic variant)
- **If FALSE**: Calculate which allele this read came from
- **Why?**: Reads matching a germline allele serve as "controls" for quality assessment

**5b. Find Closest Allele (only for non-main alleles)**
```r
c(A1len, A2len)[which.min(c(abs(midSEGlen - A1len), abs(midSEGlen - A2len)))]
```
- **What it does**:
  - Calculates distance to A1: `abs(midSEGlen - A1len)` 
  - Calculates distance to A2: `abs(midSEGlen - A2len)`
  - `which.min(...)` returns the index of the shorter distance (1 or 2)
  - Selects that allele length
- **Example**:
  - Read is 55bp, A1=54bp, A2=57bp
  - Distance to A1: |55-54|=1bp
  - Distance to A2: |55-57|=2bp
  - Closest is A1, so fromLEN=54bp
- **Biological meaning**: This read likely expanded/contracted from A1

**5c. Extract Origin Sequence**
```r
fromSEQ = ifelse(midSEGlen == A1len, A1consensus,
                ifelse(midSEGlen == A2len, A2consensus, ...))
```
- **What it does**: Returns the full consensus sequence from whichever allele the read came from
- **Why?**: Needed for detailed comparison at mismatch sites later (Step 6)
- **Example**: If fromLEN=54bp (A1), then fromSEQ = A1consensus = "AGGAGGAGCAGCAGC..."

**5d. Calculate Jump in Repeat Units**
```r
jump = (midSEGlen - fromLEN) / repLen
```
- **What it does**: Converts length difference to number of repeat units
- **Formula breakdown**:
  - `midSEGlen - fromLEN`: Length difference in basepairs
  - Divide by `repLen` (3 for AGC): Convert bp to repeat units
- **Examples** (repLen=3):
  - Read is 56bp, fromLEN=54bp: jump = (56-54)/3 = 0.67 ≈ invalid
  - Read is 57bp, fromLEN=54bp: jump = (57-54)/3 = 1 (one expansion)
  - Read is 51bp, fromLEN=54bp: jump = (51-54)/3 = -1 (one contraction)
  - Read is 48bp, fromLEN=54bp: jump = (48-54)/3 = -2 (two contractions)

**5e. Validate Jump Value**
```r
filter(jump %in% c(-2, -1, 1, 2, NA))
```
- **What it does**: Only keeps reads with jump values of exactly -2, -1, +1, +2, or NA
- **Rejects**: Any read with fractional jumps (e.g., 0.67, 1.33)
- **Why?**: 
  - In-frame jumps (even units) are far more consistent with a Mendelian somatic event
  - Fractional jumps suggest the read didn't cleanly expand/contract and may be a technical artifact
  - This filters out noise while keeping biologically plausible variants

**Summary of Step 5**:
After this step, each read is characterized by:
- **fromLEN**: Which germline allele it originated from
- **fromSEQ**: The consensus sequence of that allele
- **jump**: How many repeat units it differs (-2, -1, 1, 2, or NA for main alleles)

These values enable the critical quality assessment in Step 6.

---

### Step 6: Quality Assessment - `mismatch_bases()` Function
**Core Purpose**: Distinguish true somatic mutations from sequencing/alignment errors using multi-layer quality metrics.

This is the most sophisticated part of the algorithm. For each candidate somatic read, it evaluates four hypothetical jump scenarios (-2, -1, +1, +2) and looks for evidence that the read genuinely represents a somatic mutation rather than an error.

#### 6a. Locate Flanks in the Read
**Code Expression**:
```r
locs_flanks = as.data.frame(stringr::str_locate(som_seq, c(start, end)))
if(nrow(locs_flanks) != 2 | sum(is.na(locs_flanks)) > 0){
  # Flanks not found - likely a misaligned read
  L_FLANK_GOOD = NA; R_FLANK_GOOD = NA; ...[return NAs]
}
```
- **What it does**: 
  - Searches for STARTstr in the read sequence (`som_seq`)
  - Searches for ENDstr in the read sequence
  - Returns start/end positions of each match
- **Why it matters**: 
  - Correct alignment means both flanks are present in expected positions
  - Missing flanks = likely misalignment or extracted sequence is wrong
  - If flanks missing, mark entire read as unreliable (return NAs)
- **Example**:
  - Read sequence: "AGGAGGA[AGC AGC AGC]AGCATGAAA"
  - locs_flanks row 1 finds "AGGAGGAGC" at positions 1-9
  - locs_flanks row 2 finds "AGCATGAAA" at positions 16-24

#### 6b. Calculate Flank Matching Quality
**Code Expression**:
```r
L_FLANK = stringr::str_split(som_qual,"")[[1]][locs_flanks[1,1]:locs_flanks[1,2]]
R_FLANK = stringr::str_split(som_qual,"")[[1]][locs_flanks[2,1]:locs_flanks[2,2]]

L_FLANK_GOOD = mean(unlist(lapply(L_FLANK, FUN = function(x) utf8ToInt(x)-33)) >= 25)
R_FLANK_GOOD = mean(unlist(lapply(R_FLANK, FUN = function(x) utf8ToInt(x)-33)) >= 25)
```

**Detailed Explanation**:
- **`L_FLANK = stringr::str_split(som_qual,"")[positions]`**: 
  - Extracts quality scores for all bases in the left flank
  - `som_qual` is the ASCII quality string, split into individual characters
  - Subset positions correspond to where STARTstr was found
- **`utf8ToInt(x) - 33`**: 
  - Converts ASCII quality character to PHRED score
  - 33 is the ASCII offset (PHRED score 0 = ASCII 33)
  - PHRED score 25 = ASCII 58 (represented as "@" or ":")
- **`>= 25`**: 
  - Filters to high-quality bases only
  - PHRED ≥25 means ≤0.3% error probability
- **`mean(...)`**: 
  - Calculates the proportion of high-quality bases
  - L_FLANK_GOOD ranges from 0 to 1
  - 1.0 = all flank bases are high quality (good)
  - 0.5 = only half the bases are high quality (mediocre)
  - 0.0 = all bases are low quality (bad)

**Why This Matters**:
- **Evidence of correct alignment**: Flanks should be perfectly conserved (germline sequences) and very high quality
- **If flanks are low quality**: The flank sequence might be misread, meaning the alignment is uncertain
- **Biological validation**: True somatic mutations have accurate flanks; misaligned reads often have corrupted flanks

**Also calculated**:
```r
L_FLANK_SUB_GOOD = mean(unlist(lapply(L_FLANK_SUB, ...)) >= 25)
R_FLANK_SUB_GOOD = mean(unlist(lapply(R_FLANK_SUB, ...)) >= 25)
```
- **What is L_FLANK_SUB?**: Quality scores at "key" positions only (positions specified in L_flank_set)
- **Why subset?**: Some bases are more critical than others for counting repeat units
  - Key positions are where a sequencing error could cause you to miscount repeats
  - Example: In "AGCAGC", the base at position 2 (G) is critical; if read as C, you'd miss a repeat unit
- **Subset quality is more stringent**: It's OK if some flank bases are low quality, but the critical ones must be good

#### 6c. Calculate Mismatch Rate at Sequence Boundaries
**Code Expression** (for forward strand):
```r
locs_S = as.data.frame(stringr::str_locate(seq, start))
locs_S_som = as.data.frame(stringr::str_locate(som_seq, start))
ref = rev(stringr::str_split(seq,"")[[1]][1:locs_S[1,1]])
compar = rev(stringr::str_split(som_seq,"")[[1]][1:locs_S_som[1,1]])
compar_qual = rev(stringr::str_split(som_qual,"")[[1]][1:locs_S_som[1,1]])
set = pmin(length(ref), length(compar))
compar = compar[1:set]; ref = ref[1:set]; compar_qual = compar_qual[1:set]
compar_qual_values = unlist(lapply(compar_qual, FUN = function(x) utf8ToInt(x)-33))

START_MISMATCH_RATE = mean(compar[compar_qual_values >= 25 & ref != "X"] != 
                            ref[compar_qual_values >= 25 & ref != "X"])
START_MISMATCH_BP_CHECKED = sum(compar_qual_values >= 25 & ref != "X")
```

**Step-by-Step Breakdown**:
- **Find flank positions in both reads**: `locs_S` in reference, `locs_S_som` in somatic read
- **Extract sequence upstream of flank**: Should be identical if correctly aligned
  - For forward strand: take everything BEFORE the start flank
  - `rev(...)` reverses the sequences so we can compare base-by-base
- **Align lengths**: `set = pmin(...)` ensures we only compare bases that exist in both
- **Filter to high-quality bases**: `compar_qual_values >= 25` 
- **Calculate mismatch rate**:
  ```
  mismatches = count where somatic_read[i] != reference[i]
  high_quality_bases = count where quality >= 25
  START_MISMATCH_RATE = mismatches / high_quality_bases
  ```

**What This Detects**:
- **Good alignment**: START_MISMATCH_RATE ≈ 0 (reads match reference before repeat)
- **Misalignment**: START_MISMATCH_RATE > 0.05 (too many mismatches)
- **Diagnostic**: If this somatic read is actually a complete misalignment, you'd see high mismatch rates at the boundary

#### 6d. Predict Mismatch Regions Based on Jump Type
**Code Expression** (example for +1 expansion):
```r
# For forward strand
locs = as.data.frame(stringr::str_locate(seq, end))
end_len = stringr::str_length(seq)
sub_seq_none = stringr::str_split(seq,"")[[1]][(locs[1,1]+repLen):end_len]      # Flank with NO repeat
sub_seq_one = stringr::str_split(seq,"")[[1]][(locs[1,1]):end_len]             # Flank with ONE repeat
sub_seq_two = stringr::str_split(seq,"")[[1]][(locs[1,1]-repLen):end_len]      # Flank with TWO repeats

locs_som = as.data.frame(stringr::str_locate(som_seq, end))
som_end = stringr::str_length(som_qual)

# If jump == 1 (one expansion), predict which bases would mismatch
pred_mismatch = which(sub_seq_one[1:length(sub_seq_none)] != sub_seq_none & 
                      sub_seq_none != "X" & sub_seq_one[1:length(sub_seq_none)] != "X")
pred_match = which(sub_seq_one[1:length(sub_seq_none)] == sub_seq_none & 
                   sub_seq_none != "X" & sub_seq_one[1:length(sub_seq_none)] != "X")

sub_som_qual = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]):som_end][pred_mismatch]
sub_som_qual_hi = stringr::str_split(som_qual,"")[[1]][(locs_som[1,1]):som_end][pred_match]
```

**Detailed Explanation**:
- **Three reference sequences with different repeat counts**:
  - `sub_seq_none`: Sequence if NO repeat units (zero repeats)
  - `sub_seq_one`: Sequence if ONE repeat unit present
  - `sub_seq_two`: Sequence if TWO repeat units present
- **Compare for +1 expansion**:
  - Expected pattern: The somatic read should match `sub_seq_one` if truly expanded from the original
  - Identify positions where `sub_seq_one` ≠ `sub_seq_none` (these bases differ due to adding one repeat)
  - `pred_mismatch`: Positions where addition of one repeat changes the sequence
  - `pred_match`: Positions where the sequence stays the same
- **Extract quality scores**:
  - `sub_som_qual`: Quality scores at predicted mismatch sites
  - `sub_som_qual_hi`: Quality scores at predicted matched sites (negative control)

**Why This Design?**:
- **Specificity**: We EXPECT certain bases to mismatch if this is truly a +1 expansion
- **Quality assessment**: Are those expected mismatches at good-quality bases?
- **Negative control**: Positions that shouldn't change DO change = evidence against somatic

#### 6e. Calculate Quality Metrics at Mismatch Sites
**Code Expression**:
```r
sub_som_qual_values = unlist(lapply(sub_som_qual[!is.na(sub_som_qual)], 
                                     FUN = function(x) utf8ToInt(x)-33))
sub_som_qual_hi_values = unlist(lapply(sub_som_qual_hi[!is.na(sub_som_qual_hi)], 
                                        FUN = function(x) utf8ToInt(x)-33))

return(data.table::data.table(
  rateSTRT_mis = START_MISMATCH_RATE, 
  nSTRT_mis = START_MISMATCH_BP_CHECKED,
  flnk_l_hiP = L_FLANK_GOOD, 
  flnk_l_hiP_sub = L_FLANK_SUB_GOOD,
  flnk_r_hiP = R_FLANK_GOOD, 
  flnk_r_hiP_sub = R_FLANK_SUB_GOOD,  
  LEN = length(sub_som_qual[!is.na(sub_som_qual)]),
  Q = paste0(sub_som_qual[!is.na(sub_som_qual)], collapse=""),
  NONF_SCORE = mean(sub_som_qual_values < 30),
  HI_QUAL_Q = paste0(sub_som_qual_hi[!is.na(sub_som_qual_hi)], collapse=""),
  HI_QUAL_LEN = length(sub_som_qual_hi[!is.na(sub_som_qual_hi)]),
  NONF_HIQUAL = mean(sub_som_qual_hi_values < 30)
))
```

**Metrics Explained**:

| Metric | Calculation | Interpretation |
|--------|-------------|-----------------|
| `LEN` | Count of predicted mismatch bases | How many positions expected to differ? ≥4 is strong evidence |
| `NONF_SCORE` | Proportion of mismatch bases with Q<30 | What % are low quality? <0.2 is good (80%+ high quality) |
| `HI_QUAL_LEN` | Count of predicted match bases | How many positions should be identical? |
| `NONF_HIQUAL` | Proportion of match bases with Q<30 | Should be high quality too (negative control) |

**Key Principle**:
- **Mismatch sites** (pred_mismatch) should have HIGH quality if this is genuine
  - These bases MUST be read accurately to trust the expansion/contraction
  - Low quality here = unreliable
- **Match sites** (pred_match) should also have HIGH quality
  - If you can't read these accurately either, the whole read is suspect
  - Used as sanity check

#### 6f. Summary of mismatch_bases() Output
The function returns a data.table with 12 metrics that collectively assess:
1. **Alignment quality** (rateSTRT_mis, nSTRT_mis)
2. **Flank preservation** (flnk_l_hiP, flnk_l_hiP_sub, flnk_r_hiP, flnk_r_hiP_sub)
3. **Predicted mismatch site quality** (LEN, Q, NONF_SCORE)
4. **Predicted match site quality** (HI_QUAL_LEN, HI_QUAL_Q, NONF_HIQUAL)

---

### Step 7: Applied Four Jump-Type Analyses
For each candidate somatic read, call `mismatch_bases()` four separate times to assess all possible expansion/contraction scenarios.

#### Code Expression:
```r
tab_sum = tab_sum %>% 
  rowwise() %>% 
  mutate(j_neg2 = mismatch_bases(strand=as.character(V2), jump=-2,
                                som_seq=as.character(V5), som_qual=as.character(V6),
                                seq=as.character(fromSEQ), start=STARTstr, end=ENDstr,
                                L_flank_set=L_FLANK, R_flank_set=R_FLANK, repeatLength=repLen),
         j_neg1 = mismatch_bases(strand=as.character(V2), jump=-1, ...),
         j_pos1 = mismatch_bases(strand=as.character(V2), jump=1, ...),
         j_pos2 = mismatch_bases(strand=as.character(V2), jump=2, ...))
```

#### Detailed Explanation:

**7a. Why Call Four Times?**
- **j_neg2**: Hypothesis that read = origin_allele - 2 repeat units (deletion of 2)
- **j_neg1**: Hypothesis that read = origin_allele - 1 repeat unit (deletion of 1)
- **j_pos1**: Hypothesis that read = origin_allele + 1 repeat unit (insertion of 1)
- **j_pos2**: Hypothesis that read = origin_allele + 2 repeat units (insertion of 2)

**7b. Why Test All Hypotheses?**
- We don't know which jump (if any) is correct
- A read at length 55bp could represent a +1 from 54bp OR a -1 from 57bp
- For each read, test which scenario best fits the quality metrics
- Later (Step 8), decide which jump (if any) passes quality thresholds

**7c. Function Arguments Explained**:
```r
mismatch_bases(strand=as.character(V2),     # Read strand (forward or reverse)
               jump=-2,                      # Jump hypothesis being tested
               som_seq=as.character(V5),     # The somatic/variant read sequence
               som_qual=as.character(V6),    # Quality scores of that read
               seq=as.character(fromSEQ),    # Consensus sequence of origin allele
               start=STARTstr,               # Left flank sequence
               end=ENDstr,                   # Right flank sequence
               L_flank_set=L_FLANK,          # Key positions in left flank
               R_flank_set=R_FLANK,          # Key positions in right flank
               repeatLength=repLen)          # Repeat unit length (usually 3)
```

**7d. Output Structure**:
Each call returns a data.table with 12 quality metrics (from Step 6e):
- `j_neg2`: Metrics IF this read represents -2 jump
- `j_neg1`: Metrics IF this read represents -1 jump  
- `j_pos1`: Metrics IF this read represents +1 jump
- `j_pos2`: Metrics IF this read represents +2 jump

**Example Output** (for a single read):
| Metric | j_neg2 | j_neg1 | j_pos1 | j_pos2 |
|--------|--------|--------|--------|--------|
| LEN | 3 | 5 | 8 | 2 |
| NONF_SCORE | 0.33 | 0.2 | 0.1 | 0.5 |
| rateSTRT_mis | 0.08 | 0.02 | 0.01 | 0.06 |
| ... | ... | ... | ... | ... |

**Interpretation**:
- **j_pos1** looks best: LEN=8 good, NONF_SCORE=0.1 good, rateSTRT_mis=0.01 excellent
- **j_neg2** looks worst: LEN=3 small, NONF_SCORE=0.33 poor, rateSTRT_mis=0.08 bad
- This suggests if this read is somatic, it's a +1 expansion, not a -2 contraction

---

### Step 8: Filter for True Somatic Reads
Apply stringent quality thresholds to the jump hypotheses, separating real somatic mutations from artifacts.

#### Code Expression:
```r
tab_sum = tab_sum %>%
  rowwise() %>% 
  mutate(somatic_neg2 = ifelse(j_neg2$LEN >= 4 & j_neg2$NONF_SCORE < 0.2 & 
                               j_neg2$rateSTRT_mis < 0.05 & 
                               ((j_neg2$flnk_l_hiP_sub == 1 | is.na(j_neg2$flnk_l_hiP_sub)) & 
                                (j_neg2$flnk_r_hiP_sub == 1 | is.na(j_neg2$flnk_r_hiP_sub))), 1, 0),
         somatic_neg1 = ifelse(j_neg1$LEN >= 4 & j_neg1$NONF_SCORE < 0.2 & 
                               j_neg1$rateSTRT_mis < 0.05 & 
                               ((j_neg1$flnk_l_hiP_sub == 1 | is.na(j_neg1$flnk_l_hiP_sub)) & 
                                (j_neg1$flnk_r_hiP_sub == 1 | is.na(j_neg1$flnk_r_hiP_sub))), 1, 0),
         somatic_pos1 = ifelse(j_pos1$LEN >= 4 & j_pos1$NONF_SCORE < 0.2 & 
                               j_pos1$rateSTRT_mis < 0.05 & 
                               ((j_pos1$flnk_l_hiP_sub == 1 | is.na(j_pos1$flnk_l_hiP_sub)) & 
                                (j_pos1$flnk_r_hiP_sub == 1 | is.na(j_pos1$flnk_r_hiP_sub))), 1, 0),
         somatic_pos2 = ifelse(j_pos2$LEN >= 4 & j_pos2$NONF_SCORE < 0.2 & 
                               j_pos2$rateSTRT_mis < 0.05 & 
                               ((j_pos2$flnk_l_hiP_sub == 1 | is.na(j_pos2$flnk_l_hiP_sub)) & 
                                (j_pos2$flnk_r_hiP_sub == 1 | is.na(j_pos2$flnk_r_hiP_sub))), 1, 0))
```

#### Detailed Breakdown:

**8a. Quality Threshold 1: Depth of Mismatch Evidence**
```r
j_[jump]$LEN >= 4
```
- **What it means**: The read must have ≥4 basepairs at positions where we expect mismatches
- **Why?**: 
  - If `LEN < 4`, we have little evidence to assess quality
  - 4+ bases gives statistically meaningful signal
  - Arbitrary but reasonable threshold
  - Filters out reads where the mismatch region is too small to judge

**8b. Quality Threshold 2: Base Quality at Mismatch Sites**
```r
j_[jump]$NONF_SCORE < 0.2
```
- **What it means**: 
  - NONF_SCORE = proportion of bases with PHRED < 30 (lower quality)
  - PHRED < 30 means ~99.9% accuracy, still quite good
  - `< 0.2` = less than 20% of bases can be low quality
  - So at least 80% must have PHRED ≥ 30
- **Why this matters**:
  - Mismatch sites are CRITICAL - they define the somatic event
  - True mutations should be read at high quality
  - Sequencing errors tend to occur at low-quality positions
  - If mismatch sites are mostly low quality, probably an error
- **Interpretation**:
  - NONF_SCORE = 0.0: Perfect, all bases are high quality ✓ PASS
  - NONF_SCORE = 0.2: OK, 80% high quality ✓ PASS
  - NONF_SCORE = 0.3: Problematic, 70% high quality ✗ FAIL
  - NONF_SCORE = 1.0: All bases low quality ✗ FAIL

**8c. Quality Threshold 3: Boundary Alignment Accuracy**
```r
j_[jump]$rateSTRT_mis < 0.05
```
- **What it means**:
  - START_MISMATCH_RATE = proportion of bases that mismatch at sequence boundary
  - `< 0.05` = less than 5% mismatch rate
  - So at least 95% of boundary bases must match the reference
- **Why it matters**:
  - The boundary between read and reference should be nearly perfect
  - If misaligned, boundaries will have many mismatches
  - Good alignment = rateSTRT_mis ≈ 0.0
  - Misalignment = rateSTRT_mis could be 0.2-0.5 or higher
- **Interpretation**:
  - rateSTRT_mis = 0.01: Excellent alignment ✓ PASS
  - rateSTRT_mis = 0.04: Good alignment ✓ PASS
  - rateSTRT_mis = 0.10: Possible misalignment ✗ FAIL

**8d. Quality Threshold 4: Flank Key Base Accuracy**
```r
((j_[jump]$flnk_l_hiP_sub == 1 | is.na(j_[jump]$flnk_l_hiP_sub)) & 
 (j_[jump]$flnk_r_hiP_sub == 1 | is.na(j_[jump]$flnk_r_hiP_sub)))
```

**Detailed Explanation**:
- **`flnk_l_hiP_sub == 1`**: 
  - All key left flank bases have PHRED ≥ 25
  - Perfect score = 1.0 (means 100% of key bases are high quality)
  - `== 1` means: Either exactly 1.0 (all high quality) OR NA (no key bases to check)
  - NA is acceptable because some loci don't have key positions
  
- **`flnk_r_hiP_sub == 1`**:
  - Same logic for right flank key bases
  - Either 1.0 (perfect) or NA (not applicable)

- **`&` connective**:
  - Both flank conditions must be true
  - Left flank AND right flank must both be perfect (or N/A)
  - Even one flank with quality <1.0 fails this test

**Why Key Bases Matter**:
- Some positions in the flank are "critical" for counting repeats correctly
- A single base error at a key position can cause miscounting
- Example: In "AGC", the "G" at position 2 might be critical
  - If read as "A**C**C", you'd think it's "ACC" (wrong repeat)
  - vs. correct "AGC" (right repeat)
- These key bases MUST be high quality, with no tolerance for error

**8e. Integrated Logic**
All four conditions must be true for `somatic_[jump] = 1`:
```
somatic_[jump] = 1  if  (LEN >= 4)  AND  (NONF_SCORE < 0.2)  AND  
                         (rateSTRT_mis < 0.05)  AND  (flank quality perfect)
somatic_[jump] = 0  otherwise
```

**Decision Logic**:
- A read only gets `somatic_[jump] = 1` if it passes ALL four tests
- Failing any single test → `somatic_[jump] = 0`
- This is **AND logic**, not OR logic (highly stringent)

**8f. Interpretation of Final Output**
**For main allele reads** (jump=NA):
- `somatic_neg2, somatic_neg1, somatic_pos1, somatic_pos2` indicate how well this read would score if it WERE a somatic event
- Used as controls: Main allele reads should mostly have somatic_*=0 (because they're not actually somatic)
- Any with somatic_*=1 are outliers (possibly misclassified reads)

**For candidate somatic reads** (jump ∈ {-2,-1,1,2}):
- Only the matching somatic_* column is relevant
  - E.g., if jump=1, look at somatic_pos1 only
  - somatic_pos1=1 means "this +1 expansion passes all quality checks"
  - somatic_pos1=0 means "this read fails quality check; probably not a real somatic mutation"

**Example Interpretation**:
```
V1="EE10807Y", midSEGlen=56, fromLEN=54, jump=1
somatic_neg2=0, somatic_neg1=0, somatic_pos1=1, somatic_pos2=0
```
- Read differs by +1 repeat unit from the 54bp allele
- Only the `somatic_pos1=1` test passes
- Interpretation: "This read likely represents a true +1 somatic expansion"

```
V1="EE10807Y", midSEGlen=56, fromLEN=54, jump=1
somatic_neg2=0, somatic_neg1=0, somatic_pos1=0, somatic_pos2=0
```
- Same read, but ALL quality tests fail
- Interpretation: "This read doesn't pass quality checks; likely a sequencing error or alignment artifact"

**8g. Summary of Quality Control Philosophy**

The algorithm uses a **cascade of increasingly specific checks**:

1. **Bulk filtering** (earlier steps): Remove obvious noise (reads >±2 units away)
2. **Evidence depth** (LEN ≥ 4): Ensure enough signal to assess quality
3. **Base quality** (NONF_SCORE < 0.2): Only high-confidence bases at key sites
4. **Alignment quality** (rateSTRT_mis < 0.05): Confirm correct positioning
5. **Critical base accuracy** (flank_hiP_sub = 1): Zero tolerance for key base errors

**Biological Rationale**:
- Somatic mutations are real DNA changes, sequenced as above-normal coverage reads
- These should have high quality at mismatch sites (proven sequencing accuracy)
- Errors are typically low-quality bases or misalignments
- True somatic variants won't cluster at low-quality positions
- This layered approach makes false positives highly unlikely

---

## Output File: `summary_somatic_{ID}.txt`

### Output Columns

| Column | Type | Description |
|--------|------|-------------|
| **V1** | Character | Individual identifier |
| **V3** | Numeric | Length of this particular read's repeat (bp) |
| **V4** | Character | Repeat sequence of this read |
| **jump** | Numeric | Difference from origin allele in repeat units (-2, -1, 1, 2, NA) |
| **fromLEN** | Numeric | Length of allele this read originated from (bp) |
| **A1len** | Numeric | Length of germline allele 1 (bp) |
| **A1midconsensus** | Character | Repeat motif sequence of allele 1 |
| **A2len** | Numeric | Length of germline allele 2 (bp) |
| **A2midconsensus** | Character | Repeat motif sequence of allele 2 |
| **somatic_neg2** | Binary (0/1) | Does this read pass somatic filters for -2 jump? |
| **somatic_neg1** | Binary (0/1) | Does this read pass somatic filters for -1 jump? |
| **somatic_pos1** | Binary (0/1) | Does this read pass somatic filters for +1 jump? |
| **somatic_pos2** | Binary (0/1) | Does this read pass somatic filters for +2 jump? |

### Interpretation

Each **row** represents a single read from the CRAM file.

- **Main allele reads**: `jump = NA`, columns somatic_* indicate how they'd score under hypothetical jumps
- **Potentially somatic reads**: `jump ∈ {-2,-1,1,2}`, relevant somatic_* column indicates if it passes quality filters
- **Somatic probability**: A read with `jump = 1` and `somatic_pos1 = 1` is likely a true +1 expansion event (assuming it's real and not an artifact)

---

## Downstream Analysis (compile_summary.R)

The output file feeds into `compile_summary.R` which:

1. **Counts allele frequencies**: How many individuals carry each allele variant
   - Groups by allele length (A1len, A2len) and repeat sequence (A1midconsensus, A2midconsensus)
   - Counts unique individuals carrying each combination
   
2. **Counts somatic reads**: How many reads passed somatic filters for each jump type
   - Filters rows where jump ≠ NA (somatic candidates only)
   - For each jump value (-2, -1, 1, 2), sums the corresponding somatic_* column
   - Results grouped by jump type, origin allele length, and origin sequence
   
3. **Calculates expected rates**: Using reads without jumps as a negative control
   - For main allele reads (jump=NA), applies the same somatic filters
   - Calculates expected frequency of "false positive" somatic calls
   - Acts as denominator for rate calculations
   
4. **Computes somatic expansion rates**: Somatic events per individual per allele per read
   - Formula: (number of somatic reads passing filters) / (expected denominator based on main alleles)
   - Normalizes by allele frequency and number of individuals
   - Produces **somatic expansion rates** for each allele variant
   - Rates typically range from 0.0001 to 0.02 (0.01% to 2% per read)

**Output**: Table of somatic mutation rates by allele length and sequence, enabling:
- Comparison of mutation rates across different repeat lengths
- Identification of allele effects (some sequences/lengths mutate more frequently)
- Quantification of somatic burden in repeat expansion disorders

---

## Multi-Layer Quality Control Summary

The algorithm uses **five integrated layers of quality assessment**:

### Layer 1: Read-Level Quality (Step 1)
- PHRED quality score filtering (Q ≥ 25)
- Consensus building from multiple reads
- Effect: Noise reduction at base call level

### Layer 2: Individual-Level Validation (Step 2)
- Require heterozygous genotypes (≥2 alleles)
- Require minimum read support (≥3 reads per allele)
- Require large allele difference (≥5 repeat units)
- Effect: High-confidence germline baseline

### Layer 3: Candidate Selection (Steps 4-5)
- Restrict to ±2 repeat units from main alleles
- Calculate exact jump value in repeat units
- Require in-frame jumps (-2, -1, 1, 2)
- Effect: Plausible somatic event range only

### Layer 4: Sequence Quality (Step 6)
- Flank location verification
- Flank base quality assessment (Q ≥ 25)
- Critical flank position accuracy verification
- Boundary alignment accuracy (mismatch rate <5%)
- Effect: Confirm correct sequence alignment

### Layer 5: Mismatch Site Specificity (Steps 7-8)
- Four independent jump hypothesis testing
- Depth of evidence (≥4 affected bases)
- Quality at predicted mismatch sites (≥80% Q≥30)
- Stringent AND logic (all tests must pass)
- Effect: Only high-confidence differences flagged

**Net Effect**: A read must simultaneously score well on:
✓ Individual quality metrics
✓ Flank sequence preservation
✓ Alignment accuracy
✓ Mismatch site quality
✓ Biological plausibility

False positive probability becomes vanishingly small with all layers stacked.

---

## Key Algorithmic Insights

### 1. Why Multiple Jump Hypotheses?
- A read at intermediate length could arise from different origins
- Testing all four jumps independently avoids bias
- Quality metrics distinguish which jump (if any) is real

### 2. Why Quality Thresholds Are Stringent
- Somatic mutations are rare (0.001%-2% per read)
- Sequencing errors can appear at similar frequencies in bad runs
- Distinguish through quality patterns:
  - True mutations: High-quality bases at mismatch sites
  - Errors: Cluster at low-quality positions

### 3. Why Flanks Matter
- Flanks serve as "proof of correct alignment"
- If flanks are corrupted, alignment is uncertain
- If alignment is wrong, you can't trust the repeat length measurement
- Perfect flanks = confidence in length assignment

### 4. Why Consensus Sequences Are Built
- Multiple representations of same allele improve accuracy
- Reduces technical noise from single reads
- Enables reliable comparison between reads

### 5. Why Controls Are Necessary
- Main allele reads (jump=NA) provide negative control for quality thresholds
- If too many main alleles score as "somatic", thresholds are too loose
- This feedback ensures calibration against the actual error rate

---

## Common Interpretations

### Read with jump=NA (Main Allele)
- **Meaning**: This read matches one of the two germline alleles exactly
- **somatic_* columns**: Indicate how this read would score under hypothetical jumps
- **Expected**: Most should have somatic_*=0 (correctly identified as non-somatic)
- **Outliers**: Any with somatic_*=1 suggest possible misclassification

### Read with jump=+1, somatic_pos1=1
- **Meaning**: This read is 1 repeat unit longer than its origin allele and passes all quality tests
- **Interpretation**: Likely a true somatic expansion event
- **Frequency**: Rare (0.1%-2% depending on allele)
- **Biological significance**: Evidence of ongoing mutagenesis in that cell lineage

### Read with jump=+1, somatic_pos1=0, somatic_[other]=0
- **Meaning**: This read differs by +1 units but fails quality checks
- **Interpretation**: Sequencing error, misalignment, or technical artifact
- **Action**: Discard; not counted as somatic event

### Read with jump=1 but multiple somatic_*=1
- **Meaning**: Read passes quality checks for multiple jump hypotheses
- **Interpretation**: Ambiguous; likely borderline quality
- **Action**: May be included in somatic count if jump=1 specifically tested
- **Note**: Ideally want only a single somatic_*=1 (unambiguous)

---

## Final Notes

This algorithm represents a sophisticated multi-layered approach to distinguishing true somatic repeat mutations from technical artifacts. Each layer independently contributes to confidence, and only reads passing all layers are marked as somatic. The approach leverages both sequence quality metrics and biological expectations to achieve high specificity while maintaining reasonable sensitivity for detecting real somatic events.
