#!/usr/bin/env python3
"""
Python translation of short_somatic_perIndividual.R
Produces identical results to the R version.
"""

import sys
import re
import math
import csv
from collections import Counter

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
ID  = sys.argv[1]
DAT = sys.argv[2]

# ---------------------------------------------------------------------------
# Parse DAT (mirrors R: stringr::str_split(DAT,"_")[[1]])
# R uses 1-based indexing: [5]=idx4, [6]=idx5, [7]=idx6, [8]=idx7, [9]=idx8
# ---------------------------------------------------------------------------
parts = DAT.split("_")
STARTstr = parts[4]                                      # R: [5]
repLen   = len(parts[6])                                 # R: str_length([7])
ENDstr   = parts[5]                                      # R: [6]

def _parse_numeric_list(s):
    """Parse a comma-separated string (possibly containing 'NA') into a list.
    NA values are represented as float('nan')."""
    result = []
    for x in s.split(","):
        x = x.strip()
        if not x:
            continue
        if x.upper() == "NA":
            result.append(float("nan"))
        else:
            result.append(float(x))
    return result

# str_remove removes the first occurrence; equivalent to str.replace(..., 1)
L_FLANK = _parse_numeric_list(parts[7].replace("L", "", 1))   # R: str_remove([8],"L")
R_FLANK = _parse_numeric_list(parts[8].replace("R", "", 1))   # R: str_remove([9],"R")

# ---------------------------------------------------------------------------
# Load individual read data  (R: fread("IID_{ID}.txt", h=F, sep=" "))
# ---------------------------------------------------------------------------
tab = pd.read_csv(f"IID_{ID}.txt", header=None, sep=" ",
                  names=["V1","V2","V3","V4","V5","V6"],
                  dtype=str)
print(f"DEBUG: Loaded {len(tab)} reads for individual {ID}", file=sys.stderr)
print(f"DEBUG: Command line args - ID: {ID} DAT: {DAT}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Helper: str_locate  (mirrors stringr::str_locate)
# Returns (start, end) 1-indexed inclusive, or (None, None) if no match.
# Uses plain string search (re.escape) to match literal patterns.
# ---------------------------------------------------------------------------
def str_locate(string, pattern):
    m = re.search(re.escape(pattern), string)
    if m:
        return (m.start() + 1, m.end())   # 1-indexed
    return (None, None)


# ---------------------------------------------------------------------------
# consensus()
# Mirrors the R consensus() function exactly.
# Aligns sequences by STARTstr, accumulates bases/qualities per position
# (SPOT range -156..156), then at each position takes the most frequent
# base among those with ASCII(SEQ_char)-33 >= 25 (which is always true for
# A/C/G/T/N, replicating the R behaviour).
# ---------------------------------------------------------------------------
def consensus(sequences, qualities, start):
    if len(sequences) == 1:
        return sequences[0]

    # SEQdf: dict  SPOT -> accumulated SEQ chars, QUAL chars
    spots = list(range(-156, 157))       # seq(-156, 156, by=1)  => 313 values
    seq_acc  = {s: [] for s in spots}
    qual_acc = {s: [] for s in spots}

    start_locs = [str_locate(s, start) for s in sequences]
    seq_split  = [list(s) for s in sequences]
    qual_split = [list(q) for q in qualities]

    for i, (sl, sc, qc) in enumerate(zip(start_locs, seq_split, qual_split)):
        if sl[0] is None:
            continue
        _, sl_end = sl                      # 1-indexed end of match
        for j, (base, qual) in enumerate(zip(sc, qc)):
            spot = (j + 1) - sl_end - 2    # R: seq(1,len) - START_LOCS[i,2] - 2
            if -156 <= spot <= 156:
                seq_acc[spot].append(base)
                qual_acc[spot].append(qual)

    consensus_seq = []
    for s in spots:
        chars = seq_acc[s]
        if not chars:
            continue
        # R: opts = split(SEQdf$SEQ[i])[ utf8ToInt(SEQdf$SEQ[i]) - 33 >= 25 ]
        # utf8ToInt on SEQ chars (A/C/G/T) always >= 25+33=58 since ord('A')=65
        # so every base passes; we replicate that behaviour literally
        opts = [c for c in chars if ord(c) - 33 >= 25]
        if opts:
            counts     = Counter(opts)
            top        = counts.most_common(2)
            if len(top) == 1 or top[0][1] != top[1][1]:
                consensus_seq.append(top[0][0])
            else:
                consensus_seq.append("X")   # tie
        else:
            consensus_seq.append("X")

    return "".join(consensus_seq)


# ---------------------------------------------------------------------------
# mismatch_bases()
# Mirrors the R mismatch_bases() function exactly.
# Returns a dict with the same keys as the R data.table.
# ---------------------------------------------------------------------------
def _qual_values(chars):
    """ord(c) - 33 for a list of quality characters."""
    return [ord(c) - 33 for c in chars]


def _safe_mean(values):
    """mean() that returns NaN for empty input (matching R behaviour)."""
    if not values:
        return float("nan")
    return sum(values) / len(values)


def mismatch_bases(strand, jump, som_seq, som_qual,
                   seq, start, end,
                   L_flank_set, R_flank_set, repeat_length):

    NA_RESULT = dict(
        rateSTRT_mis=float("nan"), nSTRT_mis=float("nan"),
        flnk_l_hiP=float("nan"),  flnk_l_hiP_sub=float("nan"),
        flnk_r_hiP=float("nan"),  flnk_r_hiP_sub=float("nan"),
        LEN=float("nan"),         Q=float("nan"),
        NONF_SCORE=float("nan"),
        HI_QUAL_Q=float("nan"),   HI_QUAL_LEN=float("nan"),
        NONF_HIQUAL=float("nan"),
    )

    # ---- 1. Locate flanks in the somatic read --------------------------------
    ls_s, ls_e = str_locate(som_seq, start)
    le_s, le_e = str_locate(som_seq, end)

    if ls_s is None or le_s is None:
        return NA_RESULT

    # locs_flanks has 2 rows: row0 = start-flank location, row1 = end-flank location
    locs_flanks = [(ls_s, ls_e), (le_s, le_e)]

    som_qual_chars = list(som_qual)

    # Quality chars at left and right flanks (1-indexed inclusive → Python slice)
    L_FLANK_chars = som_qual_chars[locs_flanks[0][0] - 1 : locs_flanks[0][1]]
    R_FLANK_chars = som_qual_chars[locs_flanks[1][0] - 1 : locs_flanks[1][1]]

    # Sub-flanks at the pre-specified positions (remove NaN entries first)
    lf_idx = [int(x) - 1 for x in L_flank_set if not math.isnan(x)]   # 0-indexed
    rf_idx = [int(x) - 1 for x in R_flank_set if not math.isnan(x)]

    L_FLANK_SUB = [L_FLANK_chars[i] for i in lf_idx if i < len(L_FLANK_chars)]
    R_FLANK_SUB = [R_FLANK_chars[i] for i in rf_idx if i < len(R_FLANK_chars)]

    def _frac_hq(chars):
        if not chars:
            return float("nan")
        return _safe_mean([v >= 25 for v in _qual_values(chars)])

    L_FLANK_GOOD     = _frac_hq(L_FLANK_chars)
    R_FLANK_GOOD     = _frac_hq(R_FLANK_chars)
    L_FLANK_SUB_GOOD = _frac_hq(L_FLANK_SUB)
    R_FLANK_SUB_GOOD = _frac_hq(R_FLANK_SUB)

    # ---- 2. START_MISMATCH (mismatch rate at the "outside" part of the read) --
    seq_chars     = list(seq)
    som_seq_chars = list(som_seq)

    if strand == "reverse":
        le_seq_s, le_seq_e = str_locate(seq, end)
        le_som_s, le_som_e = str_locate(som_seq, end)
        if le_seq_s is None:
            START_MISMATCH_RATE = START_MISMATCH_BP_CHECKED = -9
        else:
            # R: [locs_E[1,2] : str_length(seq)]  (1-indexed, inclusive both ends)
            ref         = seq_chars[le_seq_e - 1:]
            compar      = som_seq_chars[le_som_e - 1:]
            compar_qual = som_qual_chars[le_som_e - 1:]
            set_len     = min(len(ref), len(compar))
            ref         = ref[:set_len]
            compar      = compar[:set_len]
            compar_qual = compar_qual[:set_len]
            cqv         = _qual_values(compar_qual)
            valid       = [i for i in range(set_len)
                           if cqv[i] >= 25 and ref[i] != "X"]
            if valid:
                START_MISMATCH_RATE       = sum(compar[i] != ref[i] for i in valid) / len(valid)
                START_MISMATCH_BP_CHECKED = len(valid)
            else:
                START_MISMATCH_RATE       = float("nan")
                START_MISMATCH_BP_CHECKED = 0
    else:  # forward
        ls_seq_s, ls_seq_e = str_locate(seq, start)
        if ls_seq_s is None:
            START_MISMATCH_RATE = START_MISMATCH_BP_CHECKED = -9
        else:
            ls_som_s, ls_som_e = str_locate(som_seq, start)
            # R: rev(split(seq)[1:locs_S[1,1]])  (1-indexed inclusive)
            ref         = seq_chars[:ls_seq_s][::-1]
            compar      = som_seq_chars[:ls_som_s][::-1]
            compar_qual = som_qual_chars[:ls_som_s][::-1]
            set_len     = min(len(ref), len(compar))
            ref         = ref[:set_len]
            compar      = compar[:set_len]
            compar_qual = compar_qual[:set_len]
            cqv         = _qual_values(compar_qual)
            valid       = [i for i in range(set_len)
                           if cqv[i] >= 25 and ref[i] != "X"]
            if valid:
                START_MISMATCH_RATE       = sum(compar[i] != ref[i] for i in valid) / len(valid)
                START_MISMATCH_BP_CHECKED = len(valid)
            else:
                START_MISMATCH_RATE       = float("nan")
                START_MISMATCH_BP_CHECKED = 0

    # ---- 3. Predicted mismatch positions per jump value ----------------------
    sub_som_qual          = None
    sub_som_qual_hi_chars = None

    if strand == "reverse":
        # Find START flank in reference sequence
        locs_s, locs_e = str_locate(seq, start)
        if locs_s is None:
            sub_som_qual          = [None]
            sub_som_qual_hi_chars = [None]
        else:
            # R: [1:(locs[1,2]-rl)]  1-indexed inclusive → Python [:locs_e - rl]
            sub_seq_none = seq_chars[:locs_e - repeat_length]
            sub_seq_one  = seq_chars[:locs_e]
            sub_seq_two  = seq_chars[:locs_e + repeat_length]

            locs_som_s, locs_som_e = str_locate(som_seq, start)

            def _rev_comparison(sA, sB, l):
                """Return (pred_mismatch_idx, pred_match_idx) for reversed arrays."""
                n = min(len(sA), len(sB), l) if l > 0 else min(len(sA), len(sB))
                rA = sA[:n][::-1]
                rB = sB[:n][::-1]
                pm  = [i for i in range(len(rA))
                       if rA[i] != rB[i] and rA[i] != "X" and rB[i] != "X"]
                phi = [i for i in range(len(rA))
                       if rA[i] == rB[i] and rA[i] != "X" and rB[i] != "X"]
                return pm, phi

            if jump < 0:
                if jump == -1:
                    # rev(sub_seq_none) vs rev(sub_seq_one[(rl+1):end])
                    part = sub_seq_one[repeat_length:]
                    pm, phi = _rev_comparison(sub_seq_none, part,
                                              min(len(sub_seq_none), len(part)))
                    sq_end = locs_som_e - repeat_length
                elif jump == -2:
                    part = sub_seq_two[2 * repeat_length:]
                    pm, phi = _rev_comparison(sub_seq_none, part,
                                              min(len(sub_seq_none), len(part)))
                    sq_end = locs_som_e - repeat_length
                else:
                    pm, phi, sq_end = [], [], locs_som_e - repeat_length

                rev_q = som_qual_chars[:sq_end][::-1]
                sub_som_qual          = [rev_q[i] for i in pm  if i < len(rev_q)]
                sub_som_qual_hi_chars = [rev_q[i] for i in phi if i < len(rev_q)]

            elif jump == 1:
                part = sub_seq_one[repeat_length:]
                pm, phi = _rev_comparison(part, sub_seq_none,
                                          min(len(sub_seq_none), len(part)))
                sq_end = locs_som_e
                rev_q = som_qual_chars[:sq_end][::-1]
                sub_som_qual          = [rev_q[i] for i in pm  if i < len(rev_q)]
                sub_som_qual_hi_chars = [rev_q[i] for i in phi if i < len(rev_q)]

            elif jump == 2:
                part = sub_seq_two[2 * repeat_length:]
                pm, phi = _rev_comparison(part, sub_seq_none,
                                          min(len(sub_seq_none), len(part)))
                sq_end = locs_som_e + repeat_length
                rev_q = som_qual_chars[:sq_end][::-1]
                sub_som_qual          = [rev_q[i] for i in pm  if i < len(rev_q)]
                sub_som_qual_hi_chars = [rev_q[i] for i in phi if i < len(rev_q)]
            else:
                sub_som_qual          = [None]
                sub_som_qual_hi_chars = [None]

    else:  # forward
        # Find END flank in reference sequence
        locs_s, locs_e = str_locate(seq, end)
        end_len = len(seq)
        if locs_s is None:
            sub_som_qual          = [None]
            sub_som_qual_hi_chars = [None]
        else:
            # R: [(locs[1,1]+rl):end_len]  1-indexed inclusive
            # Python: [locs_s - 1 + rl : end_len]  (locs_s is 1-indexed)
            sub_seq_none = seq_chars[locs_s - 1 + repeat_length:]
            sub_seq_one  = seq_chars[locs_s - 1:]
            sub_seq_two  = seq_chars[locs_s - 1 - repeat_length:]

            locs_som_s, locs_som_e = str_locate(som_seq, end)
            som_end = len(som_qual)

            def _fwd_comparison(sA, sB, n):
                """Return (pred_mismatch_idx, pred_match_idx) for forward arrays."""
                rA = sA[:n]
                rB = sB[:n]
                pm  = [i for i in range(len(rA))
                       if rA[i] != rB[i] and rA[i] != "X" and rB[i] != "X"]
                phi = [i for i in range(len(rA))
                       if rA[i] == rB[i] and rA[i] != "X" and rB[i] != "X"]
                return pm, phi

            if jump < 0:
                if jump == -1:
                    # sub_seq_none vs sub_seq_one[0:len(sub_seq_none)]
                    n = min(len(sub_seq_none), len(sub_seq_one))
                    pm, phi = _fwd_comparison(sub_seq_none, sub_seq_one, n)
                    sq_start = locs_som_s - 1 + repeat_length   # 0-indexed
                elif jump == -2:
                    n = min(len(sub_seq_none), len(sub_seq_two))
                    pm, phi = _fwd_comparison(sub_seq_none, sub_seq_two, n)
                    sq_start = locs_som_s - 1 + repeat_length
                else:
                    pm, phi, sq_start = [], [], locs_som_s - 1 + repeat_length

                q_slice = som_qual_chars[sq_start:som_end]
                sub_som_qual          = [q_slice[i] for i in pm  if i < len(q_slice)]
                sub_som_qual_hi_chars = [q_slice[i] for i in phi if i < len(q_slice)]

            elif jump == 1:
                n = min(len(sub_seq_none), len(sub_seq_one))
                pm, phi = _fwd_comparison(sub_seq_one, sub_seq_none, n)
                sq_start = locs_som_s - 1
                q_slice = som_qual_chars[sq_start:som_end]
                sub_som_qual          = [q_slice[i] for i in pm  if i < len(q_slice)]
                sub_som_qual_hi_chars = [q_slice[i] for i in phi if i < len(q_slice)]

            elif jump == 2:
                n = min(len(sub_seq_none), len(sub_seq_two))
                pm, phi = _fwd_comparison(sub_seq_two, sub_seq_none, n)
                sq_start = locs_som_s - 1 + repeat_length
                q_slice = som_qual_chars[sq_start:som_end]
                sub_som_qual          = [q_slice[i] for i in pm  if i < len(q_slice)]
                sub_som_qual_hi_chars = [q_slice[i] for i in phi if i < len(q_slice)]
            else:
                sub_som_qual          = [None]
                sub_som_qual_hi_chars = [None]

    # ---- 4. Compute quality metrics ------------------------------------------
    clean_sq   = [c for c in (sub_som_qual          or []) if c is not None]
    clean_sq_h = [c for c in (sub_som_qual_hi_chars or []) if c is not None]

    sqv  = _qual_values(clean_sq)
    sqhv = _qual_values(clean_sq_h)

    return dict(
        rateSTRT_mis = START_MISMATCH_RATE,
        nSTRT_mis    = START_MISMATCH_BP_CHECKED,
        flnk_l_hiP     = L_FLANK_GOOD,
        flnk_l_hiP_sub = L_FLANK_SUB_GOOD,
        flnk_r_hiP     = R_FLANK_GOOD,
        flnk_r_hiP_sub = R_FLANK_SUB_GOOD,
        LEN        = len(clean_sq),
        Q          = "".join(clean_sq),
        NONF_SCORE = _safe_mean([v < 30 for v in sqv]),
        HI_QUAL_Q  = "".join(clean_sq_h),
        HI_QUAL_LEN  = len(clean_sq_h),
        NONF_HIQUAL  = _safe_mean([v < 30 for v in sqhv]),
    )


# ---------------------------------------------------------------------------
# Main pipeline
# mirrors the dplyr chain in the R script step-by-step
# ---------------------------------------------------------------------------

# ---- Step 1: group_by(V1, V3) -----------------------------------------------
# For each (V1, V3) group compute:
#   nREADs, midSEG (most common V4), midSEGlen, SEQ (consensus), QUAL

def most_common_str(series):
    return series.mode().iloc[0] if len(series) > 0 else None

grp1 = tab.groupby(["V1", "V3"])

midSEG_map = grp1["V4"].agg(most_common_str).rename("midSEG")
nREADs_map = grp1["V4"].transform("count")

tab_sum_temp = tab.copy()
tab_sum_temp["nREADs"] = nREADs_map.values

midSEG_full = grp1["V4"].transform(most_common_str)
tab_sum_temp["midSEG"]    = midSEG_full.values
tab_sum_temp["midSEGlen"] = tab_sum_temp["midSEG"].apply(len)

# SEQ: consensus per (V1,V3) group; QUAL: V6 if single read else "Consensus"
def _group_consensus(grp):
    seqs  = grp["V5"].astype(str).tolist()
    quals = grp["V6"].astype(str).tolist()
    n     = len(seqs)
    cons  = consensus(seqs, quals, STARTstr)
    qual_val = grp["V6"].iloc[0] if n == 1 else "Consensus"
    return pd.DataFrame({"SEQ": [cons] * n, "QUAL": [qual_val] * n},
                        index=grp.index)

seq_qual = tab.groupby(["V1", "V3"], group_keys=False).apply(_group_consensus)
tab_sum_temp["SEQ"]  = seq_qual["SEQ"].values
tab_sum_temp["QUAL"] = seq_qual["QUAL"].values

# ---- Step 2: group_by(V1) → nALLELES; filter nALLELES >= 2 ----------------
tab_sum_temp["nALLELES"] = (
    tab_sum_temp.groupby("V1")["midSEGlen"]
    .transform(lambda x: x.nunique())
)
tab_sum_temp = tab_sum_temp[tab_sum_temp["nALLELES"] >= 2].copy()

print(f"DEBUG: After nALLELES >= 2 filter: {len(tab_sum_temp)} rows", file=sys.stderr)
if len(tab_sum_temp) > 0:
    print(f"DEBUG: Unique allele lengths: {sorted(tab_sum_temp['V3'].unique())}", file=sys.stderr)
    print(tab_sum_temp.groupby("V3").size().reset_index(name="count").to_string(), file=sys.stderr)

# ---- Step 3: Identify A1 and A2 alleles per (V1, nALLELES) -----------------
# A1: allele with most reads; A2: allele with most reads that differs from A1

# Build allele lookup with an explicit loop to avoid pandas-version-dependent
# groupby.apply behaviour when the function returns a pd.Series.
_allele_rows = []
for (v1, nalleles), _grp in tab_sum_temp.groupby(["V1", "nALLELES"]):
    _ordered = _grp.sort_values("nREADs", ascending=False)
    _a1len   = _ordered["midSEGlen"].iloc[0]
    _a2_rows = _ordered[_ordered["midSEGlen"] != _a1len]
    if len(_a2_rows) == 0:
        continue
    _allele_rows.append({
        "V1": v1, "nALLELES": nalleles,
        "A1len": _a1len,
        "A2len": _a2_rows["midSEGlen"].iloc[0],
        "A1consensus": _ordered["SEQ"].iloc[0],
        "A2consensus": _a2_rows["SEQ"].iloc[0],
        "A1midconsensus": _ordered["midSEG"].iloc[0],
        "A2midconsensus": _a2_rows["midSEG"].iloc[0],
    })
allele_df = pd.DataFrame(_allele_rows)

tab_sum_temp = tab_sum_temp.merge(allele_df, on=["V1", "nALLELES"], how="left")

# ---- Step 4: Apply read-count and allele-difference filters -----------------
# R filter conditions (per V1 group):
#   1) second allele supported by >=3 reads
#   2) both alleles supported by >2 reads  (sum(...>2)==2)
#   3) abs(A1len - A2len) >= 5*repLen

def _passes_filter(grp):
    a1len = grp["A1len"].iloc[0]
    a2len = grp["A2len"].iloc[0]
    # per-(V1,V3) read counts
    by_allele = grp.groupby("midSEGlen")["nREADs"].first()
    if a2len not in by_allele.index:
        return False
    n2 = by_allele[a2len]
    both_ok = (by_allele > 2).sum() == 2
    diff_ok = abs(a1len - a2len) >= 5 * repLen
    return (n2 >= 3) and bool(both_ok) and diff_ok

keep_mask = (
    tab_sum_temp.groupby("V1", group_keys=False)
    .apply(lambda g: pd.Series(_passes_filter(g), index=g.index))
)
tab_sum_temp = tab_sum_temp[keep_mask].copy()

print(f"DEBUG: After all filtering tab_sum_temp has {len(tab_sum_temp)} rows", file=sys.stderr)
if len(tab_sum_temp) > 0:
    print(f"DEBUG: A1len: {tab_sum_temp['A1len'].unique()} bp", file=sys.stderr)
    print(f"DEBUG: A2len: {tab_sum_temp['A2len'].unique()} bp", file=sys.stderr)
    a1 = tab_sum_temp["A1len"].iloc[0]
    a2 = tab_sum_temp["A2len"].iloc[0]
    print(f"DEBUG: Allele difference: {abs(a1 - a2)} bp", file=sys.stderr)
    print(f"DEBUG: Required difference: {5 * repLen} bp", file=sys.stderr)

# ---- Step 5: Build tab_sum if non-empty -------------------------------------
if len(tab_sum_temp) > 0:

    # SOM_POT: read within ±2 repeat units of either allele
    tab_sum = tab_sum_temp.copy()
    tab_sum["SOM_POT"] = (
        (abs(tab_sum["midSEGlen"] - tab_sum["A1len"]) <= 2 * repLen) |
        (abs(tab_sum["midSEGlen"] - tab_sum["A2len"]) <= 2 * repLen)
    ).astype(int)
    tab_sum = tab_sum[tab_sum["SOM_POT"] == 1].copy()

    # fromLEN / fromSEQ / jump  (rowwise)
    def _from_len(row):
        msl = row["midSEGlen"]
        if msl == row["A1len"] or msl == row["A2len"]:
            return float("nan")
        diffs = [abs(msl - row["A1len"]), abs(msl - row["A2len"])]
        return row["A1len"] if diffs[0] <= diffs[1] else row["A2len"]

    def _from_seq(row):
        msl = row["midSEGlen"]
        if msl == row["A1len"]:
            return str(row["A1consensus"])
        if msl == row["A2len"]:
            return str(row["A2consensus"])
        diffs = [abs(msl - row["A1len"]), abs(msl - row["A2len"])]
        return str(row["A1consensus"]) if diffs[0] <= diffs[1] else str(row["A2consensus"])

    tab_sum["fromLEN"] = tab_sum.apply(_from_len, axis=1)
    tab_sum["fromSEQ"] = tab_sum.apply(_from_seq, axis=1)

    tab_sum["jump"] = tab_sum.apply(
        lambda r: float("nan") if math.isnan(r["fromLEN"])
        else (r["midSEGlen"] - r["fromLEN"]) / repLen,
        axis=1,
    )

    # Keep only jumps in {-2,-1,1,2,NA}
    valid_jumps = {-2, -1, 1, 2}
    tab_sum = tab_sum[
        tab_sum["jump"].apply(lambda j: math.isnan(j) or j in valid_jumps)
    ].copy()

    # ---- Call mismatch_bases for each of the 4 jump values (rowwise) --------
    def _mb(row, jmp):
        return mismatch_bases(
            strand        = str(row["V2"]),
            jump          = jmp,
            som_seq       = str(row["V5"]),
            som_qual      = str(row["V6"]),
            seq           = str(row["fromSEQ"]),
            start         = STARTstr,
            end           = ENDstr,
            L_flank_set   = L_FLANK,
            R_flank_set   = R_FLANK,
            repeat_length = repLen,
        )

    for col, jmp in [("j_neg2", -2), ("j_neg1", -1),
                     ("j_pos1",  1), ("j_pos2",  2)]:
        _jmp = jmp
        tab_sum[col] = pd.Series(
            [_mb(row, _jmp) for _, row in tab_sum.iterrows()],
            dtype=object, index=tab_sum.index
        )

    # ---- somatic_* flags  ---------------------------------------------------
    def _somatic_flag(row, jcol):
        j = row[jcol]
        l_ok = j["flnk_l_hiP_sub"] == 1 or math.isnan(float(j["flnk_l_hiP_sub"]
                                                              if j["flnk_l_hiP_sub"] is not None
                                                              else float("nan")))
        r_ok = j["flnk_r_hiP_sub"] == 1 or math.isnan(float(j["flnk_r_hiP_sub"]
                                                              if j["flnk_r_hiP_sub"] is not None
                                                              else float("nan")))
        try:
            len_ok  = j["LEN"] >= 4
            nonf_ok = j["NONF_SCORE"] < 0.2
            smr     = j["rateSTRT_mis"]
            smr_ok  = (not math.isnan(smr)) and smr < 0.05
            return int(len_ok and nonf_ok and smr_ok and l_ok and r_ok)
        except (TypeError, ValueError):
            return 0

    for col, jcol in [("somatic_neg2", "j_neg2"), ("somatic_neg1", "j_neg1"),
                      ("somatic_pos1", "j_pos1"), ("somatic_pos2", "j_pos2")]:
        _jcol = jcol
        tab_sum[col] = pd.Series(
            [_somatic_flag(row, _jcol) for _, row in tab_sum.iterrows()],
            dtype=object, index=tab_sum.index
        )

    # ---- Select output columns (mirrors R select()) -------------------------
    print(f"DEBUG tab_sum columns: {tab_sum.columns.tolist()}", file=sys.stderr)
    out_cols = [
        "V1", "V3", "V4", "jump", "fromLEN",
        "A1len", "A1midconsensus", "A2len", "A2midconsensus",
        "somatic_neg2", "somatic_neg1", "somatic_pos1", "somatic_pos2",
    ]
    tab_out = tab_sum[out_cols].copy()

    # ---- Write output (mirrors write.table with row.names=F, quote=F) -------
    tab_out.to_csv(
        f"summary_somatic_{ID}_py.txt",
        sep=" ",
        index=False,
        quoting=csv.QUOTE_NONE,
        escapechar="\\",
        na_rep="NA",
    )
