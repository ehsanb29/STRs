#!/usr/bin/env bash
# compare_bam_vs_txt.sh
#
# For each (sample_id, loci_id) row in a manifest file, compare:
#   - read sequences in the BAM file          (SAM column 10)
#   - read sequences in one or more TXT files (column 5, rows where column 1 == sample_id)
# and produce per-category output files plus a summary report.
#
# TXT file format (space-separated, no mandatory header):
#   col 1 : sample_id
#   col 2 : strand (forward/reverse)
#   col 3 : repeat length
#   col 4 : repeat+flank sequence
#   col 5 : full read sequence  ← compared against BAM column 10
#   col 6 : base-quality string
#
# Usage:
#   bash compare_bam_vs_txt.sh MANIFEST TXT_FILE_LIST [OUT_DIR]
#
# Arguments:
#   MANIFEST      : manifest.tsv from extract_loci_reads_from_EH_bams.sh
#                   (columns: sample_id  loci_id  bam_path)
#   TXT_FILE_LIST : plain-text file with one TXT file path per line
#   OUT_DIR       : output directory (default: current directory)
#
# Output per matched (sample_id, loci_id) pair:
#   <sample>.<locus>.common.bam     BAM reads whose sequence is also in the TXT file(s)
#   <sample>.<locus>.common.txt     TXT rows whose sequence is also in the BAM
#   <sample>.<locus>.only_bam.bam   BAM reads whose sequence is NOT in the TXT file(s)
#   <sample>.<locus>.only_txt.txt   TXT rows whose sequence is NOT in the BAM
#   (BAM outputs are indexed)
#
# Summary files written to OUT_DIR:
#   bam_vs_txt_report.tsv    read counts + paths for every matched pair
#   bam_vs_txt_manifest.tsv  paths to the four output files per pair
#
# Dependencies: samtools ≥1.10 (for -N flag), GNU awk, GNU coreutils

set -euo pipefail

# ── argument handling ─────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: bash $(basename "$0") MANIFEST TXT_FILE_LIST [OUT_DIR]" >&2
    exit 1
fi

MANIFEST="$1"
TXT_FILE_LIST="$2"
OUT_DIR="${3:-.}"

# ── sanity checks ─────────────────────────────────────────────────────────────
for f in "$MANIFEST" "$TXT_FILE_LIST"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: file not found: $f" >&2
        exit 1
    fi
done

if ! command -v samtools &>/dev/null; then
    echo "ERROR: samtools not found on PATH" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

# ── output summary files ──────────────────────────────────────────────────────
REPORT="${OUT_DIR}/bam_vs_txt_report.tsv"
OUT_MANIFEST="${OUT_DIR}/bam_vs_txt_manifest.tsv"

printf 'sample_id\tloci_id\tseqs_eh_bam\tseqs_somatic_txt\tseqs_common\tseqs_only_eh_bam\tseqs_only_somatic_txt\teh_bam_common\teh_bam_only_bam\tsomatic_txt_common\tsomatic_txt_only_txt\n' \
    > "$REPORT"

printf 'sample_id\tloci_id\teh_bam_common\teh_bam_only_bam\tsomatic_txt_common\tsomatic_txt_only_txt\n' \
    > "$OUT_MANIFEST"

# ── helper: write a BAM keeping only reads whose SEQ (col 10) is in SEQ_FILE ─
# Usage: filter_bam_by_seq  SEQ_FILE  SOURCE_BAM  OUT_BAM
# If SEQ_FILE is empty the output BAM contains only the header (0 reads).
filter_bam_by_seq() {
    local seq_file="$1" src_bam="$2" out_bam="$3"
    local names_tmp="${WORK_DIR}/$(basename "$out_bam").names"

    if [[ ! -s "$seq_file" ]]; then
        # No target sequences → header-only BAM
        samtools view -H "$src_bam" | samtools view -b -o "$out_bam"
    else
        # Collect read names (QNAME) whose sequence matches
        samtools view "$src_bam" | awk \
            -v seqfile="$seq_file" \
            'BEGIN { while ((getline s < seqfile) > 0) seqs[s]=1 }
             $10 in seqs { print $1 }' \
            | sort -u > "$names_tmp"

        if [[ ! -s "$names_tmp" ]]; then
            samtools view -H "$src_bam" | samtools view -b -o "$out_bam"
        else
            samtools view -h -N "$names_tmp" "$src_bam" | samtools view -b -o "$out_bam"
        fi
    fi
    samtools index "$out_bam"
}

# ── pre-collect sample IDs present across all TXT files ──────────────────────
# Build a lookup: which txt files contain each sample_id?
# We only care about sample IDs that also appear in the manifest.

# First pass: load sample IDs from manifest
declare -A MANIFEST_SAMPLES   # sample_id -> "seen"
while IFS=$'\t' read -r SAMPLE LOCI N_READS BAM_PATH; do
    [[ "$SAMPLE" == "sample_id" || -z "$SAMPLE" ]] && continue
    MANIFEST_SAMPLES["$SAMPLE"]=1
done < "$MANIFEST"

# ── main loop: iterate over manifest rows ─────────────────────────────────────
PROCESSED=0
SKIPPED=0

while IFS=$'\t' read -r SAMPLE LOCI N_READS BAM_PATH; do
    [[ "$SAMPLE" == "sample_id" || -z "$SAMPLE" ]] && continue

    # Verify BAM exists
    if [[ ! -f "$BAM_PATH" ]]; then
        echo "WARNING: BAM not found, skipping: $BAM_PATH" >&2
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    echo "=== ${SAMPLE} / ${LOCI} ==="

    # ── extract unique sequences from BAM (column 10) ─────────────────────────
    BAM_SEQS="${WORK_DIR}/${SAMPLE}.${LOCI}.bam.seqs"
    samtools view "$BAM_PATH" | awk '{print $10}' | sort -u > "$BAM_SEQS"
    N_BAM=$(wc -l < "$BAM_SEQS")
    echo "  BAM sequences : ${N_BAM}"

    # ── collect sequences from all TXT files for this sample (column 5) ───────
    TXT_SEQS="${WORK_DIR}/${SAMPLE}.${LOCI}.txt.seqs"
    # Also keep full rows (all columns) for txt-output files
    TXT_ROWS="${WORK_DIR}/${SAMPLE}.${LOCI}.txt.rows"
    > "$TXT_SEQS"
    > "$TXT_ROWS"

    while IFS= read -r TXT_FILE; do
        [[ -z "$TXT_FILE" || "$TXT_FILE" =~ ^[[:space:]]*# ]] && continue
        if [[ ! -f "$TXT_FILE" ]]; then
            echo "  WARNING: txt file not found, skipping: $TXT_FILE" >&2
            continue
        fi
        # Only use TXT files whose filename contains this locus ID.
        # The somatic caller names files IID_<sample>_<gene>_<rep>.txt so the
        # locus (gene) name is embedded in the filename.  Without this filter,
        # reads from other loci for the same sample would contaminate the
        # comparison against the locus-specific BAM from script 1.
        TXT_BASENAME=$(basename "$TXT_FILE")
        if [[ "$TXT_BASENAME" != *"${LOCI}"* ]]; then
            continue
        fi
        # Rows matching this sample_id; skip any header where col1 != a real ID
        awk -v sid="$SAMPLE" '$1==sid { print $5 }' "$TXT_FILE" >> "$TXT_SEQS"
        awk -v sid="$SAMPLE" '$1==sid { print $0  }' "$TXT_FILE" >> "$TXT_ROWS"
    done < "$TXT_FILE_LIST"

    sort -u "$TXT_SEQS" -o "$TXT_SEQS"
    N_TXT=$(wc -l < "$TXT_SEQS")
    echo "  TXT sequences : ${N_TXT}"

    if [[ $N_TXT -eq 0 ]]; then
        echo "  No sequences found in any TXT file for ${SAMPLE} — skipping." >&2
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    # ── set operations (both files are already sorted) ────────────────────────
    SEQS_COMMON="${WORK_DIR}/${SAMPLE}.${LOCI}.common.seqs"
    SEQS_ONLY_BAM="${WORK_DIR}/${SAMPLE}.${LOCI}.only_eh_bam.seqs"
    SEQS_ONLY_TXT="${WORK_DIR}/${SAMPLE}.${LOCI}.only_somatic_txt.seqs"

    comm -12 "$BAM_SEQS" "$TXT_SEQS" > "$SEQS_COMMON"   # in both
    comm -23 "$BAM_SEQS" "$TXT_SEQS" > "$SEQS_ONLY_BAM" # only in BAM
    comm -13 "$BAM_SEQS" "$TXT_SEQS" > "$SEQS_ONLY_TXT" # only in TXT

    N_COMMON=$(wc -l < "$SEQS_COMMON")
    N_ONLY_BAM=$(wc -l < "$SEQS_ONLY_BAM")
    N_ONLY_TXT=$(wc -l < "$SEQS_ONLY_TXT")

    echo "  common=${N_COMMON}  only_bam=${N_ONLY_BAM}  only_txt=${N_ONLY_TXT}"

    # ── write BAM outputs ─────────────────────────────────────────────────────
    OUT_COMMON_BAM="${OUT_DIR}/${SAMPLE}.${LOCI}.common.bam"
    OUT_ONLY_BAM_BAM="${OUT_DIR}/${SAMPLE}.${LOCI}.only_eh_bam.bam"

    echo "  Writing common BAM   -> $(basename "$OUT_COMMON_BAM")"
    filter_bam_by_seq "$SEQS_COMMON"   "$BAM_PATH" "$OUT_COMMON_BAM"

    echo "  Writing only_bam BAM -> $(basename "$OUT_ONLY_BAM_BAM")"
    filter_bam_by_seq "$SEQS_ONLY_BAM" "$BAM_PATH" "$OUT_ONLY_BAM_BAM"

    # ── write TXT outputs (full original rows) ────────────────────────────────
    OUT_COMMON_TXT="${OUT_DIR}/${SAMPLE}.${LOCI}.common.txt"
    OUT_ONLY_TXT="${OUT_DIR}/${SAMPLE}.${LOCI}.only_somatic_txt.txt"

    echo "  Writing common TXT   -> $(basename "$OUT_COMMON_TXT")"
    awk -v seqfile="$SEQS_COMMON" \
        'BEGIN { while ((getline s < seqfile) > 0) seqs[s]=1 }
         $5 in seqs' \
        "$TXT_ROWS" > "$OUT_COMMON_TXT"

    echo "  Writing only_txt TXT -> $(basename "$OUT_ONLY_TXT")"
    awk -v seqfile="$SEQS_ONLY_TXT" \
        'BEGIN { while ((getline s < seqfile) > 0) seqs[s]=1 }
         $5 in seqs' \
        "$TXT_ROWS" > "$OUT_ONLY_TXT"

    # ── append to report and manifest ─────────────────────────────────────────
    printf '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE" "$LOCI" \
        "$N_BAM" "$N_TXT" "$N_COMMON" "$N_ONLY_BAM" "$N_ONLY_TXT" \
        "$(realpath "$OUT_COMMON_BAM")" \
        "$(realpath "$OUT_ONLY_BAM_BAM")" \
        "$(realpath "$OUT_COMMON_TXT")" \
        "$(realpath "$OUT_ONLY_TXT")" \
        >> "$REPORT"

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE" "$LOCI" \
        "$(realpath "$OUT_COMMON_BAM")" \
        "$(realpath "$OUT_ONLY_BAM_BAM")" \
        "$(realpath "$OUT_COMMON_TXT")" \
        "$(realpath "$OUT_ONLY_TXT")" \
        >> "$OUT_MANIFEST"

    PROCESSED=$(( PROCESSED + 1 ))

done < "$MANIFEST"

# ── final summary ─────────────────────────────────────────────────────────────
echo ""
echo "Processed : ${PROCESSED} pair(s)"
[[ $SKIPPED -gt 0 ]] && echo "Skipped   : ${SKIPPED} pair(s) (BAM or TXT not found)"
echo ""
echo "Report   : ${REPORT}"
echo "Manifest : ${OUT_MANIFEST}"
