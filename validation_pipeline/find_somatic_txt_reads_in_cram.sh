#!/usr/bin/env bash
# find_somatic_txt_reads_in_cram.sh
#
# For each (sample_id, loci_id) row in the bam_vs_txt_manifest.tsv produced by
# compare_bam_vs_txt.sh, take the sequences that were found ONLY in the somatic
# TXT file (column 5 of only_somatic_txt.txt) and retrieve the corresponding
# full reads directly from the original CRAM file.
#
# Matching is done on read sequence (SAM column 10 == TXT column 5).
# The CRAM is queried over the locus region extended by ±EXTEND bp.
#
# Usage:
#   bash find_somatic_txt_reads_in_cram.sh \
#        BAM_VS_TXT_MANIFEST  CRAM_LIST  LOCI_FILE  REF_FASTA  [OUT_DIR]
#
# Arguments:
#   BAM_VS_TXT_MANIFEST : bam_vs_txt_manifest.tsv from compare_bam_vs_txt.sh
#                         columns: sample_id  loci_id  eh_bam_common  eh_bam_only_bam
#                                  somatic_txt_common  somatic_txt_only_txt
#   CRAM_LIST           : plain-text file, one CRAM path (or URL) per line;
#                         sample_id is derived by stripping path and
#                         .final.cram / .cram suffix (same as somatic.sh)
#   LOCI_FILE           : tab-separated file with header: loci_id  chr  start  end
#                         (same file used with extract_loci_reads_from_EH_bams.sh)
#   REF_FASTA           : GRCh38 reference FASTA (needed to decode CRAM)
#   OUT_DIR             : output directory (default: current directory)
#
# Output per matched (sample_id, loci_id) pair:
#   <sample>.<locus>.somatic_txt_in_cram.bam   reads from CRAM whose sequence
#                                               matches a somatic-txt-only sequence
#   (each BAM is indexed)
#
# Summary files:
#   somatic_in_cram_report.tsv    counts + paths
#   somatic_in_cram_manifest.tsv  sample_id / loci_id / output BAM path
#
# Dependencies: samtools (with CRAM support), GNU awk, GNU coreutils

set -euo pipefail

EXTEND=4000

# ── argument handling ─────────────────────────────────────────────────────────
if [[ $# -lt 4 ]]; then
    echo "Usage: bash $(basename "$0") BAM_VS_TXT_MANIFEST CRAM_LIST LOCI_FILE REF_FASTA [OUT_DIR]" >&2
    exit 1
fi

INPUT_MANIFEST="$1"
CRAM_LIST="$2"
LOCI_FILE="$3"
REF_FASTA="$4"
OUT_DIR="${5:-.}"

# ── sanity checks ─────────────────────────────────────────────────────────────
for f in "$INPUT_MANIFEST" "$CRAM_LIST" "$LOCI_FILE" "$REF_FASTA"; do
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
REPORT="${OUT_DIR}/somatic_in_cram_report.tsv"
OUT_MANIFEST="${OUT_DIR}/somatic_in_cram_manifest.tsv"

printf 'sample_id\tloci_id\tseqs_queried\treads_found_in_cram\toutput_bam\n' > "$REPORT"
printf 'sample_id\tloci_id\toutput_bam\n' > "$OUT_MANIFEST"

# ── build CRAM lookup: sample_id → cram_path ─────────────────────────────────
declare -A CRAM_MAP
while IFS= read -r CRAM_FILE; do
    [[ -z "$CRAM_FILE" || "$CRAM_FILE" =~ ^[[:space:]]*# ]] && continue
    # Mirror the ID extraction used in somatic.sh
    SID=$(basename "$CRAM_FILE" | sed 's/\.final\.cram$//' | sed 's/\.cram$//')
    CRAM_MAP["$SID"]="$CRAM_FILE"
done < "$CRAM_LIST"

echo "Loaded ${#CRAM_MAP[@]} CRAM file(s) from list."

# ── build loci lookup: loci_id → chr, start, end ─────────────────────────────
declare -A LOCI_CHR LOCI_START LOCI_END
while IFS=$'\t' read -r LID LCHR LSTART LEND REST; do
    [[ "$LID" == "loci_id" || -z "$LID" || "$LID" =~ ^[[:space:]]*# ]] && continue
    LOCI_CHR["$LID"]="$LCHR"
    LOCI_START["$LID"]="$LSTART"
    LOCI_END["$LID"]="$LEND"
done < "$LOCI_FILE"

echo "Loaded ${#LOCI_CHR[@]} locus/loci from loci file."

# ── main loop: iterate over manifest rows ─────────────────────────────────────
PROCESSED=0
SKIPPED=0

while IFS=$'\t' read -r SAMPLE LOCI EH_COMMON EH_ONLY_BAM TXT_COMMON TXT_ONLY_TXT; do
    [[ "$SAMPLE" == "sample_id" || -z "$SAMPLE" ]] && continue

    echo "=== ${SAMPLE} / ${LOCI} ==="

    # ── resolve CRAM file ─────────────────────────────────────────────────────
    if [[ ! -v CRAM_MAP["$SAMPLE"] ]]; then
        echo "  WARNING: no CRAM found for sample '${SAMPLE}' — skipping." >&2
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi
    CRAM_FILE="${CRAM_MAP[$SAMPLE]}"
    echo "  CRAM : ${CRAM_FILE}"

    # ── resolve locus coordinates ─────────────────────────────────────────────
    if [[ ! -v LOCI_CHR["$LOCI"] ]]; then
        echo "  WARNING: locus '${LOCI}' not found in loci file — skipping." >&2
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi
    CHR="${LOCI_CHR[$LOCI]}"
    RAW_START="${LOCI_START[$LOCI]}"
    RAW_END="${LOCI_END[$LOCI]}"
    EXT_START=$(( RAW_START - EXTEND ))
    (( EXT_START < 1 )) && EXT_START=1
    EXT_END=$(( RAW_END + EXTEND ))
    REGION="${CHR}:${EXT_START}-${EXT_END}"
    echo "  Region (±${EXTEND}bp): ${REGION}"

    # ── load target sequences from only_somatic_txt file ─────────────────────
    if [[ ! -f "$TXT_ONLY_TXT" ]]; then
        echo "  WARNING: only_somatic_txt file not found: ${TXT_ONLY_TXT} — skipping." >&2
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    TARGET_SEQS="${WORK_DIR}/${SAMPLE}.${LOCI}.target.seqs"
    awk '{print $5}' "$TXT_ONLY_TXT" | sort -u > "$TARGET_SEQS"
    N_TARGET=$(wc -l < "$TARGET_SEQS")
    echo "  Target sequences from only_somatic_txt : ${N_TARGET}"

    if [[ $N_TARGET -eq 0 ]]; then
        echo "  No sequences to search for — skipping." >&2
        SKIPPED=$(( SKIPPED + 1 ))
        continue
    fi

    # ── query CRAM: extract reads whose SEQ matches a target sequence ─────────
    OUT_BAM="${OUT_DIR}/${SAMPLE}.${LOCI}.somatic_txt_in_cram.bam"
    NAMES_TMP="${WORK_DIR}/${SAMPLE}.${LOCI}.cram_match.names"

    echo "  Scanning CRAM region for matching sequences..."

    # Step 1: collect read names from CRAM whose column 10 is a target sequence
    samtools view -T "$REF_FASTA" "$CRAM_FILE" "$REGION" \
        | awk -v seqfile="$TARGET_SEQS" \
              'BEGIN { while ((getline s < seqfile) > 0) seqs[s]=1 }
               $10 in seqs { print $1 }' \
        | sort -u > "$NAMES_TMP"

    N_FOUND=$(wc -l < "$NAMES_TMP")
    echo "  Reads matched in CRAM : ${N_FOUND}"

    # Step 2: write output BAM (header + matched reads)
    if [[ ! -s "$NAMES_TMP" ]]; then
        # No matches — write header-only BAM
        samtools view -H -T "$REF_FASTA" "$CRAM_FILE" | samtools view -b -o "$OUT_BAM"
    else
        samtools view -h -T "$REF_FASTA" -N "$NAMES_TMP" "$CRAM_FILE" "$REGION" \
            | samtools view -b -o "$OUT_BAM"
    fi

    samtools index "$OUT_BAM"
    echo "  Output BAM -> $(basename "$OUT_BAM")"

    # ── append to report and manifest ─────────────────────────────────────────
    printf '%s\t%s\t%d\t%d\t%s\n' \
        "$SAMPLE" "$LOCI" "$N_TARGET" "$N_FOUND" \
        "$(realpath "$OUT_BAM")" \
        >> "$REPORT"

    printf '%s\t%s\t%s\n' \
        "$SAMPLE" "$LOCI" \
        "$(realpath "$OUT_BAM")" \
        >> "$OUT_MANIFEST"

    PROCESSED=$(( PROCESSED + 1 ))

done < "$INPUT_MANIFEST"

# ── final summary ─────────────────────────────────────────────────────────────
echo ""
echo "Processed : ${PROCESSED} pair(s)"
[[ $SKIPPED -gt 0 ]] && echo "Skipped   : ${SKIPPED} pair(s)"
echo ""
echo "Report   : ${REPORT}"
echo "Manifest : ${OUT_MANIFEST}"
