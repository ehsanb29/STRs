#!/usr/bin/env bash
# extract_loci_reads_from_EH_bams.sh
#
# For each ExpansionHunter output BAM and each locus, extract reads that
# overlap the locus region (±EXTEND bp) AND whose read name contains the
# locus ID.  Writes one indexed BAM per (sample, locus) pair.
#
# Usage:
#   bash extract_loci_reads_from_EH_bams.sh BAM_LIST LOCI_FILE [OUT_DIR]
#
# Arguments:
#   BAM_LIST  : plain-text file, one BAM path per line (can include blank lines / # comments)
#   LOCI_FILE : tab-separated file with a header line and columns:
#                 loci_id  chr  start  end
#               e.g.:
#                 TCF4   chr18  55586154  55586228
#                 HTT    chr4   3074876   3074933
#   OUT_DIR   : output directory (default: current directory)
#
# Output:
#   OUT_DIR/<sample_id>.<loci_id>.bam   (and .bam.bai index)
#
# Dependencies: samtools (must be on PATH)

set -euo pipefail

# ── argument handling ────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: bash $(basename "$0") BAM_LIST LOCI_FILE [OUT_DIR]" >&2
    exit 1
fi

BAM_LIST="$1"
LOCI_FILE="$2"
OUT_DIR="${3:-.}"
EXTEND=4000

# ── sanity checks ────────────────────────────────────────────────────────────
for f in "$BAM_LIST" "$LOCI_FILE"; do
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

# ── manifest file (created/overwritten at the start of each run) ────────────
MANIFEST="${OUT_DIR}/manifest.tsv"
echo -e "sample_id\tloci_id\tn_reads\tbam_path" > "$MANIFEST"

# ── main loop ────────────────────────────────────────────────────────────────
# Read loci file (skip header / comment lines)
while IFS=$'\t' read -r LOCI_ID CHR START END REST; do

    # Skip blank lines, comment lines, and the header row
    [[ -z "$LOCI_ID" || "$LOCI_ID" =~ ^[[:space:]]*# || "$LOCI_ID" == "loci_id" ]] && continue

    # Extended region (clamp start to ≥1)
    EXT_START=$(( START - EXTEND ))
    (( EXT_START < 1 )) && EXT_START=1
    EXT_END=$(( END + EXTEND ))
    REGION="${CHR}:${EXT_START}-${EXT_END}"

    echo "=== Locus: ${LOCI_ID}  region: ${REGION} ==="

    while IFS= read -r BAM_FILE; do
        # Skip blank lines and comment lines in the BAM list
        [[ -z "$BAM_FILE" || "$BAM_FILE" =~ ^[[:space:]]*# ]] && continue

        if [[ ! -f "$BAM_FILE" ]]; then
            echo "  WARNING: BAM not found, skipping: $BAM_FILE" >&2
            continue
        fi

        # Ensure the BAM is coordinate-sorted and indexed
        IS_SORTED=$(samtools view -H "$BAM_FILE" | grep -c "SO:coordinate" || true)
        HAS_INDEX=0
        [[ -f "${BAM_FILE}.bai" || -f "${BAM_FILE%.bam}.bai" ]] && HAS_INDEX=1

        if [[ "$IS_SORTED" -eq 0 ]]; then
            echo "  BAM is not coordinate-sorted; sorting: $(basename "$BAM_FILE")"
            samtools sort -@ 16 -o "${BAM_FILE%.bam}_sorted.bam" "$BAM_FILE"
            # mv "${BAM_FILE%.bam}_sorted.bam" "$BAM_FILE"
            # Update BAM_FILE to point to the sorted BAM for downstream processing
            BAM_FILE="${BAM_FILE%.bam}_sorted.bam"
            HAS_INDEX=0   # any existing index is now stale
        fi
        if [[ "$HAS_INDEX" -eq 0 ]]; then
            echo "  BAM is not indexed; indexing: $(basename "$BAM_FILE")"
            samtools index "$BAM_FILE"
        fi

        BAM_BASENAME=$(basename "$BAM_FILE" .bam)
        # Derive sample ID by stripping path and .bam extension
        SAMPLE_ID="$BAM_BASENAME"
        # if the samplee ID contains words like realigned or sorted, remove those too (NDAR_INVAA021FG8_realigned_sorted)
        SAMPLE_ID=$(echo "$SAMPLE_ID" | sed -E 's/(_realigned)?(_sorted)?$//')

        OUT_BAM="${OUT_DIR}/${BAM_BASENAME}.${LOCI_ID}.bam"
        TMPFILE="${OUT_BAM%.bam}.tmp.sam"

        echo "  Processing: ${SAMPLE_ID} -> $(basename "$OUT_BAM")"

        # 1) Write the BAM header
        # 2) Append reads overlapping the extended region whose name contains LOCI_ID
        # 3) Convert SAM stream to a sorted, indexed BAM
        {
            samtools view -H "$BAM_FILE"
            samtools view "$BAM_FILE" "$REGION" | grep -F "$LOCI_ID" || true
        } | samtools view -b -o "$OUT_BAM"

        # Index the output BAM
        samtools index "$OUT_BAM"

        N_READS=$(samtools view -c "$OUT_BAM")
        echo "    Reads written: ${N_READS}"

        # Append entry to manifest
        echo -e "${SAMPLE_ID}\t${LOCI_ID}\t${N_READS}\t$(realpath "$OUT_BAM")" >> "$MANIFEST"

    done < "$BAM_LIST"

done < "$LOCI_FILE"

echo ""
echo "Done. Output BAMs written to: ${OUT_DIR}"
echo "Manifest written to:          ${MANIFEST}"
