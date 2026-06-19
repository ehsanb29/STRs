#!/usr/bin/env bash
# compare_bam_sets.sh
#
# Given two manifest files produced by extract_loci_reads_from_EH_bams.sh,
# for every (sample_id, loci_id) pair present in BOTH manifests:
#   - find reads common to both BAMs
#   - find reads only in set 1
#   - find reads only in set 2
# and write the results as three indexed output BAMs plus a summary report.
#
# Matching is based on read name (QNAME, column 1 of the SAM record).
#
# Usage:
#   bash compare_bam_sets.sh MANIFEST1 MANIFEST2 [OUT_DIR]
#
# Arguments:
#   MANIFEST1  : manifest.tsv from first run  (sample_id  loci_id  bam_path)
#   MANIFEST2  : manifest.tsv from second run
#   OUT_DIR    : output directory (default: current directory)
#
# Output per matched (sample, locus) pair:
#   <sample>.<locus>.common.bam      reads present in both BAMs
#   <sample>.<locus>.only_set1.bam   reads present only in set 1
#   <sample>.<locus>.only_set2.bam   reads present only in set 2
#   (each BAM is indexed)
#
# Summary files:
#   comparison_report.tsv   counts per pair + paths to the three output BAMs
#   comparison_manifest.tsv sample_id / loci_id / paths to the three BAMs
#
# Dependencies: samtools ≥1.10 (for -N flag), standard GNU coreutils

set -euo pipefail

# ── argument handling ────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: bash $(basename "$0") MANIFEST1 MANIFEST2 [OUT_DIR]" >&2
    exit 1
fi

MANIFEST1="$1"
MANIFEST2="$2"
OUT_DIR="${3:-.}"

# ── sanity checks ────────────────────────────────────────────────────────────
for f in "$MANIFEST1" "$MANIFEST2"; do
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

# Temp directory cleaned up on exit
WORK_DIR=$(mktemp -d)
trap 'rm -rf "$WORK_DIR"' EXIT

# ── output summary files ─────────────────────────────────────────────────────
REPORT="${OUT_DIR}/comparison_report.tsv"
OUT_MANIFEST="${OUT_DIR}/comparison_manifest.tsv"

printf '%s\n' \
    "sample_id"$'\t'"loci_id"$'\t'"reads_set1"$'\t'"reads_set2"$'\t'"reads_common"$'\t'"reads_only_set1"$'\t'"reads_only_set2"$'\t'"bam_common"$'\t'"bam_only_set1"$'\t'"bam_only_set2" \
    > "$REPORT"

printf '%s\n' \
    "sample_id"$'\t'"loci_id"$'\t'"bam_common"$'\t'"bam_only_set1"$'\t'"bam_only_set2" \
    > "$OUT_MANIFEST"

# ── helper: write a BAM containing only reads whose name is in a list file ───
# Usage: filter_bam  NAMES_FILE  SOURCE_BAM  OUT_BAM
# If the names file is empty, the output BAM contains only the header.
filter_bam() {
    local names_file="$1" src_bam="$2" out_bam="$3"
    if [[ ! -s "$names_file" ]]; then
        # No matching reads – write header-only BAM
        samtools view -H "$src_bam" | samtools view -b -o "$out_bam"
    else
        samtools view -h -N "$names_file" "$src_bam" | samtools view -b -o "$out_bam"
    fi
    samtools index "$out_bam"
}

# ── load manifest 1 into an associative array  key="sample__loci" → bam_path ─
declare -A MAP1
while IFS=$'\t' read -r SAMPLE LOCI BAM_PATH; do
    [[ "$SAMPLE" == "sample_id" || -z "$SAMPLE" ]] && continue
    MAP1["${SAMPLE}__${LOCI}"]="$BAM_PATH"
done < "$MANIFEST1"

# ── load manifest 2 ──────────────────────────────────────────────────────────
declare -A MAP2
while IFS=$'\t' read -r SAMPLE LOCI BAM_PATH; do
    [[ "$SAMPLE" == "sample_id" || -z "$SAMPLE" ]] && continue
    MAP2["${SAMPLE}__${LOCI}"]="$BAM_PATH"
done < "$MANIFEST2"

# ── iterate over keys present in both manifests ──────────────────────────────
FOUND=0

for KEY in "${!MAP1[@]}"; do

    [[ -v MAP2["$KEY"] ]] || continue

    FOUND=$(( FOUND + 1 ))
    SAMPLE="${KEY%%__*}"
    LOCI="${KEY##*__}"
    BAM1="${MAP1[$KEY]}"
    BAM2="${MAP2[$KEY]}"

    echo "=== ${SAMPLE} / ${LOCI} ==="
    echo "  Set1 BAM: ${BAM1}"
    echo "  Set2 BAM: ${BAM2}"

    # Verify both BAMs exist before proceeding
    SKIP=0
    for b in "$BAM1" "$BAM2"; do
        if [[ ! -f "$b" ]]; then
            echo "  WARNING: BAM not found, skipping pair: $b" >&2
            SKIP=1
        fi
    done
    [[ $SKIP -eq 1 ]] && continue

    # ── extract sorted unique read names from each BAM ───────────────────────
    NAMES1="${WORK_DIR}/${SAMPLE}.${LOCI}.set1.names"
    NAMES2="${WORK_DIR}/${SAMPLE}.${LOCI}.set2.names"

    samtools view "$BAM1" | awk '{print $1}' | sort -u > "$NAMES1"
    samtools view "$BAM2" | awk '{print $1}' | sort -u > "$NAMES2"

    N1=$(wc -l < "$NAMES1")
    N2=$(wc -l < "$NAMES2")

    # ── set operations using comm (both inputs are already sorted) ────────────
    NAMES_COMMON="${WORK_DIR}/${SAMPLE}.${LOCI}.common.names"
    NAMES_ONLY1="${WORK_DIR}/${SAMPLE}.${LOCI}.only_set1.names"
    NAMES_ONLY2="${WORK_DIR}/${SAMPLE}.${LOCI}.only_set2.names"

    comm -12 "$NAMES1" "$NAMES2" > "$NAMES_COMMON"   # in both
    comm -23 "$NAMES1" "$NAMES2" > "$NAMES_ONLY1"    # only in set1
    comm -13 "$NAMES1" "$NAMES2" > "$NAMES_ONLY2"    # only in set2

    N_COMMON=$(wc -l < "$NAMES_COMMON")
    N_ONLY1=$(wc -l < "$NAMES_ONLY1")
    N_ONLY2=$(wc -l < "$NAMES_ONLY2")

    echo "  Reads: set1=${N1}  set2=${N2}  common=${N_COMMON}  only_set1=${N_ONLY1}  only_set2=${N_ONLY2}"

    # ── write output BAMs ─────────────────────────────────────────────────────
    OUT_COMMON="${OUT_DIR}/${SAMPLE}.${LOCI}.common.bam"
    OUT_ONLY1="${OUT_DIR}/${SAMPLE}.${LOCI}.only_set1.bam"
    OUT_ONLY2="${OUT_DIR}/${SAMPLE}.${LOCI}.only_set2.bam"

    echo "  Writing common BAM   -> $(basename "$OUT_COMMON")"
    filter_bam "$NAMES_COMMON" "$BAM1" "$OUT_COMMON"

    echo "  Writing only_set1 BAM -> $(basename "$OUT_ONLY1")"
    filter_bam "$NAMES_ONLY1" "$BAM1" "$OUT_ONLY1"

    echo "  Writing only_set2 BAM -> $(basename "$OUT_ONLY2")"
    filter_bam "$NAMES_ONLY2" "$BAM2" "$OUT_ONLY2"

    # ── append to report and manifest ─────────────────────────────────────────
    printf '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n' \
        "$SAMPLE" "$LOCI" "$N1" "$N2" "$N_COMMON" "$N_ONLY1" "$N_ONLY2" \
        "$(realpath "$OUT_COMMON")" \
        "$(realpath "$OUT_ONLY1")" \
        "$(realpath "$OUT_ONLY2")" \
        >> "$REPORT"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$SAMPLE" "$LOCI" \
        "$(realpath "$OUT_COMMON")" \
        "$(realpath "$OUT_ONLY1")" \
        "$(realpath "$OUT_ONLY2")" \
        >> "$OUT_MANIFEST"

done

# ── final summary ─────────────────────────────────────────────────────────────
if [[ $FOUND -eq 0 ]]; then
    echo "WARNING: no matching (sample_id, loci_id) pairs found between the two manifests." >&2
else
    echo ""
    echo "Processed ${FOUND} matching (sample, locus) pair(s)."
fi

echo ""
echo "Comparison report : ${REPORT}"
echo "Output manifest   : ${OUT_MANIFEST}"
