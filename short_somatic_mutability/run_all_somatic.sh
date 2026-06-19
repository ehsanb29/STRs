#!/usr/bin/env bash
# run_all_somatic.sh
# Usage: bash run_all_somatic.sh <cram_list_file> <locus_string> <OUT_DIR>
#
# Arguments:
#   cram_list_file  : Path to a text file with one CRAM file path per line
#   locus_string    : Locus descriptor passed as DAT to somatic.sh
#                     e.g. TCF4_chr18_55586154_55586228_AGGAGGAGC_AGCATGAAA_AGC_L6,9_R,
#
# Example:
#   bash run_all_somatic.sh cram_paths.txt TCF4_chr18_55586154_55586228_AGGAGGAGC_AGCATGAAA_AGC_L6,9_R,

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: bash $0 <cram_list_file> <locus_string> <OUT_DIR>" >&2
    exit 1
fi

CRAM_LIST="$1"
LOCUS="$2"
OUT_DIR="$3"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -f "$CRAM_LIST" ]; then
    echo "Error: CRAM list file not found: $CRAM_LIST" >&2
    exit 1
fi

if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

total=$(grep -c '.' "$CRAM_LIST" || true)
echo "Found $total CRAM file(s) to process."

count=0
while IFS= read -r CRAM_FILE || [ -n "$CRAM_FILE" ]; do
    # Skip empty lines and comment lines
    [[ -z "$CRAM_FILE" || "$CRAM_FILE" == \#* ]] && continue

    count=$(( count + 1 ))
    echo "=== [$count/$total] Processing: $CRAM_FILE ==="

    # Determine the file to pass to somatic.sh (may be replaced by sorted BAM)
    INPUT_FILE="$CRAM_FILE"

    # Check if an index file exists (.bai or .crai alongside the file)
    BASE="${CRAM_FILE%.*}"
    EXT="${CRAM_FILE##*.}"
    SORTED_BAM="${BASE}_sorted.bam"

    INDEX_EXISTS=false
    # check if the sorted BAM or its index already exists to avoid unnecessary sorting and indexing

    if [ -f "${SORTED_BAM}" ] && [ -f "${SORTED_BAM}.bai" ]; then
        INPUT_FILE="$SORTED_BAM"
        INDEX_EXISTS=true
    fi
    if [ -f "${CRAM_FILE}.bai" ] || [ -f "${BASE}.bai" ]  || \
       [ -f "${CRAM_FILE}.crai" ] || [ -f "${BASE}.crai" ]; then
        INDEX_EXISTS=true
    fi

    if [ "$INDEX_EXISTS" = false ]; then
        echo "  No index found for $CRAM_FILE — sorting and indexing..."
        samtools sort -o "$SORTED_BAM" "$CRAM_FILE"
        samtools index "$SORTED_BAM"
        INPUT_FILE="$SORTED_BAM"
        echo "  Using sorted file: $INPUT_FILE"
    fi

    bash "$SCRIPT_DIR/somatic.sh" "$INPUT_FILE" "$LOCUS" "$OUT_DIR"

done < "$CRAM_LIST"

echo "Done. Processed $count CRAM file(s)."
