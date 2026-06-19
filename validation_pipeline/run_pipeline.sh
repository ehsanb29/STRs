#!/usr/bin/env bash
# run_pipeline.sh
#
# Wrapper that runs the three-step somatic STR validation pipeline in order:
#   Step 1 : extract_loci_reads_from_EH_bams.sh
#   Step 2 : compare_bam_vs_txt.sh
#   Step 3 : find_somatic_txt_reads_in_cram.sh
#            (only runs if step 2 produced any somatic-only sequences)
#
# Usage:
#   bash run_pipeline.sh \
#        BAM_LIST  LOCI_FILE  TXT_FILE_LIST  CRAM_LIST  REF_FASTA  [OUT_DIR]
#
# Arguments:
#   BAM_LIST      : plain-text file, one ExpansionHunter BAM path per line
#   LOCI_FILE     : TSV with header: loci_id  chr  start  end
#   TXT_FILE_LIST : plain-text file, one somatic-caller TXT path per line
#   CRAM_LIST     : plain-text file, one original CRAM path (or URL) per line
#   REF_FASTA     : GRCh38 reference FASTA (needed for CRAM decoding)
#   OUT_DIR       : root output directory (default: current directory)
#                   Sub-directories step1/, step2/, step3/ are created inside it.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── argument handling ─────────────────────────────────────────────────────────
if [[ $# -lt 5 ]]; then
    echo "Usage: bash $(basename "$0") BAM_LIST LOCI_FILE TXT_FILE_LIST CRAM_LIST REF_FASTA [OUT_DIR]" >&2
    exit 1
fi

BAM_LIST="$1"
LOCI_FILE="$2"
TXT_FILE_LIST="$3"
CRAM_LIST="$4"
REF_FASTA="$5"
ROOT_OUT="${6:-.}"

# ── sanity checks ─────────────────────────────────────────────────────────────
for f in "$BAM_LIST" "$LOCI_FILE" "$TXT_FILE_LIST" "$CRAM_LIST" "$REF_FASTA"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: file not found: $f" >&2
        exit 1
    fi
done

mkdir -p "${ROOT_OUT}/step1" "${ROOT_OUT}/step2" "${ROOT_OUT}/step3"

# ── Step 1 : extract reads from EH BAMs ───────────────────────────────────────
echo "============================================================"
echo " STEP 1 : extract_loci_reads_from_EH_bams.sh"
echo "============================================================"
bash "${SCRIPT_DIR}/extract_loci_reads_from_EH_bams.sh" \
    "$BAM_LIST" \
    "$LOCI_FILE" \
    "${ROOT_OUT}/step1"

MANIFEST_STEP1="${ROOT_OUT}/step1/manifest.tsv"
if [[ ! -f "$MANIFEST_STEP1" ]]; then
    echo "ERROR: step 1 did not produce manifest.tsv" >&2
    exit 1
fi
echo "Step 1 complete. Manifest: ${MANIFEST_STEP1}"
echo ""

# ── Step 2 : compare EH BAM reads vs somatic TXT ──────────────────────────────
echo "============================================================"
echo " STEP 2 : compare_bam_vs_txt.sh"
echo "============================================================"
bash "${SCRIPT_DIR}/compare_bam_vs_txt.sh" \
    "$MANIFEST_STEP1" \
    "$TXT_FILE_LIST" \
    "${ROOT_OUT}/step2"

MANIFEST_STEP2="${ROOT_OUT}/step2/bam_vs_txt_manifest.tsv"
if [[ ! -f "$MANIFEST_STEP2" ]]; then
    echo "ERROR: step 2 did not produce bam_vs_txt_manifest.tsv" >&2
    exit 1
fi
echo "Step 2 complete. Manifest: ${MANIFEST_STEP2}"
echo ""

# ── Check whether any somatic-only sequences exist ────────────────────────────
# Column 6 of bam_vs_txt_manifest.tsv is the path to the only_somatic_txt file.
# Step 3 is skipped if every such file is empty (0 sequences).
HAS_SOMATIC_ONLY=0
while IFS=$'\t' read -r SAMPLE LOCI _COMMON _ONLY_BAM _TXT_COMMON ONLY_TXT; do
    [[ "$SAMPLE" == "sample_id" || -z "$SAMPLE" ]] && continue
    if [[ -f "$ONLY_TXT" && -s "$ONLY_TXT" ]]; then
        HAS_SOMATIC_ONLY=1
        break
    fi
done < "$MANIFEST_STEP2"

if [[ $HAS_SOMATIC_ONLY -eq 0 ]]; then
    echo "No somatic-only sequences found in any sample/locus pair."
    echo "Step 3 is not needed — pipeline complete."
    echo ""
    echo "Outputs:"
    echo "  Step 1 : ${ROOT_OUT}/step1/"
    echo "  Step 2 : ${ROOT_OUT}/step2/"
    exit 0
fi

# ── Step 3 : retrieve somatic-only reads from original CRAMs ──────────────────
echo "============================================================"
echo " STEP 3 : find_somatic_txt_reads_in_cram.sh"
echo "============================================================"
bash "${SCRIPT_DIR}/find_somatic_txt_reads_in_cram.sh" \
    "$MANIFEST_STEP2" \
    "$CRAM_LIST" \
    "$LOCI_FILE" \
    "$REF_FASTA" \
    "${ROOT_OUT}/step3"

echo "Step 3 complete."
echo ""
echo "Pipeline finished successfully."
echo ""
echo "Outputs:"
echo "  Step 1 : ${ROOT_OUT}/step1/"
echo "  Step 2 : ${ROOT_OUT}/step2/"
echo "  Step 3 : ${ROOT_OUT}/step3/"
