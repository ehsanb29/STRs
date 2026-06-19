# !/bin/bash
# Usage: ./extract_repeat_size.sh <input_directory> <output_prefix>

INPUT_DIR=$1
OUTPUT_PREFIX=$2

find "$INPUT_DIR" -name "*.vcf" | \
  xargs -P $(nproc) -I{} sh -c '
    sample=$(bcftools query -l "$1")
    bcftools query -f "%CHROM\t%POS\t%INFO/VARID\t${sample}\t[%REPCN]\n" "$1"
' _ {} 2>/dev/null >> "${OUTPUT_PREFIX}_all_repcn.tsv"