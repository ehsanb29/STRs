INPUT_DIR=$1
OUTPUT_PREFIX=$2

extract_repcn() {
    vcf="$1"
    sample=$(bcftools query -l "$vcf")
    bcftools query -f "%CHROM\t%POS\t%INFO/VARID\t${sample}\t[%REPCN]\n" "$vcf"
}
export -f extract_repcn

find "$INPUT_DIR" -name "*.vcf" | \
  parallel -j $(nproc) extract_repcn {} >> "${OUTPUT_PREFIX}_all_repcn.tsv"