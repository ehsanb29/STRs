set -euo pipefail
FILE="$1"; ID=$( basename "$FILE" | sed 's/\.final\.cram$//' | sed 's/\.cram$//' ); DAT=$2;
OUT_DIR=$3
#  output dir exists check na make it
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi
echo "Processing individual: $ID with data: $DAT"

GENE=$(echo $DAT | cut -d'_' -f 1); REP=$(echo $DAT | cut -d'_' -f 7);
echo "Gene: $GENE with repeat type: $REP"
CHR=$(echo $DAT | cut -d'_' -f 2); START=$(echo $DAT | cut -d'_' -f 3); END=$(echo $DAT | cut -d'_' -f 4);
STARTSEQ=$(echo $DAT | cut -d'_' -f 5); ENDSEQ=$(echo $DAT | cut -d'_' -f 6);

echo "Region: $CHR:$START-$END with flanking sequences $STARTSEQ and $ENDSEQ"

echo "samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa -F 0x4 -F 0x100 -F 0x200 -F 0x400 -F 0x800 $FILE ${CHR}:${START}-${END}"

# Extract reads from the specified region, filter by mapping quality and flags, and process with awk to prepare for R analysis
# -F 0x4: exclude unmapped reads
# -F 0x100: exclude secondary alignments
# -F 0x200: exclude reads that fail platform/vendor quality checks
# -F 0x400: exclude PCR or optical duplicates
# -F 0x800: exclude supplementary alignments
# The command keeps only primary, mapped, non-duplicate, high-quality alignments
samtools view -T GRCh38_full_analysis_set_plus_decoy_hla.fa -F 0x4 -F 0x100 -F 0x200 -F 0x400 -F 0x800 "$FILE" ${CHR}:${START}-${END} \
| awk -v s=$STARTSEQ -v e=$ENDSEQ -v iid=$ID '
BEGIN {STARTSEQ_TOADD = substr(s,7,3); ENDSEQ_TOADD=substr(e,1,3)}
$5 >= 30 && $7=="=" {
start=split($10,a,s);end=split($10,b,e); 
if (start == 2 && end == 2) {
split(a[2],segMID,e);
print iid,(and($2,16)==0?"forward":"reverse"),length(STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD), STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD,$10,$11
}
}' > "${OUT_DIR}/IID_${ID}_${GENE}_${REP}.txt"

nREADS=$(cat "${OUT_DIR}/IID_${ID}_${GENE}_${REP}.txt" | wc -l)
echo "Number of reads in region: $nREADS"

# if [ $nREADS -gt 0 ]; then
# # Rscript short_somatic_perIndividual.R $ID $DAT
# fi

# rm IID_${ID}.txt
# Clean up index files more safely (only remove the one we created)
# [ -f "${FILE}.crai" ] && rm "${FILE}.crai"

