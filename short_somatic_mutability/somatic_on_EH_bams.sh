set -euo pipefail
FILE="$1"; DAT=$2; OUT_DIR=$3

# Strip .bam, .cram, or .final.cram extension to get sample ID
ID=$( basename "$FILE" | sed 's/\.final\.cram$//' | sed 's/\.cram$//' | sed 's/\.bam$//' )

# Create output directory if it doesn't exist
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

echo "Processing individual: $ID with data: $DAT"

GENE=$(echo $DAT | cut -d'_' -f 1); REP=$(echo $DAT | cut -d'_' -f 7);
echo "Gene: $GENE with repeat type: $REP"
CHR=$(echo $DAT | cut -d'_' -f 2); START=$(echo $DAT | cut -d'_' -f 3); END=$(echo $DAT | cut -d'_' -f 4);
STARTSEQ=$(echo $DAT | cut -d'_' -f 5); ENDSEQ=$(echo $DAT | cut -d'_' -f 6);

echo "Region: $CHR:$START-$END with flanking sequences $STARTSEQ and $ENDSEQ"

# ExpansionHunter places reads at positions offset from the locus (potentially
# 1000+ bp away from the exact locus coordinates).  A buffer of 2000 bp on
# each side ensures all EH-realigned reads are captured by the region query.
# Arithmetic uses bash to avoid negative coordinates at chromosome start.
BUFFER=2000
QUERY_START=$(( START > BUFFER ? START - BUFFER : 1 ))
QUERY_END=$(( END + BUFFER ))

echo "samtools view $FILE ${CHR}:${QUERY_START}-${QUERY_END}"

# ExpansionHunter BAM-specific notes:
# - Reads overlapping the repeat are marked as unmapped (FLAG 0x4 set) with MAPQ=0
#   and CIGAR="*", so all standard mapping-quality / flag filters must be omitted.
# - RNEXT is "*" (not "="), so the $7=="=" guard used in somatic.sh is also dropped.
# - No -T reference flag: BAM files do not require a reference for decoding (unlike CRAM).
# - The raw read sequence ($10) is preserved by EH, so flanking-sequence matching
#   (STARTSEQ / ENDSEQ) works identically to the original pipeline.
samtools view "$FILE" ${CHR}:${QUERY_START}-${QUERY_END} \
| awk -v s=$STARTSEQ -v e=$ENDSEQ -v iid=$ID \
      -v locus="${GENE}" '
BEGIN {STARTSEQ_TOADD = substr(s,7,3); ENDSEQ_TOADD=substr(e,1,3)}
{
  # Keep only reads belonging to this locus via the XG tag
  # (guards against picking up reads from a nearby locus in the expanded window)
  xg_ok = 0
  for (i=12; i<=NF; i++) {
    if ($i ~ ("^XG:Z:" locus)) { xg_ok = 1; break }
  }
  if (!xg_ok) next

  start=split($10,a,s); end=split($10,b,e);
  if (start == 2 && end == 2) {
    split(a[2],segMID,e);
    print iid,(and($2,16)==0?"forward":"reverse"),length(STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD), STARTSEQ_TOADD segMID[1] ENDSEQ_TOADD,$10,$11
  }
}' > "${OUT_DIR}/IID_${ID}_${GENE}_${REP}.txt"

nREADS=$(wc -l < "${OUT_DIR}/IID_${ID}_${GENE}_${REP}.txt")
echo "Number of reads in region: $nREADS"

if [ $nREADS -gt 0 ]; then
    Rscript short_somatic_perIndividual_EH.R $ID $DAT $OUT_DIR
fi

# rm "${OUT_DIR}/IID_${ID}_${GENE}_${REP}.txt"

