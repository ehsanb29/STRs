# This script is going to read two list of BAM files 
# run samtools view on each file and get the number of reads
# report the number of reads for all file in a table
# file names are like this: NDAR_INVAA524AZX_realigned_sorted.bam
# ID would be this: INVAA524AZX 
# sample ids exists in both lists, but the number of reads may be different, so we will report the number of reads for each list separately 
# the output table will have the following columns: ID, L1_READS, L2_READS
set -euo pipefail
LIST1="$1";
LIST2="$2";
OUT_DIR=$3
#  output dir exists check na make it
if [ ! -d "$OUT_DIR" ]; then
    mkdir -p "$OUT_DIR"
fi

# Process lists and write the numbers of reads for each list in one column
# sample ids exists in both lists, but the number of reads may be different, so we will report the number of reads for each list separately 
# create a temporary file to store the combined list of BAM files
COMBINED_LIST=$(mktemp)
# combine the two lists into one file
cat "$LIST1" "$LIST2" > "$COMBINED_LIST"
# iterate over each BAM file in the combined list
while IFS= read -r FILE; do
    # extract the ID from the file name 
    ID=$(basename "$FILE" | sed 's/NDAR_//' | sed 's/_realigned_sorted\.bam$//' | sed 's/\.bam$//')
    # get the number of reads in the BAM file
    echo $(date) "Processing $FILE with ID $ID"
    NREAD=$(samtools view -c "$FILE")
    echo $(date) "Number of reads in $(basename "$FILE"): $NREAD"
    # output the result to a table: ID, L1_READS, L2_READS
    if grep -q "$FILE" "$LIST1"; then
        echo -e "${ID}\t${NREAD}\t0" >> "${OUT_DIR}/read_counts.tsv"
    else
        echo -e "${ID}\t0\t${NREAD}" >> "${OUT_DIR}/read_counts.tsv"
    fi
done < "$COMBINED_LIST"

# remove the temporary combined list
rm "$COMBINED_LIST"
