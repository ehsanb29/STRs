### samtools

The -F flag in samtools view excludes reads that have the specified SAM flag bits set. Multiple -F flags are combined with bitwise OR, so any read matching any of those bits is excluded.

The values here filter out:

Flag	Hex	Meaning
-F 0x4	4	Read is unmapped
-F 0x100	256	Secondary alignment
-F 0x200	512	Read fails platform/vendor quality checks
-F 0x400	1024	PCR or optical duplicate
-F 0x800	2048	Supplementary alignment
So the command keeps only primary, mapped, non-duplicate, high-quality alignments — the clean, uniquely-mapped reads you'd want for calling STR alleles.

This is combined with the $5 >= 30 filter in the awk block, which additionally requires mapping quality ≥ 30.


### awk command
Here's a step-by-step breakdown:

Variables passed in (-v)

s = $STARTSEQ — left flanking sequence
e = $ENDSEQ — right flanking sequence
iid = $ID — individual ID
BEGIN block

These small anchors are later added back onto the extracted allele to give it a tiny flanking context.

Row filter: $5 >= 30 && $7 == "="

$5 = MAPQ — keep reads with mapping quality ≥ 30
$7 = RNEXT — "=" means the mate maps to the same chromosome (proper pair)
Splitting the read sequence ($10) by the flanking sequences

split() returns the number of pieces. If start == 2 and end == 2, each flanking sequence appears exactly once in the read — i.e., the read spans the entire STR with unambiguous boundaries.

Extracting the STR allele

a[2] is everything after STARTSEQ. Splitting that by ENDSEQ gives segMID[1] = the sequence between the two flanks — the raw STR allele.

Output (one line per valid read)

Field	Value
1	Individual ID
2	Strand: forward if flag bit 16 unset, else reverse
3	Allele length (with 3-bp anchors added on each side)
4	Allele sequence: STARTSEQ_TOADD + STR + ENDSEQ_TOADD
5	Full read sequence ($10)
6	Base quality string ($11)
The and($2,16)==0 check uses the SAM flag bit 0x10 to determine strand direction.

In summary: the awk command finds all properly-paired, high-quality reads that contain both flanking sequences exactly once, extracts the STR repeat sequence between them, and records the allele length and sequence — forming the input for the downstream R somatic expansion analysis