# NAGATA
***Nanopore Augmented Genome and Transcriptome Annotation***

NAGATA uses Nanopore direct RNA sequencing reads aligned to a genome to produce a transcriptome annotation.
## Generating neccessary files

# Bed file
```
minimap2 -ax splice -k14 -uf --secondary=no GENOMIC.fasta dRNA-READS.fastq > dRNA-READS.GENOMIC.sam
bamToBed -bed12 -i dRNA-READS.GENOMIC.sam > dRNA-READS.GENOMIC.sam.bed
```
## Running NAGATA
### Test command
```
python3 NAGATA.py -i test-dataset/Ad5-12h-24h.bed -o test-outs/ -n test-dataset/Ad5-12h-24h.polyA.fwd.rev.tsv -cf test-dataset/sam_seq_strand_cigar.txt
```
### Required Flags
```
-i      A BED file in bed12 format
-o      Output directory 
-n      Path to file containing per read estimates of poly-A tail lengths from Nanopolish-polya
-cf     Path to file containing individual sequence names with corresponding strand information, and cigar string
```
### Optional Flags
```
-eg     Grouping value (int) for TESs to be considered within same cluster (default: 50)
-sg     Grouping value (int) for TSSs to be considered within same cluster (default: 20)
-c      Soft clipping value to filter lower quality reads (default: 12.5)
-m      Minimum number of reads a cluster must have to be considered real (default: 20)
```
