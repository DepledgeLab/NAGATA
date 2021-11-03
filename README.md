# NAGATA
***Nanopore Augmented Genome and Transcriptome Annotation***

## Running NAGATA
#Test command
python3 NAGATA.py -i test-dataset/Ad5-12h-24h.bed -o test-outs/ -n test-dataset/Ad5-12h-24h.polyA.fwd.rev.tsv -cf test-dataset/sam_seq_strand_cigar.txt
### Required Flags
```

-i      A BED file in bed12 format
-o      Output directory 
-n      Path to file containing per read estimates of poly-A tail lengths from Nanopolish-polya
-cf     Path to file containing individual sequence names with corresponding strand information, and cigar string
```
