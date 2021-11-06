# NAGATA
***Nanopore Augmented Genome and Transcriptome Annotation***

NAGATA uses Nanopore direct RNA sequencing reads aligned to a genome to produce a transcriptome annotation.
## Data Preparation
### Required software
Nanopolish v0.11.1 or higher
MiniMap2 v2.15 or higher
SAMtools v1.3 or higher
BEDtools v2.26 or v2.27

### ***BED12 file***

#### BED12 files are used to cluster sequences based on alignment similarities in an iterative manner. 
NAGATA parses read alignments to identify Transcriptional Units (TUs) by numerically sorting "start" and "end" positions and then grouping alignments with similar "start" and "end" co-ordinates. This is performed on a row-by-row basis 

with a new TU defined only if the alignment co-ordinates of a given row differ from the previous row by greater than user-defined threshold (25 nt for transcription start sites (TSS), 50 nt for cleavage and polyadenylation sites (CPAS)). 

, after an inital noise filter, NAGATA numerically sorts "start" positions and identifies row-by-row if the difference between current value and previous is < than the grouping value (sg), if this is true, the current sequence is assigned to this cluster. If the value is larger, then a new cluster is created and the process is repeated. 
![TSS-example](/modules/TSS-example.png)
![Algorithm example](/modules/Grouping-TSS.pdf)

Once clustered have been identifed using "starts", NAGATA uses a similar algorithm using "ends" within each cluster.
![TES-example](/modules/TES-example.png)
##### Generating Bed file
```
minimap2 -ax splice -k14 -uf --secondary=no "genomic".fasta "dRNA-READS".fastq > "dRNA-READS.GENOMIC".sam
bamToBed -bed12 -i "dRNA-READS.GENOMIC".sam > "dRNA-READS.GENOMIC".sam.bed
```

### ***CIGAR Report***
#### CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings for each read alignment are extracted and supplied to NAGATA to identify putative 5' alignment artefacts resulting from splice junctions. This step is critical for preventing artefact TSS identification.
note - can we not automate this part? - Yes. I will do it in a future version. This would require NAGATA to accept a BAM file instead of a BED file.

##### Generating CIGAR Report file
```
samtools view "dRNA-READS.GENOMIC".sam | cut -f1,2,6 > "SEQUENCE.STRAND.CIGAR".txt
```


### ***Nanopolish poly(A) length estimation file***
#### Nanopolish is used first to index-link reads in their fastq and fast5 formats, and subsequently to estimate the poly(A) tail length present in each read. NAGATA subsequently filters these data to exclude reads for which poly(A) tail lengths could not be reported as these produce incorrect 3' end alignments.
##### Generating Nanopolish file
```
nanopolish index -d "dRNA-FAST5-READS" "dRNA-FASTQ-READS".fastq
nanopolish polya --threads="n" --reads="dRNA-READS".fastq --bam="dRNA-READS.GENOMIC".bam --genome="genomic".fasta > "dRNA-READS".polyA.tsv
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

## Troubleshooting

