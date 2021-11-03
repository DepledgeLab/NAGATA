# NAGATA
***Nanopore Augmented Genome and Transcriptome Annotation***

NAGATA uses Nanopore direct RNA sequencing reads aligned to a genome to produce a transcriptome annotation.
## Generating neccessary files

### ***Bed file***

#### The Bed file is used to cluster sequences based on alignment similarity in an iterative manner. 
To identify Transcriptional units, after an inital noise filter, NAGATA numerically sorts "start" positions and identifies row-by-row if the difference between current value and previous is < than the grouping value (sg), if this is true, the current sequence is assigned to this cluster. If the value is larger, then a new cluster is created and the process is repeated. 
![TSS-example](/modules/TSS-example.png)
![Algorithm example](/modules/Grouping-TSS.pdf)

Once clustered have been identifed using "starts", NAGATA uses a similar algorithm using "ends" within each cluster.
![TES-example](/modules/TES-example.png)
##### Generating Bed file
```
minimap2 -ax splice -k14 -uf --secondary=no "genomic".fasta "dRNA-READS".fastq > "dRNA-READS.GENOMIC".sam
bamToBed -bed12 -i "dRNA-READS.GENOMIC".sam > "dRNA-READS.GENOMIC".sam.bed
```

### ***Cigar file***
#### A cigar file is used to filter out 5' misaligned sequences prior to clustering. This step is important in preventing false TSSs in the final annotation.
##### Generating Cigar file
```
samtools view "dRNA-READS.GENOMIC".sam | cut -f1,2,6 > "SEQUENCE.STRAND.CIGAR".txt
```
### ***Nanopolish file***
#### A nanopolish file containing per read estimates of poly-A tails is used to filter out reads with incorrectly mapped 3' ends due to technical errors.
##### Generating Nanopolish file
```
nanopolish index -d "FAST5-READS" "dRNA-READS".fastq
nanopolish polya --threads=8 --reads="dRNA-READS".fastq --bam="dRNA-READS.GENOMIC".bam --genome="genomic".fasta > "dRNA-READS".polyA.tsv
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
