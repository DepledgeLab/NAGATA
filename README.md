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


## Running NAGATA
### Test command
```
python3 NAGATA.py -i test-dataset/Ad5-12h-24h.bed -o test-outs/ -n test-dataset/Ad5-12h-24h.polyA.fwd.rev.tsv -cf test-dataset/sam_seq_strand_cigar.txt
NEED TO CHANGE INPUT TO SORTED.BAM AND REMOVE REQUIREMENT FOR CIGAR FILE (HARDCODE INSTEAD)
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


## Data preparation
In order to reconstruct a transcriptome, NAGATA requires, at minimum, a sorted BAM file containing DRS reads aligned against a reference genome. For optimal performance, we recommend (i) filtering aligned sequence reads to retain only primary alignments, and (ii) using nanopolish to exclude alignments originating from reads without defined polyA sequences.

### Basecalling
We recommend that DRS basecalling is performed with Guppy v4.2.2 or higher and the following parameters:
```
guppy_basecaller -i /nanopore/fast5/data -s /output/dir -c rna_r9.4.1_70bps_hac.cfg -r --calib_detect --trim_strategy rna --reverse_sequence true 
```

### Alignment
In general (inc. for most DNA viruses), reads can be aligned to a reference genome using minimap2 with the following parameters:
```
minimap2 -ax splice -k14 -uf ref.fasta DRS.reads.fastq > all.sam
```
For many RNA viruses, particular those which form subgenomicRNAs (e.g. coronaviruses), we recommend alternative parameters optimized for detecting non-splice junctions e.g:
```
minimap2 -ax splice -k 8 -w 3 -g 30000 -G 30000 -C0 -uf --no-end-flt --splice-flank=no ref.fasta DRS.reads.fastq > all.sam
```
Following alignment, filtering to retain only primary alignments in a sorted BAM file should be be performed as follows:
```
samtools view -b -F2308 -o primary.aligned.bam all.sam
samtools sort -o primary.aligned.sorted.bam primary.aligned.bam
samtools index primary.aligned.sorted.bam
```

### Read filtering
To identify and remove reads that might introduce artefacts in NAGATA, we perform a filtering step prior to alignment. Here, nanopolish is used to index-link reads in their fastq and fast5 formats, and subsequently to estimate the poly(A) tail length present in each read. 
```
nanopolish index -d DRS.reads.fast5/ DRS.reads.fastq
nanopolish polya --threads="n" --reads=DRS.reads.fastq --bam=primary.aligned.bam --genome=ref.fasta > DRS.polyA.tsv
```




## Troubleshooting

