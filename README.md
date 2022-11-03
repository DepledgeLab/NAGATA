# NAGATA
***Nanopore Guided Annotation of Transcriptome Architectures***

NAGATA uses Nanopore direct RNA sequencing reads aligned to a genome to produce a transcriptome annotation.
## Data Preparation
### Required software
Nanopolish v0.11.1 or higher\
MiniMap2 v2.15 or higher\
SAMtools v1.3 or higher\
BEDtools v2.26 or v2.27

### ***BAM file***

#### BAM files are used to cluster sequences based on alignment similarities in an iterative manner. 
NAGATA parses read alignments to identify Transcriptional Units (TUs) by internally converting BAM file into BED12 followed by numerically sorting "start" and "end" positions and then grouping alignments with similar "start" and "end" co-ordinates. This is performed on a row-by-row basis with a new TU defined only if the alignment co-ordinates of a given row differ from the previous row by greater than user-defined threshold (20 nt for transcription start sites (TSS), 50 nt for cleavage and polyadenylation sites (CPAS)).  TSS threshold can be tuned using the -sg flag while the CPAS threshold can be control using the -eg flag.

![TSS-example](/modules/TSS-example.png)
![Algorithm example](/modules/Grouping-TSS.pdf)

Once clustered have been identifed using "starts", NAGATA uses a similar algorithm using "ends" within each cluster.
![TES-example](/modules/TES-example.png)

## ***NAGATA outputs***
Using the commands detailed in the "Running NAGATA" section, running NAGATA using the test data should produce:
(1) Final_cluster.fwd.bed file with 22 distinct transcripts
(2) Final_cluster.rev.bed file with 4 distinct transcripts 
among the rest of the files indicated below with a total runtime of <1min

```
In the directory specified by the -o flag, the following files should be produced for each available strand

1: Filtering-counts.'strand'.tsv                How many reads are being filtered at each filter
2: Final_cluster.'strand'.bed                   Main output of NAGATA, regular BED file with 4th column being the NAGATA-name and 5th column being the transcript abundance
3: Final_cluster.NAGATA.'strand'.gff3           Final_cluster.'strand'.bed converted in GFF3 file
4: Final_cluster.precollapsed.'strand'.tsv      A precollapsed version of Final_cluster.'strand'.bed which 
5: NAGATA-parameters.tsv                        List of all parameters used for this run
6: TMP directory                                For testing purposes - intermediate files are saved here
```

## Running NAGATA
### Environment setup
For now creating an environment using directions detailed in [DRUMMER](https://github.com/DepledgeLab/DRUMMER) should suffice. 

### Test command
```
python3 NAGATA.py -i test-dataset/Ad5.combined.05.subsampled.sorted.bam -n test-dataset/Ad5.subsample.nanopolish.tsv -o test-outs

```
### Required Flags
```
-i      BAM file containing reads aligned to reference genome
-o      Output directory 
-n      Path to file containing per read estimates of poly-A tail lengths from Nanopolish-polya
```
### Optional Flags
```
-eg     Grouping value (int) for TESs to be considered within same cluster (default: 50)
-sg     Grouping value (int) for TSSs to be considered within same cluster (default: 20)
-c      Soft clipping value to filter lower quality reads (default: 12.5)
-m      Minimum number of reads a Transcriptional unit (TU) must have to be considered "real" (default: 20)
-mi     Minimum number of reads a transcript isoforms needs to be considered "real" (default: 5)
-pc     Within a TU, if the most abundant TSS does not have at least (pc) of reads, this TU is ignored (default: .25)
-TSS_pad     TSS peak is defined by most abundant TSS in TU, (+/- TSS_pad) are used to include surrounding reads and subsequently corrected for more reads(default 12)
TES_pad      CPAS peak is defined by most abundant CPAS in TU, (+/- TES_pad) are used to include surrounding reads and subsequently corrected for more reads(default 25)
-b     Output bedgraph coverage data for each transcript isoform to file (default: False) 
-s     Initially filter out TSS sites that have a count <s (default: 3)
-p     Initially filter out CPAS sites that have a count <c (default: 1)
-f     Value to distinguish isoforms of the same TU, i.e. if sizes of transcripts (of the same exon count) fall within f, these transcripts are considered the same (default: 10)
-t     How to apply nanopolish filter, retain reads P (PASS), N (NO_PASS), A (N + P, all reads) (default: P)
```


## Data preparation
In order to reconstruct a poly(A)+ transcriptome, NAGATA requires, at minimum, a sorted BAM file containing DRS reads aligned against a reference genome. For optimal performance, we recommend (i) filtering aligned sequence reads to retain only primary alignments, and (ii) using nanopolish polya to exclude alignments originating from reads without well-supported polyA sequences.

### Alignment
In general (inc. for most DNA viruses), reads can be aligned to a reference genome using minimap2 with the following parameters:
```
minimap2 -ax splice -k14 -uf ref.fasta DRS.reads.fastq > all.sam
```
For RNA viruses which generate subgenomic RNAs by discountinuous transcription (e.g. coronaviruses), we  we recommend alternative parameters for alignment that are optimized for detecting non-canonical junctions
```
minimap2 -ax splice -k 8 -w 3 -g 30000 -G 30000 -C0 -uf --no-end-flt --splice-flank=no ref.fasta DRS.reads.fastq > all.sam
```

### Filtering SAM files to retain only primary alignments (alignment flag 0 [forward strand] and 16 [reverse strand]) prior to co-ordinate sorting and indexing
```
samtools view -b -F2308 -o primary.aligned.bam all.sam
samtools sort -o primary.aligned.sorted.bam primary.aligned.bam
samtools index primary.aligned.sorted.bam
```

### Filtering for reads with robust poly(A) tails
To identify and remove reads that might introduce artefacts in NAGATA, we perform a filtering step prior to alignment. Here, nanopolish is used to index-link reads in their fastq and fast5 formats, and subsequently to estimate the length and quality of the poly(A) tail present in each read. By default, NAGATA only considers reads for which nanopolish reports the poly(A) tail as 'PASS'
```
nanopolish index -d DRS.reads.fast5/ DRS.reads.fastq
nanopolish polya --threads="n" --reads=DRS.reads.fastq --bam=primary.aligned.bam --genome=ref.fasta > DRS.polyA.tsv
```

## Troubleshooting







## Wisdom
Im ta nating
