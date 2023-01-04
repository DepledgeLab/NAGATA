# NAGATA
**Nanopore Guided Annotation of Transcriptome Architectures**

NAGATA uses Nanopore direct RNA sequencing reads aligned to a genome to produce a transcriptome annotation.

## Required software
Nanopolish v0.11.1 or higher\
MiniMap2 v2.15 or higher\
SAMtools v1.3 or higher\
BEDtools v2.26 or v2.27

## Required data inputs 
####(for detailed information on how to generate required data inputs, please read the data preparation section)
### ***BAM file***
#### BAM files are used to cluster sequences based on alignment similarities in an iterative manner. 
NAGATA parses read alignments to identify Transcriptional Units (TUs) by internally converting BAM file into BED12 followed by numerically sorting "start" and "end" positions and then grouping alignments with similar "start" and "end" co-ordinates. This is performed on a row-by-row basis with a new TU defined only if the alignment co-ordinates of a given row differ from the previous row by greater than user-defined threshold (20 nt for transcription start sites (TSS), 50 nt for cleavage and polyadenylation sites (CPAS)).  TSS threshold can be tuned using the -sg flag while the CPAS threshold can be control using the -eg flag.

![TSS-example](/modules/TSS-example.png)
![Algorithm example](/modules/Grouping-TSS.pdf)

Once clustered have been identifed using "starts", NAGATA uses a similar algorithm using "ends" within each cluster.
![TES-example](/modules/TES-example.png)


## Running NAGATA
### Environment setup
For now creating an environment using directions detailed in [DRUMMER](https://github.com/DepledgeLab/DRUMMER) should suffice. 

### Test command
```
python3 NAGATA.py -i test-dataset/Ad5.combined.05.subsampled.sorted.bam -n test-dataset/Ad5.subsample.nanopolish.tsv -o test-outs

```
### Required Flags
```
-i      --input_file            BAM file containing reads aligned to reference genome
-o      --output_directory      Output directory 
-p      --polya                 Nanopolish output file containing poly(A) length estimates
```
### Optional Flags
```
-c      --CPAS_noise_filter     Ignore CPAS that have an abundance count < c (default: XX)
-cg     --CPAS_clustering       Grouping value (int) for CPASs to be considered within same cluster (default: 50)
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


## ***NAGATA outputs***
Using the commands detailed in the "Running NAGATA" section on the test data should produce:

(1) Final_cluster.strand'.bed                    BED12 format file detailing all transcript isoforms identified
(2) Final_cluster.NAGATA.'strand'.gff3           Final_cluster.'strand'.bed converted in GFF3 file
(3) NAGATA-parameters.tsv                        List of all parameters used for this run
(4) Filtering-counts.txt                         Details on how many alignments/reads are being filtered at each step
(5) Final_cluster.precollapsed.'strand'.tsv      A precollapsed version of Final_cluster.'strand'.bed which 

When used with the test dataset and default parameters, NAGATA should identify 26 transcripts (22 forward strand, 4 reverse strand) with a total runtime of < 1 min
```
NAGATA also produces a series of intermediary files to aid in optimisation/troubleshooting. These are stored in a tmp/ directory within the main output directory specified by the -o flag.

1.raw-alignment.bed
2.filter_nanopolish.'strand'.bed 
3.filter_cigar.'strand'.bed
4.CPAS-grouping.'strand'.bed
5.TSS-grouping.'strand'.bed 
6.columns-fully-corrected.'strand'.bed 
7.Isoform-deconvolution.'strand'.bed 

Correct.blocksizes.'strand'.bed
Correct.TSS.CPAS.'strand'.bed
parsed.cigar.'strand'.tsv
Passed_NANOPOLISH.bed
seq-cigar-orient.tmp
```


## Data preparation 
#### (generating high-quality genome-level alignments and divining estimates of poly(A) tail lengths)
In order to construct a transcriptome annotation, NAGATA requires, at minimum, a sorted BAM file containing DRS reads aligned against a reference genome. For optimal performance, we recommend (i) filtering aligned sequence reads to retain only primary alignments, and (ii) using nanopolish polya to exclude alignments originating from reads without well-supported polyA sequences.

### Alignment
In general (inc. for most DNA viruses), reads can be aligned to a reference genome using minimap2 with the following parameters:
```
minimap2 -ax splice -k14 -uf ref.fasta DRS.reads.fastq > all.sam
```
For RNA viruses which generate subgenomic RNAs by discountinuous transcription (e.g. coronaviruses), we recommend alternative parameters for alignment that are optimized for detecting non-canonical junctions
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


## Parameter optimisation
In most cases, running NAGATA with default input values will produce a high quality transcriptome. However, several parameters are sensitive to depth, particularly -t, -c, and -m. Here we outline strategies for visually inspecting NAGATA outputs to assist with parameter optimisation.

### Choosing a value for -t

### Choosing a value for -c 

### Choosing a value for -m


## Troubleshooting
Over time we will populate this section with advice on troubleshooting. For now, if you have a problem, please re-read the documentation carefully and if that doesn't help, raise an issue and we will respond as soon as possible.




## Wisdom
Im ta nating
