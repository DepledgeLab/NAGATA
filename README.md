# NAGATA
**Nanopore Guided Annotation of Transcriptome Architectures**

NAGATA uses Nanopore direct RNA sequencing reads aligned to a genome to produce a transcriptome annotation. NAGATA functions by parsing read alignments (sorted BAM files) to identify Transcription Units (TUs) by internally converting BAM file into BED12 followed by numerically sorting "start" and "end" positions and then grouping alignments with similar "start" and "end" co-ordinates. This is performed on a row-by-row basis with a new TU defined only if the alignment co-ordinates of a given row differ from the previous row by greater than user-defined threshold (20 nt for transcription start sites (TSS), 50 nt for cleavage and polyadenylation sites (CPAS)).

![schema](/docs/Schema.jpg)

## Installing NAGATA
### Requirements
Nanopolish v0.11.1 or higher\
MiniMap2 v2.15 or higher\
SAMtools v1.3 or higher\
BEDtools v2.26 or higher\
Python v3.7.0 or higher\

### Installing  with git
```
git clone https://github.com/DepledgeLab/NAGATA
```

### Setting a python environment
NAGATA requires multiple python packages and benefits from setting up a dedicted environment

```
python3 -m venv NAGATA
# Activate environment 
source NAGATA/bin/activate

# Install dependencies 
pip install seaborn scipy pandas numpy biopython matplotlib

# Deactivate environment 
deactivate
```

Conda enthusiasts may use the environment-setup.yml file and instructions detailed in [DRUMMER](https://github.com/DepledgeLab/DRUMMER)


## Running NAGATA
### Testing NAGATA
```
python3 NAGATA.py -i test-dataset/Ad5.combined.05.subsampled.sorted.bam -n test-dataset/Ad5.subsample.nanopolish.tsv -o test-outs
```
Using this command, NAGATA should identify XX transcripts (XX forward strand, X reverse strand) with a total runtime of < 1 min

### Required arguments
```
-i      --input_file            BAM file containing reads aligned to reference genome
-o      --output_directory      Output directory 
```
### Optional arguments (basic)
```
-p      --polya                     Nanopolish file containing poly(A) length estimates
-nt     --nanopolish_tag            How to filter nanopolish file: P - retain PASS reads only, N - retain only reads that do not PASS, A - retain all reads (default: P)
-s      --soft_clip_filter          Remove alignments with 5' soft-clipping values greater than specified value (default: 3)
-c      --CPAS_noise_filter         Ignore CPAS that have an abundance count < c (default: 20)
-t      --TSS_noise_filter          Ignore TSS that have an abundance count < t (default: 4)
-m      --min_transcript_abundance  minimum transcript abundance (default: 3)
-r1     --reference_bed_f           BED12 file containing existing annotation of forward strand (see XXX)
-r2     --reference_bed_r           BED12 file containing existing annotation of reverse strand (see XXX)
-d      --strand                    specify strand as '+', '-', or 'both' (default: both)
```
### Optional arguments (advanced)
```
-cg     --CPAS_clustering           Grouping value for clustering CPAS (default: 50)
-tg     --TSS_clustering            Grouping value for clustering TSS (default: 20)
-iso    --isoform_clustering        Grouping value for blockSize and blockStarts (default: 50)
-a      --TSS_abundance_per_TU      TSS abundance per transcription unit (default: 0.1)
-b      --blocksize_noise_filter    Prior to isoform deconvolution - filter out low abundant blocksize sums to prevent incorrect daisy chaining, (default: 3)
```


## ***NAGATA outputs***
Using the commands detailed in the "Running NAGATA" section on the test data should produce:
```
(1) Final_cluster.strand'.bed                    BED12 format files detailing all transcript isoforms identified
(2) Final_cluster.NAGATA.'strand'.gff3           Final_cluster.'strand'.bed converted into GFF3 files
(3) NAGATA-parameters.tsv                        List of all parameters used for this run
(4) Filtering-counts.txt                         Details on how many alignments/reads are being filtered at each step
(5) Final_cluster.precollapsed.'strand'.tsv      A precollapsed version of Final_cluster.'strand'.bed which 
```


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
#### Generating high-quality genome-level alignments and divining estimates of poly(A) tail lengths
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
To identify and remove reads that might introduce artefacts in NAGATA, we perform a filtering step prior to alignment. Here, nanopolish is used to index-link reads in their fastq and fast5 formats, and subsequently to estimate the length and quality of the poly(A) tail present in each read. By default, NAGATA only considers reads for which nanopolish reports the poly(A) tail as 'PASS'. Applying the nanopolish filter to retain only 'failed' reads (-nt N) may be useful for interrogating artefacts.
```
nanopolish index -d DRS.reads.fast5/ DRS.reads.fastq
nanopolish polya --threads="n" --reads=DRS.reads.fastq --bam=primary.aligned.bam --genome=ref.fasta > DRS.polyA.tsv
```
### Preparing BED12 annotations files from existing GTF/GFF3 files


If you are experiencing errors during conversion from GFF3 to BED12 then please validate you GFF3 file using a [GFF3validator](http://genometools.org/cgi-bin/gff3validator.cgi). 

## Parameter optimisation
In many cases, running NAGATA with default input values will produce a high quality transcriptome. However, several parameters are sensitive to depth, particularly -t, -c, and -m. Here we outline strategies for visually inspecting NAGATA outputs to assist with parameter optimisation.

### Choosing a value for -t/-c
The optimal values for these parameters may be influenced by sequencing depth and/or the overall transcriptome structure. For instance, the -t parameter which influences the TSS noise threshold prior to TSS grouping may lead to incorrectly grouping adjacent TSSs at too low of a value or eliminated lesser abundant TSSs at higher values. The -t parameter has a direct impact on the -tg parameter, which defines how adjacent TSS positions are grouped. At too low of a value, there is an increased chance NAGATA defines incorrect TSS values, while a high value may lead to the merging of adjacent TSSs. By initially visualizing the data, the user may identify the robust TSS/CPAS position, and define noise/clustering parameters appropriately. 

![NoiseFilter](/docs/NoiseFilter.jpg)

Note that the behaviour of -t/-c is linked to -tg/-cg. In the below example, (1) the abundance of reads with 5' alignments at the same position are counted and filtered according to the specified -t value to remove noise and identify clusters. Subsequently (2), the highest value within a cluster is identified and, (3) for each cluster, all positions within a range specified by -tg are corrected to the position of the most abundant position. The same behaviour applies to the -c/-cg flags, the only difference being this is performed at the 3' end of read alignments. 

![Collapse](/docs/Collapse.jpg)


### Choosing a value for -m


## Troubleshooting
Over time we will populate this section with advice on troubleshooting. For now, if you have a problem, please re-read the documentation carefully and if that doesn't help, raise an issue and we will respond as soon as possible.




## Wisdom
Im ta nating
