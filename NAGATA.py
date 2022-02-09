import pandas as pd
import modules.nanopolish_filter as nanopolish_filter
import modules.TSS_soft_clip_filter as TSS_soft_clip_filter
import modules.TSS_grouping as TSS_grouping
import modules.TES_grouping as TES_grouping
import modules.isoform_grouping as iso_grouping
import modules.help_functions as help_functions
from collections import Counter
import os
import subprocess

import sys

strand_map = {'+':'fwd','-':'rev'}

def create_tmp_files(input_bam:'string',output_location:'string'):
    """Takes in an imput file and creates a temp file used for actual processing
    """
    os.makedirs(output_location,exist_ok=True)
    ## Samtools to get read-ids, cigar strings, orientation into a temporary file
    samtools_p1 = subprocess.Popen(f"samtools view {input_bam}".split(' '),stdout=subprocess.PIPE)
    fout_samtools = open(output_location +'/seq-cigar-orient.tmp', 'wb')
    samtools_p2 = subprocess.run(['cut','-f1,2,6'], stdin=samtools_p1.stdout, stdout=fout_samtools)
    ## Bedtools to convert BAM file to BED
    fout_bedtools = open(output_location +'/BED-file.tmp','wb')
    bedtools_p1 = f"bedtools bamtobed -bed12 -i {input_bam}".split(" ")
    bedtools_p2 = subprocess.run(bedtools_p1,stdout=fout_bedtools)
    print('temp files created')

def noise_filter(df,noise_value,start_or_end):
    """Takes in dataframe and filters start or end by the occurence 
    noise_value = 1, start_or_end = 'start', filters for 'start' that has occurrence > 1
    """
    filtered_positions = [positions for positions,occurences in Counter(df[start_or_end]).items() if occurences > noise_value]
    filtered_df = df[df[start_or_end].isin(filtered_positions)]
    return filtered_df
    
def read_in_bed(bedfile,strand,retain_seqs):
    if strand == '+':
        names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    else: 
        names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed12_np = pd.read_csv(bedfile,sep = '\t',header = None,names = names_list)
    original_bed_size =  bed12_np.shape[0]
    bed12_np = bed12_np[bed12_np['seq-name'].isin(retain_seqs)]
    print(f'Reads lacking readable polyA tail:\t',original_bed_size-bed12_np.shape[0])
    polyA_reads = bed12_np.shape[0]
    bed12_np = noise_filter(bed12_np,3,'start')
    bed12_np = noise_filter(bed12_np,1,'end')
    print(f'Reads rejected by alignment location:\t',polyA_reads-bed12_np.shape[0])
    return bed12_np

def create_bedgraph(input_bam,tmp_folder,outfile,TU_iso,TU_iso_read_ids):
    temp_output_location = tmp_folder + '/' + TU_iso + '.bam'
    infile = pysam.AlignmentFile(input_bam)
    tmp_out = pysam.AlignmentFile(temp_output_location, template= infile, mode= 'wb')
    for aln in infile:
        if aln.query_name in TU_iso_read_ids:
            tmp_out.write(aln)
    tmp_out.close()
    infile.close()
    samtools_sort = f"samtools sort -o {temp_output_location}.sorted.bam {temp_output_location}".split(' ')
    subprocess.run(samtools_sort)
    
    samtools_idx = f"samtools index {temp_output_location}.sorted.bam".split(' ')
    subprocess.run(samtools_idx)
    fout_bedgraph = open(outfile,'w')

    genomecoverage_p2 = subprocess.run(f"bedtools genomecov -ibam {temp_output_location}.sorted.bam -bg -split".split(' '),stdout=fout_bedgraph)


def filter_clusters(df,min_abundance):
    '''Takes in dataframe, goes through the clusters and filters using TSS-abundance% > min_abundance
    '''
    retain_clusters_df = pd.DataFrame()
    for clusters in df['Final_clusters'].unique():
        current_df = df[df['Final_clusters'] == clusters]
        top_TSS = current_df['start'].value_counts().head(1).index[0]
        top_count = list(current_df['start']).count(top_TSS)
        total = sum(Counter(list(current_df['start'])).values())
        if top_count/total > min_abundance:
            keep_unique_start_df = current_df[current_df['start']==top_TSS]
            retain_clusters_df = pd.concat([retain_clusters_df,keep_unique_start_df])
    return retain_clusters_df    
def get_individual_clusters(full_df):
    """Takes in the full dataframe and collapses its
    """
    cluster_indv = pd.DataFrame()
    for clust in set(full_df['TU.#-iso.#']):
        current_cluster = full_df[full_df['TU.#-iso.#']==clust]
        top_blocksize_count = current_cluster.groupby('blockSize-sums').size().sort_values(ascending = False).head(1).index[0]
        single_row_df = current_cluster[current_cluster['blockSize-sums'] == top_blocksize_count].sample(1)
        cluster_indv = pd.concat([cluster_indv,single_row_df])
    return cluster_indv
def formatting_output(cluster_indv):
    """Returns a formatted bed12 file
    """
    cols = list(cluster_indv)
    cols[3], cols[13] = cols[13], cols[3]
    cols[4], cols[14] = cols[14], cols[4]
    cluster_indv = cluster_indv.loc[:,cols]
    cluster_indv = cluster_indv[cluster_indv.columns[:len(cluster_indv.columns)-3]]
    cluster_indv.rename(columns={'TU.#-iso.#': 'name', 'TU.#-iso.#-count': 'score'}, inplace=True)
    return cluster_indv
if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser(description = 'Takes in input dataset')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="BAM file input location")
    requiredGrp.add_argument("-o",'--output_directory', required=False, help="output file location",default = '.',type = str)
    requiredGrp.add_argument("-eg",'--TES_grouping', required=False, help="Grouping value for clustering (TES), default 50",default = 50.0,type = float)
    requiredGrp.add_argument("-sg",'--TSS_grouping', required=False, help="Grouping value for clustering (TSS), default 20",default = 20.0,type = float)
    requiredGrp.add_argument("-n",'--nanopolish', required=True, help="nanopolish location output location (requires concatenated fwd and rev)")
    requiredGrp.add_argument("-c",'--cigar_filter', required=False, help="cigar filter on TSS strand, default 12",default=12,type = float)
    requiredGrp.add_argument("-m",'--minimum_cluster_seqs', required=False, help="Retains clusters that have a minimum m (default 20)",default=20,type = float)
    requiredGrp.add_argument("-pc",'--post_clustering_filter', required=False, help="Within a TU, if the top TSS does not have at least (pc) of the reads, this cluster is ignored (default .25)",default=.25,type = float)
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_directory']
    TES_group_val = int(args['TES_grouping'])
    TSS_group_val = int(args['TSS_grouping'])
    nanopolish_path = args['nanopolish']
#    cigar_file_path = args['cigar_file']
    cigar_filter_val = int(args['cigar_filter'])
    min_clust = int(args['minimum_cluster_seqs'])
    pc_filter = float(args['post_clustering_filter'])
#     help_functions.check_bedtools()
#     if help_functions.check_samtools() == True and help_functions.check_bedtools() == True:
    create_tmp_files(input_file,output_file +'/tmp')
    cigar_file_path = output_file +'/tmp/' + 'seq-cigar-orient.tmp'
    polyA = pd.read_csv(nanopolish_path,sep = '\t')
    raw_bed12 = pd.read_csv(output_file +'/tmp/' + 'BED-file.tmp',sep = '\t',header = None,usecols = [3,5],names = ['readname','strand'])
    #os.makedirs(output_file,exist_ok=True)
    final_counts = pd.DataFrame()
    full_cluster_df = pd.DataFrame()
    for strand in ['+','-']:
#         print('STRAND',strand)
#         print(f'Currently processing {strand_map[strand]} strand...')
        os.system(f'echo Currently processing {strand_map[strand]} strand...')
        sys.stdout = open(output_file +'/' +f'Filtering-counts.{strand_map[strand]}.tsv', 'w')

        strand_sequences = nanopolish_filter.return_nanopolish_dedup(raw_bed12,polyA,strand)
        print(f'Total input reads:\t',raw_bed12.shape[0])
        bed12_nano_dedup = read_in_bed(output_file +'/tmp/' + 'BED-file.tmp',strand,strand_sequences)
#         print('Nanopolish FILTER',bed12_nano_dedup.shape)
        
        filter_cigar = TSS_soft_clip_filter.filter_sequences(cigar_file_path,cigar_filter_val,strand)
        cigar_retain_lst = bed12_nano_dedup[bed12_nano_dedup['seq-name'].isin(list(filter_cigar['sequence']))]
        print(f'Reads with soft-clipping exceeding filter:\t',bed12_nano_dedup.shape[0]-cigar_retain_lst.shape[0])
        TSS_df = TSS_grouping.TSS_grouping(cigar_retain_lst,TSS_group_val)
        final_TU_cluster_df,TES_counts = TES_grouping.TES_grouping(TSS_df,TES_group_val)
#         print('PC-FILTER',pc_filter)
        final_filtered_clusters =  filter_clusters(final_TU_cluster_df,pc_filter)
        print(f'Reads available for transcriptome construction:\t',final_filtered_clusters.shape[0],'\n')


        print(f'Number of transcription units identified:\t',len(set(final_filtered_clusters['Final_clusters'])))
        filtered_min_clusters = [TU for TU,count in Counter(final_filtered_clusters['Final_clusters']).items() if count > min_clust]
        final_filtered_clusters = final_filtered_clusters[final_filtered_clusters['Final_clusters'].isin(filtered_min_clusters)]
        print(f'Number of transcription units passing filter:\t',len(set(final_filtered_clusters['Final_clusters'])),'\n')
        final_filtered_clusters['blockSize-sums'] = iso_grouping.get_blocksize_length(final_filtered_clusters['blockSizes'])  
        final_filtered_clusters_iso = iso_grouping.run_isoform(final_filtered_clusters)
#         print(final_filtered_clusters.head())
        print(f'Number of isoforms identified:\t',len(set(final_filtered_clusters_iso['TU.#-iso.#'])))
        final_filtered_clusters_iso_bed = get_individual_clusters(final_filtered_clusters_iso)
        final_filtered_clusters_iso_bed = formatting_output(final_filtered_clusters_iso_bed)
        print(f'Number of isoforms passing filter:\t',final_filtered_clusters_iso_bed.shape[0])
#         pd.DataFrame(final_filtered_clusters['Final_clusters'].value_counts()).to_csv(output_file+'/Final_cluster.counts.' + strand_map[strand]+'.txt',header =None,sep = '\t')
        final_filtered_clusters_iso_bed.to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.bed',sep ='\t',index = None,header = None)
        final_filtered_clusters_iso.to_csv(output_file+'/Final_cluster.precollapsed.' + strand_map[strand]+'.tsv',sep ='\t',index = None)
        bedgraph_outs = output_file+'/bedgraphs-'+strand_map[strand]
        os.makedirs(bedgraph_outs,exist_ok=True)
        sys.stdout.close()
        import pysam
        for TU_isoforms in set(final_filtered_clusters_iso['TU.#-iso.#']):
            os.system(f'echo {TU_isoforms}')
            current_seqs = list(final_filtered_clusters_iso[final_filtered_clusters_iso['TU.#-iso.#'] == TU_isoforms]['seq-name'])
            #create_bedgraph(infile,tmp_folder,outfile,TU_iso,TU_iso_read_ids)
            create_bedgraph(input_file,output_file +'/tmp',bedgraph_outs+'/'+TU_isoforms +'.'+strand_map[strand]+ '.bedgraph',TU_isoforms,current_seqs)
#         TES_counts.to_csv(output_file+'/counts.TES.' + strand_map[strand]+'.txt',sep ='\t',index = None)
#         cluster_outs = output_file + '/clust_sequences-'+strand_map[strand]+'/'
#         os.makedirs(cluster_outs)
#         for i in filtered_min_clusters:
#             current_df = final_filtered_clusters[final_filtered_clusters['Final_clusters'] == i]
#             current_df[['seq-name']].to_csv(cluster_outs+i+'.txt',header = None,index = None)
      
        
#     if strand == '+':
#         names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
#     else: 
#         names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
#     return_nanopolish_dedup(raw_bed12,polyA)