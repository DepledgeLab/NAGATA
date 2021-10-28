import pandas as pd
import modules.nanopolish_filter as nanopolish_filter
import modules.TSS_soft_clip_filter as TSS_soft_clip_filter
import modules.TSS_grouping as TSS_grouping
import modules.TES_grouping as TES_grouping
from collections import Counter
import os
strand_map = {'+':'fwd','-':'rev'}

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
    bed12_np = pd.read_csv(input_file,sep = '\t',header = None,names = names_list)
    bed12_np = bed12_np[bed12_np['seq-name'].isin(retain_seqs)]
    bed12_np = noise_filter(bed12_np,3,'start')
    bed12_np = noise_filter(bed12_np,1,'end')
    return bed12_np

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
            # retain_clusters.append(clusters)
#     clust_filt_df = df[df['Final_clusters'].isin(retain_clusters)]
    return retain_clusters_df    

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input dataset')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="Bed12 file input location")
    requiredGrp.add_argument("-o",'--output_directory', required=False, help="output file location",default = '.',type = str)
    requiredGrp.add_argument("-eg",'--TES_grouping', required=False, help="Grouping value for clustering (TES)",default = 50.0,type = float)
    requiredGrp.add_argument("-sg",'--TSS_grouping', required=False, help="Grouping value for clustering (TSS)",default = 20.0,type = float)
    requiredGrp.add_argument("-n",'--nanopolish', required=True, help="nanopolish location output location (requires concatenated fwd and rev)")
    requiredGrp.add_argument("-cf",'--cigar_file', required=True, help="cigar file containing [read_name,['fwd|reverse' alignment],cigar]")
    requiredGrp.add_argument("-c",'--cigar_filter', required=False, help="cigar filter on TSS strand",default=12.5,type = float)
    requiredGrp.add_argument("-m",'--minimum_cluster_seqs', required=False, help="Retains clusters that have a minimum m (default 20)",default=20,type = float)
#     requiredGrp.add_argument("-s",'--strand', required=True, help="Strand of bed file")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_directory']
    TES_group_val = int(args['TES_grouping'])
    TSS_group_val = int(args['TSS_grouping'])
    nanopolish_path = args['nanopolish']
    cigar_file_path = args['cigar_file']
    cigar_filter_val = int(args['cigar_filter'])
    min_clust = int(args['minimum_cluster_seqs'])
    polyA = pd.read_csv(nanopolish_path,sep = '\t')
    raw_bed12 = pd.read_csv(input_file,sep = '\t',header = None,usecols = [3,5],names = ['readname','strand'])
    os.makedirs(output_file,exist_ok=True)
    final_counts = pd.DataFrame()
    full_cluster_df = pd.DataFrame()
    for strand in ['+','-']:
#         print('STRAND',strand)
        strand_sequences = nanopolish_filter.return_nanopolish_dedup(raw_bed12,polyA,strand)
        bed12_nano_dedup = read_in_bed(input_file,strand,strand_sequences)
#         print('Nanopolish FILTER',bed12_nano_dedup.shape)

        filter_cigar = TSS_soft_clip_filter.filter_sequences(cigar_file_path,cigar_filter_val,strand)
        cigar_retain_lst = bed12_nano_dedup[bed12_nano_dedup['seq-name'].isin(list(filter_cigar['sequence']))]
        TSS_df = TSS_grouping.TSS_grouping(cigar_retain_lst,TSS_group_val)
        final_TU_cluster_df,TES_counts = TES_grouping.TES_grouping(TSS_df,TES_group_val)
        final_filtered_clusters =  filter_clusters(final_TU_cluster_df,.25)
        
        filtered_min_clusters = [TU for TU,count in Counter(final_filtered_clusters['Final_clusters']).items() if count > min_clust]
        final_filtered_clusters = final_filtered_clusters[final_filtered_clusters['Final_clusters'].isin(filtered_min_clusters)]
        
        pd.DataFrame(final_filtered_clusters['Final_clusters'].value_counts()).to_csv(output_file+'/Final_cluster.counts.' + strand_map[strand]+'.txt',header =None,sep = '\t')
        final_filtered_clusters.to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.txt',sep ='\t',index = None)
        TES_counts.to_csv(output_file+'/counts.TES.' + strand_map[strand]+'.txt',sep ='\t',index = None)
        cluster_outs = output_file + '/clust_sequences-'+strand_map[strand]+'/'
        os.makedirs(cluster_outs)
        for i in filtered_min_clusters:
            current_df = final_filtered_clusters[final_filtered_clusters['Final_clusters'] == i]
            current_df[['seq-name']].to_csv(cluster_outs+i+'.txt',header = None,index = None)
        
#     if strand == '+':
#         names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
#     else: 
#         names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
#     return_nanopolish_dedup(raw_bed12,polyA)