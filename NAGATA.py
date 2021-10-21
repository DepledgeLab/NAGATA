import pandas as pd
import modules.nanopolish_filter as nanopolish_filter
import modules.TSS_soft_clip_filter as TSS_soft_clip_filter
import modules.TSS_grouping as TSS_grouping
import modules.TES_grouping as TES_grouping
from collections import Counter
strand_map = {'+':'fwd','-':'rev'}
def read_in_bed(bedfile,strand,retain_seqs):
    if strand == '+':
        names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    else: 
        names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed12_np = pd.read_csv(input_file,sep = '\t',header = None,names = names_list)
    bed12_np = bed12_np[bed12_np['seq-name'].isin(retain_seqs)]
    return bed12_np
    

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
        
        filtered_min_clusters = [TU for TU,count in Counter(final_TU_cluster_df['Final_clusters']).items() if count > min_clust]
        final_filtered_clusters = final_TU_cluster_df[final_TU_cluster_df['Final_clusters'].isin(filtered_min_clusters)]
        
        final_filtered_clusters.to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.txt',sep ='\t',index = None)
        TES_counts.to_csv(output_file+'/counts.TES.' + strand_map[strand]+'.txt',sep ='\t',index = None)
        
#     if strand == '+':
#         names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
#     else: 
#         names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
#     return_nanopolish_dedup(raw_bed12,polyA)