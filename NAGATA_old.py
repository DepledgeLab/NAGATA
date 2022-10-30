import pandas as pd
import modules.nanopolish_filter as nanopolish_filter
import modules.TSS_soft_clip_filter as TSS_soft_clip_filter
import modules.TSS_grouping as TSS_grouping
import modules.TES_grouping as TES_grouping
import modules.isoform_grouping as iso_grouping
import modules.help_functions as help_functions
import modules.BED2GFF3 as bed_convert
from collections import Counter
import os
import subprocess

import sys
#python3 NAGATA.py -i /gpfs/home/ja3539/depledge/Jonathan/NAGATA/NAGATA-testing/Ad5-12h-24h.raw.AD5.genome.sorted.bam -n test-dataset/Ad5-12h-24h.polyA.fwd.rev.tsv -o isoform-filter_test/$i.isoform_results
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
    fout_bedtools = open(output_location +'/raw-alignment.2.bed','wb')
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
    
def read_in_bed(bedfile,strand,retain_seqs,noise_filt_TSS,noise_filt_CPAS):
    """Input: Takes in raw bed file, current strand, list of sequences that (1. pass nanopolish filter, 2. Not duplicated).
    Does further filtering for 1. Each TSS must have an abundance of >3, 2. Each TES must have an abundance of > 1 
    """
    if strand == '+':
        names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    else: 
        names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    bed12_np = pd.read_csv(bedfile,sep = '\t',header = None,names = names_list)
    original_bed_size =  bed12_np.shape[0]
    bed12_np = bed12_np[bed12_np['seq-name'].isin(retain_seqs)]
    print(f'Reads lacking readable polyA tail:\t',original_bed_size-bed12_np.shape[0])
    polyA_reads = bed12_np.shape[0]
    bed12_np = noise_filter(bed12_np,noise_filt_TSS,'start')
    bed12_np = noise_filter(bed12_np,noise_filt_CPAS,'end')
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


def filter_clusters(df,min_abundance,pad):
    '''Takes in dataframe, goes through the clusters and filters using TSS-abundance% > min_abundance
    '''
    retain_clusters_df = pd.DataFrame()
    for clusters in df['Final_clusters'].unique():
        current_df = df[df['Final_clusters'] == clusters]
        top_TSS = current_df['start'].value_counts().head(1).index[0]
        top_count = list(current_df['start']).count(top_TSS)
        total = sum(Counter(list(current_df['start'])).values())
#         keep_unique_start_df = current_df[(current_df['start'] < current_df['start'].max()) & (current_df['start'] > current_df['start'].min())]
        
#         keep_unique_start_df = current_df[(current_df['start'] < top_TSS+pad) & (current_df['start'] > top_TSS-pad)]


#         keep_unique_start_df = df_correction(keep_unique_start_df)
#         retain_clusters_df = pd.concat([retain_clusters_df,keep_unique_start_df])
#         keep_unique_start_df['pc-value'] = round(top_count/total,2) * 100 
        pc_val_return  = round(top_count/total,2) * 100
#         print(top_TSS,pc_val_return)
#         retain_clusters_df = pd.concat([retain_clusters_df,keep_unique_start_df])
#     retain_clusters_df.to_csv(output_file +f'/tmp/PC_ANALYSIS.{strand_map[strand]}.6.bed',sep = '\t',index = None)
        if top_count/total > min_abundance:
#             keep_unique_start_df = current_df[current_df['start']== top_TSS]
            keep_unique_start_df = current_df[(current_df['start'] < top_TSS+pad) & (current_df['start'] > top_TSS-pad)]
            keep_unique_start_df = df_correction(keep_unique_start_df)
            retain_clusters_df = pd.concat([retain_clusters_df,keep_unique_start_df])
#             keep_unique_start_df = current_df[(current_df['start'] < top_TSS+pad) & (current_df['start'] > top_TSS-pad)]
            # keep_unique_start_df = df_correction(keep_unique_start_df)
#             retain_clusters_df = pd.concat([retain_clusters_df,keep_unique_start_df])
    return retain_clusters_df 
    
def df_correction(current_df_TSS):
    highest_TSS = current_df_TSS['start'].value_counts().head(1).index[0]

    current_df_TSS['val-correct'] = highest_TSS - current_df_TSS['start']
#     print(current_df_TSS['start'].value_counts())
    current_df_TSS['start'] = current_df_TSS['start'] + current_df_TSS['val-correct']
    current_df_TSS['thickStart'] = current_df_TSS['thickStart'] + current_df_TSS['val-correct']
#     print(current_df_TSS['val-correct'].value_counts())

#     blocksize_correction_lst = []   
#     for correct,blocksizes in zip(current_df_TSS['val-correct'],current_df_TSS['blockSizes']):
#         old_blockSize_ind1 = int(blocksizes.split(',')[0])
#         correction = str(old_blockSize_ind1-correct)
# #         print(correct)
#         blocksize_correction_lst.append(blocksizes.replace(str(old_blockSize_ind1),correction))
#     current_df_TSS['blockSizes'] = blocksize_correction_lst
    current_df_TSS  = current_df_TSS.drop(['val-correct'],axis = 1)
    return current_df_TSS 
######### TESTING ########## 

def get_individual_clusters(full_df):
    """Takes in the full dataframe and collapses it
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
    requiredGrp.add_argument("-c",'--cigar_filter', required=False, help="cigar filter on TSS strand, default 12",default=12.5,type = float)
    requiredGrp.add_argument("-m",'--minimum_cluster_seqs', required=False, help="Retains clusters that have a minimum m (default 20)",default=20,type = float)
    requiredGrp.add_argument("-mi",'--minimum_isoform', required=False, help="Retains clusters that have a minimum m (default 5)",default=5,type = float)
    requiredGrp.add_argument("-pc",'--post_clustering_filter', required=False, help="Within a TU, if the top TSS does not have at least (pc) of the reads, this cluster is ignored (default .25)",default=.25,type = float)
    requiredGrp.add_argument("-TSS_pad",'--padding_val_TSS', required=False, help="TSS peak is defined by most abundant TSS in TU, (+/- value) are used to include surrounding reads and correct for more reads(default 12)",default=12,type = float)
    requiredGrp.add_argument("-TES_pad",'--padding_val_TES', required=False, help="TES peak is defined by most abundant TES in TU, (+/- value) are used to include surrounding reads and correct for more reads(default 25)",default=25,type = float)
    requiredGrp.add_argument("-b",'--bedgraphs', required=False, help="Output bedgraph coverage data to tmp file (default False)",default=False,type = bool)
    requiredGrp.add_argument("-s",'--noise_filter_TSS', required=False, help="Filter out TSS sites that have a count < s",default=3,type = float)
    requiredGrp.add_argument("-p",'--noise_filter_CPAS', required=False, help="Filter out CPAS sites that have a count < c",default=1,type = float)
    requiredGrp.add_argument("-f",'--isoform_grouping', required=False, help="Value to cluster isoform blocksizes",default=10,type = float)
    requiredGrp.add_argument("-t",'--nanopolish_tag', required=False, help="How to apply nanopolish tag filter, P - PASS, N - NO_PASS, A - All reads (N + P)",default = 'P',type = str)
    args = vars(ap.parse_args())
    parameters_df = pd.DataFrame.from_dict(args, orient='index')

    input_file = args['input_file']
    output_file = args['output_directory']
    TES_group_val = int(args['TES_grouping'])
    TSS_group_val = int(args['TSS_grouping'])
    nanopolish_path = args['nanopolish']
    cigar_filter_val = int(args['cigar_filter'])
    min_clust = int(args['minimum_cluster_seqs'])
    pc_filter = float(args['post_clustering_filter'])
    minimum_isoform = float(args['minimum_isoform'])
    padding_value_TSS = float(args['padding_val_TSS'])
    padding_value_TES = float(args['padding_val_TES'])
    filter_noise_TSS_val = float(args['noise_filter_TSS'])
    filter_noise_CPAS_val = float(args['noise_filter_CPAS'])
    blocksize_clustering = float(args['isoform_grouping'])
    np_tag = str(args['nanopolish_tag'])
    ##The bool() function is not recommended as a type converter. All it does is convert empty strings to False and non-empty strings to True. This is usually not what is desired.
    bedgraph_bool = args['bedgraphs']
    
    create_tmp_files(input_file,output_file +'/tmp')
    cigar_file_path = output_file +'/tmp/' + 'seq-cigar-orient.tmp'
    
    parameters_df.to_csv(output_file +'/' +'NAGATA-parameters.tsv',header = None,sep = '\t')
    polyA = pd.read_csv(nanopolish_path,sep = '\t')
    raw_bed12 = pd.read_csv(output_file +'/tmp/' + 'raw-alignment.2.bed',sep = '\t',header = None,usecols = [3,5],names = ['readname','strand'])
    final_counts = pd.DataFrame()
    full_cluster_df = pd.DataFrame()
    for strand in ['+','-']:
        os.system(f'echo Currently processing {strand_map[strand]} strand...')
        sys.stdout = open(output_file +'/' +f'Filtering-counts.{strand_map[strand]}.tsv', 'w')
        #Step3: Initial filtering using nanopolish polyA estimates and duplicated ‘readname’ information
        strand_sequences = nanopolish_filter.return_nanopolish_dedup(raw_bed12,polyA,strand,np_tag)
        print(f'Total input reads:\t',raw_bed12.shape[0])
        #Step4: Additional filtering for low abundance TSS (>3) and TES (>1)
        bed12_nano_dedup = read_in_bed(output_file +'/tmp/' + 'raw-alignment.2.bed',strand,strand_sequences,filter_noise_TSS_val,filter_noise_CPAS_val)
        bed12_nano_dedup.to_csv(output_file +f'/tmp/nanopolish.low-abundance_TSS_CPAS.{strand_map[strand]}.4.bed',sep = '\t',index = None)
        #Step5: Parses seg-cigar-orient.tmp to extract soft-clipping and filter by value specified by cigar_filter_val
        filter_cigar = TSS_soft_clip_filter.filter_sequences(cigar_file_path,cigar_filter_val,strand)
        filter_cigar.to_csv(output_file +f'/tmp/cigar-parsing.{strand_map[strand]}.bed',sep = '\t',index = None)
        cigar_retain_lst = bed12_nano_dedup[bed12_nano_dedup['seq-name'].isin(list(filter_cigar['sequence']))]
#         print(cigar_retain_lst)
        cigar_retain_lst.to_csv(output_file +f'/tmp/cigar-nanopolish.{strand_map[strand]}.5.bed',sep = '\t',index = None)
        print(f'Reads with soft-clipping (value:{cigar_filter_val}) filter:\t',bed12_nano_dedup.shape[0]-cigar_retain_lst.shape[0])
        #Step6: Takes TSS values of filtered BED file and returns grouped transcriptional unit
        TSS_df = TSS_grouping.TSS_grouping(cigar_retain_lst,TSS_group_val)
        TSS_df.to_csv(output_file +f'/tmp/TSS_clustering.{strand_map[strand]}.6.bed',sep = '\t',index = None)
        #Step7: Takes in initial TUs defined by TSSs (Step6), and defines TU further by using TES information in a similar way. 
        final_TU_cluster_df,TES_counts = TES_grouping.TES_grouping(TSS_df,TES_group_val,padding_value_TES)
        #final_TU_cluster_df.to_csv(output_file +'/tmp/initial-grouping.{strand_map[strand]}.txt',sep = '\t',index = None)
        final_TU_cluster_df.to_csv(output_file +f'/tmp/CPAS_clustering.{strand_map[strand]}.7.bed',sep = '\t',index = None)
        #Step8: Filter transcriptional units by top TSS abundance + padding correction 
        final_filtered_clusters =  filter_clusters(final_TU_cluster_df,pc_filter,padding_value_TSS)
        final_filtered_clusters.to_csv(output_file +f'/tmp/TU-top_abund_filter.padding_correction.{strand_map[strand]}.8.bed',sep = '\t',index = None)
        print(f'Reads available for transcriptome construction:\t',final_filtered_clusters.shape[0],'\n')
        print(f'Number of transcription units identified:\t',len(set(final_filtered_clusters['Final_clusters'])))
        
        #Step9: Retain TUs that have transcript count at least (-m)
        filtered_min_clusters = [TU for TU,count in Counter(final_filtered_clusters['Final_clusters']).items() if count > min_clust]
        final_filtered_clusters = final_filtered_clusters[final_filtered_clusters['Final_clusters'].isin(filtered_min_clusters)]
        final_filtered_clusters.to_csv(output_file +f'/tmp/TU-count_filtering.{strand_map[strand]}.9.bed',sep = '\t',index = None)
        print(f'Number of transcription units passing (value:{min_clust}) filter:\t',len(set(final_filtered_clusters['Final_clusters'])),'\n')
        
        #Step10: Within a TU, transcripts are initially separated by blockCounts and then clustered according to the sum of blocksizes 
        final_filtered_clusters['blockSize-sums'] = iso_grouping.get_blocksize_length(final_filtered_clusters['blockSizes'])  
        final_filtered_clusters_iso = iso_grouping.run_isoform(final_filtered_clusters,blocksize_clustering)
        final_filtered_clusters_iso.to_csv(output_file +f'/tmp/isoform-level-clustering.{strand_map[strand]}.10.bed',sep = '\t',index = None)
        print(f'Number of isoforms identified:\t',len(set(final_filtered_clusters_iso['TU.#-iso.#'])))
        #Step11: Retain transcripts that have a count of at least (-mi)
        final_filtered_clusters_iso = final_filtered_clusters_iso[final_filtered_clusters_iso['TU.#-iso.#-count']> minimum_isoform]  
        #Step12: Collapse final bed file
        final_filtered_clusters_iso_bed = get_individual_clusters(final_filtered_clusters_iso)
        #Step13: Format final bed output
        final_filtered_clusters_iso_bed = formatting_output(final_filtered_clusters_iso_bed)

        print(f'Number of isoforms passing (value:{int(minimum_isoform)}) filter:\t',final_filtered_clusters_iso_bed.shape[0])
        #Step14: Print out final collapse bed file and convert to GFF3 file
        final_filtered_clusters_iso_bed.to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.bed',sep ='\t',index = None,header = None)
        gff3_file = bed_convert.run_BED2GFF3(final_filtered_clusters_iso_bed)
        #Step15: Print out final PRE-collapse bed and GFF3 to output directory
        gff3_file.to_csv(output_file+'/Final_cluster.NAGATA.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
        final_filtered_clusters_iso.to_csv(output_file+'/Final_cluster.precollapsed.' + strand_map[strand]+'.tsv',sep ='\t',index = None)
        #Output each 
        bedgraph_outs = output_file+'/bedgraphs-'+strand_map[strand]
        os.makedirs(bedgraph_outs,exist_ok=True)
        sys.stdout.close()
        import pysam
        if bedgraph_bool == True:
            for TU_isoforms in set(final_filtered_clusters_iso['TU.#-iso.#']):
                os.system(f'echo {TU_isoforms}')
                current_seqs = list(final_filtered_clusters_iso[final_filtered_clusters_iso['TU.#-iso.#'] == TU_isoforms]['seq-name'])
                create_bedgraph(input_file,output_file +'/tmp',bedgraph_outs+'/'+TU_isoforms +'.'+strand_map[strand]+ '.bedgraph',TU_isoforms,current_seqs)

