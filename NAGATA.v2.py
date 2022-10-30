from collections import Counter
import os
import subprocess
import pandas as pd

import modules.CPAS_processing as CPAS_processing
import modules.help_scripts as helpers
import modules.TSS_processing as TSS_processing
import modules.blocksize_processing as blocksize_processing
import modules.isoform_deconvolution as iso_deconv
import modules.post_processing as post_pro
import modules.TSS_soft_clip_filter as TSS_soft_clip_filter
import modules.BED2GFF3 as bed_convert
import modules.post_intersect_processing_v3_alt as post_scoring
import warnings
## Ignore warning dealing from "A value is trying to be set on a copy of a slice from a DataFrame..." 
warnings.filterwarnings("ignore")


strand_map = {'+':'fwd','-':'rev'}
colnames_fwd = ['chrom','start','end','name','score','strand','thickstart','thickend','itemRGB','blockcount','blocksizes','blockstarts']
reference_files_scoring = {'+':'reference_bed_f','-':'reference_bed_r'}




if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser(description = 'Takes in input dataset')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="BAM file input location")
    requiredGrp.add_argument("-o",'--output_directory', required=True, help="output file location",default = '.',type = str)
    requiredGrp.add_argument("-n",'--nanopolish', required=False, help="nanopolish location output location (requires concatenated fwd and rev)",default = None)
    requiredGrp.add_argument("-c",'--cigar_filter', required=False, help="cigar filter on TSS strand, default 12",default=3,type = float)
    requiredGrp.add_argument("-t",'--nanopolish_tag', required=False, help="How to apply nanopolish tag filter, P - PASS, N - NO_PASS, A - All reads (N + P)",default = 'P',type = str)
    requiredGrp.add_argument("-p",'--noise_filter_CPAS', required=False, help="Filter out CPAS sites that have a count < c",default=20,type = int)
    requiredGrp.add_argument("-s",'--noise_filter_TSS', required=False, help="Filter out TSS sites that have a count < c",default=4,type = int)
    requiredGrp.add_argument("-g",'--CPAS_grouping', required=False, help="Grouping value for clustering (CPAS), default 50",default = 50.0,type = int)
    requiredGrp.add_argument("-d",'--TSS_grouping', required=False, help="Grouping value for clustering (TSS), default 20",default = 20.0,type = int)
    requiredGrp.add_argument("-f",'--isoform_grouping', required=False, help="Grouping value for grouping blocksizes, default 50",default = 50.0,type = int)
    requiredGrp.add_argument("-m",'--min_transcript_abundance', required=False, help="minimum transcript abundance",default = 3,type = int)
    requiredGrp.add_argument("-a",'--TSS_abundance_per_TU', required=False, help="TSS abundance per transcriptional unit (%)",default = 2,type = float)
    requiredGrp.add_argument("-b",'--noise_filter_blocksize', required=False, help="Prior to isoform deconvolution - filter out low abundant blocksize sums to prevent incorrect daisy chainging ",default = 3,type = float)
    requiredGrp.add_argument("-r1",'--reference_bed_f', required=False, help="reference beds for scoring",default = None,type = str)
    requiredGrp.add_argument("-r2",'--reference_bed_r', required=False, help="reference beds for scoring",default = None,type = str)
    args = vars(ap.parse_args())
    parameters_df = pd.DataFrame.from_dict(args, orient='index')
    #print(parameters_df)
    
    input_file = args['input_file']
    output_file = args['output_directory']
    nanopolish_path = args['nanopolish']
    cigar_filter_val = int(args['cigar_filter'])
    CPAS_noise_filter = int(args['noise_filter_CPAS'])
    TSS_noise_filter = int(args['noise_filter_TSS'])
    CPAS_grouping_val = int(args['CPAS_grouping'])
    TSS_grouping_val = int(args['TSS_grouping'])
    isoform_grouping_val = int(args['isoform_grouping'])
    min_transcript_abundance = int(args['min_transcript_abundance'])
    TSS_abundance_per_TU = float(args['TSS_abundance_per_TU'])
    np_tag = str(args['nanopolish_tag'])
    blocksizesum_noise_filter = int(args['noise_filter_blocksize'])
#     known_beds_f = args['reference_beds_f']
#     known_beds_r = args['reference_bed_r']
    ##The bool() function is not recommended as a type converter. All it does is convert empty strings to False and non-empty strings to True. This is usually not what is desired.
    print(input_file.split('.')[-1])
    is_bam = True
    if input_file.split('.')[-1] == 'bed':
        is_bam = False
    if is_bam:
        helpers.create_tmp_files(input_file,output_file +'/tmp')
        cigar_file_path = output_file +'/tmp/' + 'seq-cigar-orient.tmp'
    

    parameters_df.to_csv(output_file +'/' +'NAGATA-parameters.tsv',header = None,sep = '\t')
    if nanopolish_path != None:
        polyA = pd.read_csv(nanopolish_path,sep = '\t',usecols =[0,9] )
    
    final_counts = pd.DataFrame()
    #full_cluster_df = pd.DataFrame()
    
    
    #df = pd.read_csv(input_file,sep ='\t',header = None,names = colnames)
    for strand in ['+','-']:
        
        if strand == '-':
            colnames = ['chrom','start','end','name','score','strand','thickend','thickstart','itemRGB','blockcount','blocksizes','blockstarts']
            if is_bam:
                df = pd.read_csv(output_file +'/tmp/' + 'raw-alignment.2.bed',sep = '\t',header = None,names = colnames)
            else: 
                df = pd.read_csv(input_file,sep = '\t',header = None,names = colnames)
                print('HERE -')
        else:
            colnames = ['chrom','start','end','name','score','strand','thickstart','thickend','itemRGB','blockcount','blocksizes','blockstarts']
            if is_bam:
                df = pd.read_csv(output_file +'/tmp/' + 'raw-alignment.2.bed',sep = '\t',header = None,names = colnames_fwd)
            else: 
                df = pd.read_csv(input_file,sep = '\t',header = None,names = colnames)
                print('HERE +')
        print(f'Currently processing {strand_map[strand]} strand...\n')
        filtering_counts = ""
            
#         df = pd.read_csv(output_file +'/tmp/' + 'raw-alignment.2.bed',sep = '\t',header = None,names = colnames_fwd)
        ##1 Filter for relevant strand
        df = df[df['strand'] ==strand]
        filtering_counts += f'Total input reads:\t {df.shape[0]}\n'
        ##2 Filter for sequences that contain a PASS for nanopolish
        if nanopolish_path != None:
            passed_polyA = polyA[polyA['qc_tag'] =='PASS']
            passed_polyA.to_csv(output_file +f'/tmp/Passed_NANOPOLISH.{strand_map[strand]}.2.bed',sep = '\t',index = None)
            df = df[df['name'].isin(passed_polyA['readname'].to_list())]
        filtering_counts += f'Reads with readable polyA tail:\t {df.shape[0]}\n'
        ##3 Cigar string filter
        if is_bam:
            filter_cigar = TSS_soft_clip_filter.filter_sequences(cigar_file_path,cigar_filter_val,strand)
            filter_cigar.to_csv(output_file +f'/tmp/filter_cigar.{strand_map[strand]}.3.bed',sep = '\t',index = None)
            df = df[df['name'].isin(filter_cigar['sequence'].to_list())]
        filtering_counts += f"Reads with acceptable 5\' alignment:\t {df.shape[0]}\n\n"
        ##4 Filter low abundant CPAS followed by CPAS daisy chaining
        print(df.head())
        CPAS_pass = CPAS_processing.filter_CPAS_noise(df['end'],CPAS_noise_filter)
        CPAS_groups = helpers.identify_daisy_chain_groups(sorted(map(int,CPAS_pass.keys())),CPAS_grouping_val)
        swapped_ends = helpers.swap_key_vals(dict(enumerate(CPAS_groups,1)))
        df['Transcriptional-unit'] = df['end'].astype(int).map(swapped_ends)
        df = df[df['Transcriptional-unit'].notna()]
        df['Transcriptional-unit'] = 'TU.' + df['Transcriptional-unit'].astype(int).astype(str)
        most_abund_CPAS_in_groups = dict(df.groupby('Transcriptional-unit')['end'].agg(lambda handle_multiple_modes: pd.Series.mode(handle_multiple_modes)[0]))
        df['most_abund_CPAS'] = df['Transcriptional-unit'].map(most_abund_CPAS_in_groups)
        
        filtering_counts += f"Reads remaining after CPAS noise filter:\t {df.shape[0]} \n"
        filtering_counts += f"Number of transcription units (CPAS) identified:\t {len(set(df['Transcriptional-unit']))}\n\n"
        
        df.to_csv(output_file +f'/tmp/CPAS-grouping.{strand_map[strand]}.4.bed',sep = '\t',index = None)
        df['TU-count']= df['Transcriptional-unit'].map(df['Transcriptional-unit'].value_counts())
        df.to_csv(output_file +f'/tmp/CPAS-grouping.{strand_map[strand]}.4.bed',sep = '\t',index = None)
        
        ##5 Add columns giving unique TSS names and corresponding counts for each within each TU
        
        df['TSS.unique'] = df['Transcriptional-unit']+ '-TSS.' + df['start'].astype(str)
        filtering_counts += f"Number of unique TSS values:\t {len(set(df['TSS.unique']))}\n"
        df['TSS-count'] = df['TSS.unique'].map(df['TSS.unique'].value_counts())
        ## Filter unique TSSs by noise filter
        df = df[df['TSS-count']>TSS_noise_filter]
        filtering_counts += f"Number of unique TSS values after noise filtering:\t {len(set(df['TSS.unique']))}\n\n"

        ## Daisy chaining of TSSs within TUs
        TU_with_TSS_df = TSS_processing.parse_TSS_within_TU(df,TSS_grouping_val)
        TU_with_TSS_df.to_csv(output_file +f'/tmp/TSS-grouping.{strand_map[strand]}.5.bed',sep = '\t',index = None)

        ##6 BED file correction of start, end, and blocksizes columns # df[abs(df['CPAS-diff']) > 10]
        TU_with_TSS_df['CPAS-diff'] = TU_with_TSS_df['most_abund_CPAS'] - TU_with_TSS_df['end']
        TU_with_TSS_df = TU_with_TSS_df[abs(TU_with_TSS_df['CPAS-diff']) < 25 ]
        TU_with_TSS_df['end'] = TU_with_TSS_df['most_abund_CPAS']
        
        TU_with_TSS_df['TSS-diff'] = TU_with_TSS_df['most_abund_TSS_in_TU'] - TU_with_TSS_df['start']
        TU_with_TSS_df = TU_with_TSS_df[abs(TU_with_TSS_df['TSS-diff']) < 25 ]
        TU_with_TSS_df['start'] = TU_with_TSS_df['most_abund_TSS_in_TU']
        
        ### APPLY TU based cigar filter
        
        # TU_with_TSS_df = TU_with_TSS_df.merge(filter_cigar[['sequence','soft_clip_values']],left_on = 'name',right_on='sequence')
#         TU['TU_soft_clipping'] = 
        ####
        ## Correct 1st and last blocksize values using TSS-diff and CPAS-diff values respectively 
#         print(len(set(TU_with_TSS_df['blocksizes'])))
#         TU_with_TSS_df[TU_with_TSS_df['blocksizes']]
#         TU_with_TSS_df = TU_with_TSS_df[~TU_with_TSS_df['blocksizes'].duplicated()]
#         TU_with_TSS_df[colnames].to_csv(f'Ad5.nocorrection.{strand_map[strand]}.bed',sep = '\t',header = None,index = None)
#         print(TU_with_TSS_df.head(),TU_with_TSS_df.columns)
        TU_with_TSS_df = blocksize_processing.correct_last(TU_with_TSS_df) ###

        TU_with_TSS_df = blocksize_processing.correct_first(TU_with_TSS_df)  ###
        df_full_correct = TU_with_TSS_df
#         print(df_full_correct.columns)
#         print(len(set(df_full_correct['blocksizes.new'])))
#         colnames = ['chrom','start','end','name','score','strand','thickstart','thickend','itemRGB','blockcount','blocksizes.new','blockstarts']
#         df_full_correct = TU_with_TSS_df[~TU_with_TSS_df['blocksizes.new'].duplicated()]
        
        TU_with_TSS_df.to_csv(output_file +f'/tmp/columns-fully-corrected.{strand_map[strand]}.6.bed',sep = '\t',index = None)
#         df_full_correct = TU_with_TSS_df
        ##7 Isoform deconvolution
        df_full_correct['new-name'] = TU_with_TSS_df['Transcriptional-unit'] +'-TSS_group.'+TU_with_TSS_df['TSS-group'].astype(str)
        ## New unique name for reads

        
        ## Get sum of each blocksizes and abundance of each 
        df_full_correct['blocksize-sum'] = blocksize_processing.get_blocksize_length(df_full_correct['blocksizes.new'])

        df_final = iso_deconv.iso_deconv(df_full_correct,isoform_grouping_val,blocksizesum_noise_filter)

        # colnames = ['chrom','start','end','name','score','strand','thickstart','thickend','itemRGB','blockcount','blocksizes.new','blockstarts']
#         df_final = df_final[~df_final['blocksizes.new'].duplicated()]
#         df_final[colnames].to_csv(f'Ad5.nocorrection.{strand_map[strand]}.bed',sep = '\t',header = None,index = None)
#         print(len(set(df_final['blocksizes.new'])))
        
        ######### TU BASED FILTERING 
#         df_final =df_final.merge(filter_cigar[['sequence','soft_clip_values']],left_on = 'name',right_on='sequence')
#         median_vals = df_final.groupby('new-name').soft_clip_values.agg('median')

        
#         df_final = df_final.merge(median_vals,left_on='new-name',right_on='new-name')

        df_final.to_csv(output_file +f'/tmp/Isoform-deconvolution.{strand_map[strand]}.7.bed',sep = '\t',index = None)

        Final_df = post_pro.dataframe_editing(df_final)
        print(len(set(Final_df['blocksizes.new'])))
        Final_df = Final_df.sort_values(by=['end','start','blocksizes'])
        
        Final_df.to_csv(output_file+'/Final_cluster.precollapsed.' + strand_map[strand]+'.tsv',sep ='\t',index = None)

        most_common_df = pd.DataFrame()
        for i in set(Final_df['full-id']):
            current_df = Final_df[Final_df['full-id']==i]

            current_df = current_df.sort_values(by='splicing-count',ascending = False).head(1)
            most_common_df = pd.concat([most_common_df,current_df])
        most_common_df[colnames].to_csv(f'most-common-df.{strand_map[strand]}.7.bed',sep ='\t',index = None)
        filtering_counts += f"Number of isoforms identified:\t {len(set(most_common_df['full-id']))}\n"

        most_common_df = most_common_df[most_common_df['TSS-abundance-per-TU'] > TSS_abundance_per_TU]
        filtering_counts += f"Number of isoforms identified which pass TSS abundance per TU filter:\t {len(set(most_common_df['full-id']))}\n"
        most_common_df = most_common_df[most_common_df['score'] > min_transcript_abundance]
        filtering_counts += f"Number of isoforms identified which pass minimum transcript count:\t {len(set(most_common_df['full-id']))}\n\n"
        
        most_common_df['name'] = most_common_df['name'] + '--' + most_common_df['TSS-abundance-per-TU'].astype(str)
        most_common_df = most_common_df.sort_values( by=['start','end','blockcount'])
        most_common_df[colnames_fwd].to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.bed',sep ='\t',index = None,header = None)
        gff3_file = bed_convert.run_BED2GFF3(most_common_df[colnames_fwd])
        gff3_file.to_csv(output_file+'/Final_cluster.NAGATA.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
        
        processing_file = open(f"{output_file}/Filtering-counts.{strand_map[strand]}.txt","w")


        processing_file.write(filtering_counts)

        processing_file.close() #to change file access modes
        known_bed = args[reference_files_scoring[strand]]
        if known_bed != None:
            print(known_bed)
            nagata_annot = most_common_df[colnames_fwd]
            known_df = pd.read_csv(known_bed,sep ='\t')
#             print(known_df.head())
            nagata_file =output_file+'/Final_cluster.' + strand_map[strand]+'.bed'
            output_file_by_strand = output_file + '/overlap.' + strand_map[strand]
            post_scoring.run_overlap_scoring(output_file_by_strand,nagata_file,known_bed,50,20)
    if known_bed != None:
        forward_overlap_score = pd.read_csv(output_file+'/overlap.' + 'fwd/' +'overlap-performance.txt',sep = '\t',names = ['class','forward'])
        reverse_overlap_score = pd.read_csv(output_file+'/overlap.' + 'rev/' +'overlap-performance.txt',sep = '\t',names = ['class','reverse'])
    #     print(forward_overlap_score)
    #     print(reverse_overlap_score)
        merged_scores = forward_overlap_score.merge(reverse_overlap_score,left_on='class',right_on='class')
        merged_scores['total-overlap'] = merged_scores['forward'] + merged_scores['reverse']
        merged_scores.to_csv(output_file+'/merged_scoring.txt',sep = '\t',index = None)
        print(merged_scores)

        
        
        