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
import modules.post_intersect_processing_v4.1 as post_scoring
import warnings
import shutil
## Ignore warning dealing from "A value is trying to be set on a copy of a slice from a DataFrame..." 
warnings.filterwarnings("ignore")


strand_map = {'+':'fwd','-':'rev'}
#colnames_fwd = ['chrom','start','end','name','score','strand','thickstart','thickend','itemRGB','blockcount','blocksizes','blockstarts']
reference_files_scoring = {'+':'r1','-':'r2'}

parsing_dictionary ={'input_file':'-i','output_directory':'-o', 'polya':'-p','soft_clip_filter':'-s','nanopolish_tag':'-nt','CPAS_noise_filter':'-c',
'TSS_noise_filter':'-t', 'CPAS_clustering':'-cg','TSS_clustering':'-tg','isoform_clustering':'-iso','min_transcript_abundance':'-m','TSS_abundance_per_TU':'-a',
'blocksize_noise_filter':'-b','reference_bed_f':'-r1','reference_bed_r':'-r2','max_TSS_per_CPAS':'-pm','padding_TSS':'-pt','padding_CPAS':'-pc','-tmp':'temporary_files'}


if __name__ == '__main__':
    import argparse

    ap = argparse.ArgumentParser(description = 'Takes in input dataset')
#help
    
    optional = ap._action_groups.pop()
    requiredGrp = ap.add_argument_group('required arguments')
    optional_basic = ap.add_argument_group('Optional (basic) arguments')
    optional_advanced = ap.add_argument_group('Optional (advanced) arguments')
    requiredGrp.add_argument("-i", metavar= 'input_file', required=True, help="BAM file input location")
    requiredGrp.add_argument("-o",metavar='output_directory', required=True, help="output file location",default = '.',type = str)
    optional_basic.add_argument("-p",metavar= 'polya', required=False, help="nanopolish location output location, (default None)",default = None)
    optional_basic.add_argument("-nt",metavar='nanopolish_tag', required=False, help="How to apply nanopolish tag filter, P - PASS, N - NO_PASS, A - All reads (N + P), (default P) [Not available in version 2.2]",default = 'P',type = str)
    optional_basic.add_argument("-s",metavar='soft_clip_filter', required=False, help="cigar filter on TSS strand, (default 3)",default=3,type = float)
    optional_basic.add_argument("-c",metavar='CPAS_noise_filter', required=False, help="Filter out CPAS with count < c, (default 20)",default=20,type = int)
    optional_basic.add_argument("-t",metavar='TSS_noise_filter', required=False, help="Filter out TSS with count < t, (default 4)",default=4,type = int)
    optional_advanced.add_argument("-cg",metavar='CPAS_clustering', help="Grouping value for clustering (CPAS), (default 50)",default = 50.0,type = int)
    optional_advanced.add_argument("-tg",metavar='TSS_clustering', help="Grouping value for clustering (TSS), (default 20)",default = 20.0,type = int)
    optional_advanced.add_argument("-iso",metavar='isoform_clustering', help="Grouping value for blockSize and blockStarts, (default 50)",default = 50.0,type = int)
    optional_basic.add_argument("-m",metavar='min_transcript_abundance', required=False, help="minimum transcript abundance, (default 3 )",default = 3,type = int)
    optional_advanced.add_argument("-a",metavar='TSS_abundance_per_TU', help="TSS abundance per transcriptional unit, (default 0.1)",default = 0.1,type = float)
    optional_advanced.add_argument("-b",metavar='blocksize_noise_filter', help="Prior to isoform deconvolution - filter out low abundant blocksize sums to prevent incorrect daisy chainging, (default 3)",default = 3,type = float)
    optional_basic.add_argument("-r1",metavar='reference_bed_f', required=False, help="BED file containing existing annotation (forward strand)",default = None,type = str)
    optional_basic.add_argument("-r2",metavar='reference_bed_r', required=False, help="BED file containing existing annotation (reverse strand)",default = None,type = str)
    optional_basic.add_argument("-d",metavar='strand', help="specify if interested in only forward (+) or reverse (-) strand, (default both)",default = ['+','-'],type = list)
    optional_basic.add_argument("-tmp",metavar='temporary_files', help="output intermediate files to a tmp directory (default False)",default = False,type = str)

    
    args = vars(ap.parse_args())
    parameters_df = pd.DataFrame.from_dict(args, orient='index')

    parameters_df[1] = parameters_df.index
    parameters_df[2] = parameters_df[1].map(parsing_dictionary)
    
    input_file = args['i']
    output_file = args['o']
    nanopolish_path = args['p']
    cigar_filter_val = int(args['s'])
    CPAS_noise_filter = int(args['c'])
    TSS_noise_filter = int(args['t'])
    CPAS_grouping_val = int(args['cg'])
    TSS_grouping_val = int(args['tg'])
    isoform_grouping_val = int(args['iso'])
    min_transcript_abundance = int(args['m'])
    TSS_abundance_per_TU = float(args['a'])
    np_tag = str(args['nt'])
    blocksizesum_noise_filter = int(args['b'])
    strands = list(args['d'])
    known_beds_f = args['r1']
    known_beds_r = args['r2']
    output_tmp_files = args["tmp"]
    ##The bool() function is not recommended as a type converter. All it does is convert empty strings to False and non-empty strings to True. This is usually not what is desired.
    print(input_file.split('.')[-1])
    is_bam = True
    if input_file.split('.')[-1] == 'bed':
        is_bam = False
    if is_bam:
        helpers.create_tmp_files(input_file,output_file +'/tmp')
        cigar_file_path = output_file +'/tmp/' + 'seq-cigar-orient.tmp'
    

    parameters_df[[2,1,0]].to_csv(output_file +'/' +'parameters.tsv',header = ['flag','parameter_name','value'],sep = '\t',index = None)
    if nanopolish_path != None:
        polyA = pd.read_csv(nanopolish_path,sep = '\t',usecols =[0,9] )
        passed_polyA = polyA[polyA['qc_tag'] =='PASS']
        passed_polyA.to_csv(output_file +f'/tmp/Passed_NANOPOLISH.bed',sep = '\t',index = None)
    final_counts = pd.DataFrame()

    filtering_counts = ""

    for strand in strands:
        if strand == '-':
            colnames = ['chrom','end','start','name','score','strand','thickend','thickstart','itemRGB','blockcount','blocksizes','blockstarts']
            if is_bam:
                df = pd.read_csv(output_file +'/tmp/' + '1.raw-alignment.bed',sep = '\t',header = None,names = colnames)
            else: 
                df = pd.read_csv(input_file,sep = '\t',header = None,names = colnames)
#                 print('HERE -')
        else:
            colnames = ['chrom','start','end','name','score','strand','thickstart','thickend','itemRGB','blockcount','blocksizes','blockstarts']
            if is_bam:
                df = pd.read_csv(output_file +'/tmp/' + '1.raw-alignment.bed',sep = '\t',header = None,names = colnames)
            else: 
                df = pd.read_csv(input_file,sep = '\t',header = None,names = colnames)
#                 print('HERE +')
        print(f'Currently processing {strand_map[strand]} strand...\n')

        ##1 Filter for relevant strand
        df = df[df['strand'] ==strand]
        # print("raw-count",df.shape)
        filtering_counts += f'{strand_map[strand]} strand...\nTotal input reads:\t {df.shape[0]}\n'
        ##2 Filter for sequences that contain a PASS for nanopolish
        if nanopolish_path != None:
            
            df = df.merge(passed_polyA,left_on ="name", right_on ="readname",how = "left").drop(["readname"],axis =1)
            
            df["qc_tag"] = df["qc_tag"].fillna("fail-nanopolish")
            df.to_csv(output_file +f'/tmp/2.filter_nanopolish.{strand_map[strand]}.bed',sep = '\t',index = None)
        if nanopolish_path == None:
            print(nanopolish_path)
            df["qc_tag"] = "PASS"

        filtering_counts += f'Reads with readable polyA tail:\t {df[df["qc_tag"]!="fail-nanopolish"].shape[0]}\n'
        ##3 Cigar string filter
        if is_bam:
            filter_cigar,cigar_df = TSS_soft_clip_filter.filter_sequences(cigar_file_path,cigar_filter_val,strand)
            # print(filter_cigar.shape)
            cigar_df.to_csv(output_file +f'/tmp/parsed.cigar.{strand_map[strand]}.tsv',sep = '\t',index = None)
            df = df.merge(filter_cigar[["sequence","soft_clip_values"]],left_on ="name", right_on ="sequence",how = "left").drop(["sequence"],axis =1)
            df["soft_clip_values"] = df["soft_clip_values"].fillna("fail-soft-clip")

        print("total size",df.shape)

        df_fail3 = df[(df["qc_tag"]== "fail-nanopolish") & (df["soft_clip_values"]!= "fail-soft-clip")]
        df_fail5 = df[(df["qc_tag"]!= "fail-nanopolish") & (df["soft_clip_values"]== "fail-soft-clip")]
        df_failboth = df[(df["qc_tag"]== "fail-nanopolish") & (df["soft_clip_values"]== "fail-soft-clip")]
        df_passboth = df[(df["qc_tag"]!= "fail-nanopolish") & (df["soft_clip_values"]!= "fail-soft-clip")]

        df.to_csv(output_file +f'/tmp/3.filter_cigar.{strand_map[strand]}.bed',sep = '\t',index = None)
        filtering_counts += f"Reads with acceptable 5\' alignment:\t {df.shape[0]}\n\n"

        ##4 Filter low abundant CPAS followed by CPAS daisy chaining
        df_3_define = pd.concat([df_passboth,df_fail5])


        CPAS_pass = CPAS_processing.filter_CPAS_noise(df_3_define['end'],CPAS_noise_filter)
        CPAS_groups = helpers.identify_daisy_chain_groups(sorted(map(int,CPAS_pass.keys())),CPAS_grouping_val)

        new_CPAS_groups = []
        for i in CPAS_groups:
            df_identified_groups = df_3_define[(df_3_define['end'].isin(i)) ]
            most_abundant_cpas = df_identified_groups['end'].value_counts().head(1).index[0]
            tmp_df = df_3_define.copy()
            tmp_df['absolute_end'] = abs(tmp_df['end']-most_abundant_cpas)

            tmp_df = tmp_df[tmp_df['absolute_end'] <= CPAS_grouping_val/2]
            new_CPAS_groups.append(sorted(set(tmp_df['end'])))
        swapped_ends = helpers.swap_key_vals(dict(enumerate(new_CPAS_groups,1)))
        df_3_define['Transcriptional-unit'] = df_3_define['end'].astype(int).map(swapped_ends)
        df_3_define = df_3_define[df_3_define['Transcriptional-unit'].notna()]
        df_3_define['Transcriptional-unit'] = 'TU.' + df_3_define['Transcriptional-unit'].astype(int).astype(str)
        most_abund_CPAS_in_groups = dict(df_3_define.groupby('Transcriptional-unit')['end'].agg(lambda handle_multiple_modes: pd.Series.mode(handle_multiple_modes)[0]))
        df_3_define['most_abund_CPAS'] = df_3_define['Transcriptional-unit'].map(most_abund_CPAS_in_groups)       
        filtering_counts += f"Reads remaining after CPAS noise filter:\t {df.shape[0]} \n"
        filtering_counts += f"Number of transcription units (CPAS) identified:\t {len(set(df_3_define['Transcriptional-unit']))}\n\n"
        df_3_define['TU-count']= df_3_define['Transcriptional-unit'].map(df_3_define['Transcriptional-unit'].value_counts())
        df_3_define.to_csv(output_file +f'/tmp/4.CPAS-grouping.{strand_map[strand]}.bed',sep = '\t',index = None)


        df_5_define = pd.concat([df_passboth,df_fail3])
        TSS_pass = CPAS_processing.filter_CPAS_noise(df_5_define['start'],TSS_noise_filter)
        TSS_groups = helpers.identify_daisy_chain_groups(sorted(map(int,TSS_pass.keys())),TSS_grouping_val)

        new_TSS_groups = []
        for i in TSS_groups:
            df_identified_groups = df_5_define[df_5_define['start'].isin(i)]
            most_abundant_tss = df_identified_groups['start'].value_counts().head(1).index[0]
            tmp_df = df_5_define.copy()
            tmp_df['absolute_start'] = abs(df_5_define['start']-most_abundant_tss)

            tmp_df = tmp_df[tmp_df['absolute_start'] <= TSS_grouping_val/2]
            new_TSS_groups.append(sorted(set(tmp_df['start'])))        
        swapped_starts = helpers.swap_key_vals(dict(enumerate(new_TSS_groups,1)))
        df_5_define['TSS-unit'] = df_5_define['start'].astype(int).map(swapped_starts)
        df_5_define = df_5_define[df_5_define['TSS-unit'].notna()]
        df_5_define['TSS-unit'] = 'TSS.' + df_5_define['TSS-unit'].astype(int).astype(str)
        most_abund_TSS_in_groups = dict(df_5_define.groupby('TSS-unit')['start'].agg(lambda handle_multiple_modes: pd.Series.mode(handle_multiple_modes)[0]))
        df_5_define['most_abund_TSS'] = df_5_define['TSS-unit'].map(most_abund_TSS_in_groups)
        df_5_define['TSS-count']= df_5_define['TSS-unit'].map(df_5_define['TSS-unit'].value_counts()) 
        df_5_define.to_csv(output_file +f'/tmp/5.TSS-grouping.{strand_map[strand]}.bed',sep = '\t',index = None)   

        ### TESTING new functionality
        ## (1) Get most abundant TSS (2) Look into which intersect TSS from raw BED file 

        df_5_define_abundant_TSS = df_5_define[["chrom","most_abund_TSS"]].drop_duplicates()
        df_5_define_abundant_TSS["start"] = df_5_define_abundant_TSS["most_abund_TSS"] - 25
        df_5_define_abundant_TSS["end"] = df_5_define_abundant_TSS["most_abund_TSS"] + 25
        df_5_define_abundant_TSS = df_5_define_abundant_TSS.drop(["most_abund_TSS"],axis = 1)
        df_5_define_abundant_TSS.to_csv(output_file +f'/tmp/most_abundant.TSSs.{strand_map[strand]}.bed',sep = '\t',index = None,header = None)
        most_abundant_TSS_path = output_file +f'/tmp/most_abundant.TSSs.{strand_map[strand]}.bed'
        df_most_abundant_TSS_path = open(output_file +'/tmp/' + '1.raw-alignment.TSS.bed','wb')


        df_starts = df[['chrom','start']].drop_duplicates()
        df_starts['end'] = df_starts['start']
        df_starts.to_csv(output_file +f'/tmp/df_raw_tss.{strand_map[strand]}.bed',sep="\t",header = None,index = None)
        raw_file_path = output_file +f'/tmp/df_raw_tss.{strand_map[strand]}.bed'
        
        bedtools_p1 = f"bedtools intersect -wb -a {raw_file_path} -b {most_abundant_TSS_path}".split(" ")
        bedtools_p2 = subprocess.run(bedtools_p1,stdout=df_most_abundant_TSS_path)

        # print(df_3_define)
        df_3_define_abundant_cpas = df_3_define[["chrom","most_abund_CPAS"]].drop_duplicates()
        df_3_define_abundant_cpas["start"] = df_3_define_abundant_cpas["most_abund_CPAS"] - 25
        df_3_define_abundant_cpas["end"] = df_3_define_abundant_cpas["most_abund_CPAS"] + 25
        df_3_define_abundant_cpas = df_3_define_abundant_cpas.drop(["most_abund_CPAS"],axis = 1)
        most_abundant_cpas_path = output_file +f'/tmp/most_abundant.CPASs.{strand_map[strand]}.bed'
        df_3_define_abundant_cpas.to_csv(most_abundant_cpas_path,sep = '\t',index = None,header = None)
        
        df_most_abundant_CPAS_path = open(output_file +'/tmp/' + '1.raw-alignment.CPAS.bed','wb')

        df_end = df[['chrom','end']].drop_duplicates()

        df_end['start'] = df_end['end']
        raw_file_path_cpas = output_file +f'/tmp/df_raw_cpas.{strand_map[strand]}.bed'
        df_end.to_csv(raw_file_path_cpas,sep="\t",header = None,index = None)
        
        bedtools_p1 = f"bedtools intersect -wb -a {raw_file_path_cpas} -b {most_abundant_cpas_path}".split(" ")
        bedtools_p2 = subprocess.run(bedtools_p1,stdout=df_most_abundant_CPAS_path)


        overlaps_cpas = pd.read_csv(output_file +'/tmp/' + '1.raw-alignment.CPAS.bed',sep ="\t",usecols = [1,4],header = None,names = ["end","most_abund_CPAS"])
        overlaps_tss = pd.read_csv(output_file +'/tmp/' + '1.raw-alignment.TSS.bed',sep ="\t",usecols = [1,4],header = None,names = ["start","most_abund_TSS"])
        df = df.merge(overlaps_tss[["start"]],on = "start",how ="inner")
        df = df.merge(overlaps_cpas[["end"]],on = "end",how ="inner")

        df = df.merge(df_3_define[["end","Transcriptional-unit","most_abund_CPAS","TU-count"]].drop_duplicates(),on = "end",how = "inner")
        df = df.merge(df_5_define[["start","most_abund_TSS","TSS-unit","TSS-count"]].drop_duplicates(),on = "start",how = "inner")





        ###
        print(df[(df["TSS-unit"]=="TSS.9") & (df["Transcriptional-unit"]=="TU.6") &(df["blockcount"] ==1) ])
        df.to_csv(output_file +f'/tmp/Delete.Non-Corrected.TSS.CPAS.{strand_map[strand]}.bed',sep = '\t',index = None)  
        df['CPAS-diff'] = df['most_abund_CPAS'] - df['end']
        df['end'] = df['most_abund_CPAS']
        
        df['TSS-diff'] = df['most_abund_TSS'] - df['start']
        df['start'] = df['most_abund_TSS']
 

        df.to_csv(output_file +f'/tmp/Correct.TSS.CPAS.{strand_map[strand]}.bed',sep = '\t',index = None)  


        df = blocksize_processing.correct_last(df,strand) ###
        df = blocksize_processing.correct_first(df,strand)  ###

        df_full_correct = df


        df_full_correct['new-name'] = df_full_correct['Transcriptional-unit'] +'-'+ df_full_correct['TSS-unit'].astype(str)
        df_full_correct.to_csv(output_file +f'/tmp/6.columns-fully-corrected.{strand_map[strand]}.bed',sep = '\t',index = None)



        df_full_correct['blocksize-sum'] = blocksize_processing.get_blocksize_length(df_full_correct['blocksizes.new'])
        df_full_correct['new-name.ex'] = df_full_correct['new-name'] +'.'+ df_full_correct['blockcount'].astype(str)
        
        
        df_final = iso_deconv.iso_deconv(df_full_correct,isoform_grouping_val,blocksizesum_noise_filter)
        df_final.to_csv(output_file +f'/tmp/7.Isoform-deconvolution.{strand_map[strand]}.bed',sep = '\t',index = None)
        Final_df = post_pro.dataframe_editing(df_final)

        Final_df = Final_df.sort_values(by=['end','start','blocksizes'])
        

        most_common_df = pd.DataFrame()
        for i in set(Final_df['full-id']):
            current_df = Final_df[Final_df['full-id']==i]
            top_blocksize = current_df.value_counts('blocksizes.new').head(1).index[0]
            current_df = current_df[current_df['blocksizes.new'] == top_blocksize]
            top_blockstart = current_df.value_counts('blockstarts').head(1).index[0]
            current_df = current_df[current_df['blockstarts'] == top_blockstart].head(1)
            most_common_df = pd.concat([most_common_df,current_df])

        filtering_counts += f"Number of isoforms identified:\t {len(set(most_common_df['full-id']))}\n"

        most_common_df = most_common_df[most_common_df['TSS-abundance-per-TU'] >= TSS_abundance_per_TU]
        filtering_counts += f"Number of isoforms identified which pass TSS abundance per TU filter:\t {len(set(most_common_df['full-id']))}\n"
        most_common_df = most_common_df[most_common_df['score'] >= min_transcript_abundance]
        filtering_counts += f"Number of isoforms identified which pass minimum transcript count:\t {len(set(most_common_df['full-id']))}\n\n\n"
        
        most_common_df['name'] = most_common_df['name'] + '--' + most_common_df['TSS-abundance-per-TU'].astype(str)
        most_common_df = most_common_df.sort_values( by=['start','end','blockcount'])
        most_common_df[colnames].to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.bed',sep ='\t',index = None,header = None)
        gff3_file = bed_convert.run_BED2GFF3(most_common_df[colnames],None)
        gff3_file['feature-start'] = gff3_file['feature-start'] + 1 

        gff3_file.sort_values(by = ['feature-start','feature-end']).to_csv(output_file+'/Final_cluster.NAGATA.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
        

        known_bed = args[reference_files_scoring[strand]]
        if known_bed != None:
            print(known_bed)
            nagata_annot = most_common_df[colnames]
            known_df = pd.read_csv(known_bed,sep ='\t')
            nagata_file =output_file+'/Final_cluster.' + strand_map[strand]+'.bed'
            output_file_by_strand = output_file + '/overlap.' + strand_map[strand]
            final_output_overlap, nagata_specific, annotation_specific = post_scoring.run_overlap_scoring(output_file_by_strand,nagata_file,known_bed,50,20)
            NAGATA_known_mapping = dict(zip(final_output_overlap[3],final_output_overlap[15]))
            test_outputs = {k:k.replace('--',f'--{v}--')for k,v in NAGATA_known_mapping.items()}
            nagata_annot = nagata_annot.replace({'name':test_outputs})
            nagata_annot['name'] = nagata_annot['name'] + '--' + nagata_annot['score'].astype(str)
            nagata_annot.to_csv(output_file+'/Final_cluster.' + strand_map[strand]+'.bed',sep ='\t',index = None,header = None)

            final_output_overlap = final_output_overlap.replace({3:test_outputs})
            final_output_overlap[3] = final_output_overlap[3] + '--' + final_output_overlap[4].astype(str)
            final_output_overlap.sort_values(by=[1,2,9]).to_csv(output_file_by_strand + '/NAGATA-Annotation.overlaps.bed',sep ='\t',index = None,header = None)
            
            gff3_file_overlap = bed_convert.run_BED2GFF3(final_output_overlap[[0,1,2,3,4,5,6,7,8,9,10,11]],'4C33FF')
            gff3_file_nagata = bed_convert.run_BED2GFF3(nagata_specific,'FF3333')
            gff3_file_annotation = bed_convert.run_BED2GFF3(annotation_specific,'B0AEAE')
            
            final_overlap_gff = pd.concat([gff3_file_overlap,gff3_file_nagata])
            final_overlap_gff['feature-start'] = final_overlap_gff['feature-start'] + 1 
            final_overlap_gff.sort_values(by = ['feature-start','feature-end']).to_csv(output_file+'/Final_cluster.NAGATA.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
            
            gff3_file_overlap.sort_values(by = ['feature-start','feature-end']).to_csv(output_file_by_strand+'/NAGATA.annotation.overlap.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
            gff3_file_nagata.sort_values(by = ['feature-start','feature-end']).to_csv(output_file_by_strand+'/NAGATA.specific.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
            gff3_file_annotation.sort_values(by = ['feature-start','feature-end']).to_csv(output_file_by_strand+'/Annotation.specific.' + strand_map[strand]+'.gff3',sep = '\t',index = None,header = None)
            print(output_file_by_strand)

    if known_bed != None and len(strands) == 2:
        forward_overlap_score = pd.read_csv(output_file+'/overlap.' + 'fwd/' +'overlap-performance.txt',sep = '\t',names = ['class','forward'])
#         print(forward_overlap_score)
        reverse_overlap_score = pd.read_csv(output_file+'/overlap.' + 'rev/' +'overlap-performance.txt',sep = '\t',names = ['class','reverse'])

        merged_scores = forward_overlap_score.merge(reverse_overlap_score,left_on='class',right_on='class')
        merged_scores['total-overlap'] = merged_scores['forward'] + merged_scores['reverse']

        merged_scores.to_csv(output_file+'/merged_scoring.txt',sep = '\t',index = None)
    elif strands == ['+']:
        forward_overlap_score = pd.read_csv(output_file+'/overlap.' + 'fwd/' +'overlap-performance.txt',sep = '\t',names = ['class','forward'])
        print(forward_overlap_score)
    elif strands == ['-']:
        reverse_overlap_score = pd.read_csv(output_file+'/overlap.' + 'rev/' +'overlap-performance.txt',sep = '\t',names = ['class','reverse'])
        print(forward_overlap_score)
processing_file = open(f"{output_file}/Filtering-counts.txt","w")


processing_file.write(filtering_counts)

processing_file.close() #to change file access modes

if not output_tmp_files:
    shutil.rmtree(output_file +f'/tmp/', ignore_errors=True)

        
