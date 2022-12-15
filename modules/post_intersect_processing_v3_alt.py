import pandas as pd
import subprocess
import os 
import numpy as np
from collections import Counter
def get_blocksize_length(lst):
    """Uses bed file blockSize and returns overall length
    """
    return_lst = []
    for i in lst:
#         print(i,type(i))
        if len(str(i).split(',')) >1:
            ints = sum([int(item) for item in i.split(',') if len(item) != 0][1:])
            return_lst.append(ints)
        else:
            return_lst.append(int(i))
    return return_lst
def compare_ind_blocks(nagata_blockSizes,annotate_blockSizes,some_threshold,blockCount):
    pass_threshold = []
    # all([True if abs(nagata_block-anno_block) < some_threshold else False  for nagata_block,anno_block in zip(x[1:],y[1:])])
    if blockCount > 1:
        nagata_blocksize_list = list(map(int,str(nagata_blockSizes).strip(',').split(',')))
        
        annotate_blocksize_list = list(map(int,str(annotate_blockSizes).strip(',').split(',')))
        for nagata_block,anno_block in zip(nagata_blocksize_list,annotate_blocksize_list):
            if abs(nagata_block-anno_block) < some_threshold:
                pass_threshold.append(True)
            else:
                pass_threshold.append(False)
        return all(pass_threshold)
    else: 
        if abs(int(str(nagata_blockSizes).strip(','))-int(str(annotate_blockSizes).strip(','))) < some_threshold:
            return True
        else:
            return False

def run_overlap_scoring(output_file,nagata_file,annotation_file,overlap_parameter=50,ind_overlap_parameter=20):
    os.makedirs(output_file,exist_ok=True)
    fout_bedtools = open(output_file +"/Raw_BEDTOOLS_intersect_results.tmp",'wb')
    bedtools_p1 = f"bedtools intersect -wo -split -a {nagata_file} -b {annotation_file}".split(" ")
    bedtools_p2 = subprocess.run(bedtools_p1,stdout=fout_bedtools)
    df = pd.read_csv(output_file +"/Raw_BEDTOOLS_intersect_results.tmp",sep = '\t',header = None)
    nagata_df = pd.read_csv(nagata_file,sep = '\t',header = None)
    annotation_df = pd.read_csv(annotation_file,sep ='\t',header =None)
    df = df[df[9]==df[21]]  ## NAGATA-isoform and annotation-isoform need to have same number of exons
    df[25] = get_blocksize_length(df[10]) ## NAGATA-isoform blocksize sums
    df[26] = get_blocksize_length(df[22]) ## Annotation-isoform blocksize sums
    df[27] = abs(df[25]-df[24]) ## Absolute value of NAGATA-isoform - bedtools.intersect command
    df[28] = abs(df[26]-df[25]) ## Absolute value of Annotation-isoform - NAGATA-isoform)
    df.to_csv(output_file + '/Before-filtering.overlaps.bed',sep ='\t',index = None)

    full_df = df[df[28]< overlap_parameter]
    full_df = full_df[full_df[27]< overlap_parameter]  

    keep_rows = []
    for _,rows in df.iterrows():
        current_indv = compare_ind_blocks(rows[10],rows[22],ind_overlap_parameter,rows[9])
        keep_rows.append(current_indv)
    df[29] = keep_rows
    df = df[df[29] == True]
    
    keep_rows = []
    for _,rows in df.iterrows():
        current_indv = compare_ind_blocks(rows[11],rows[23],ind_overlap_parameter,rows[9])
        keep_rows.append(current_indv)
    df[29] = keep_rows
    df = df[df[29] == True]
    
    greater_than_1_transcript = [k for k,v in Counter(df[15]).items() if v>1]
    found_nagata_trans = []
    max_abundance_df = pd.DataFrame()
    for known_tran in greater_than_1_transcript:
        current_df = df[df[15] == known_tran]
#         print(found_nagata_trans, known_tran)
#         print(current_df)
        current_df = current_df[~current_df[3].isin(found_nagata_trans)]
        current_df = current_df.sort_values(by = 4,ascending = False).head(1)
        if current_df.shape[0] == 0:
            current_df = df[df[15] == known_tran]
            current_df = current_df.sort_values(by = 28).head(1)
            max_abundance_df = pd.concat([max_abundance_df,current_df])
        else: 
            max_abundance_df = pd.concat([max_abundance_df,current_df])
            found_nagata_trans.append(current_df[3].to_list()[0])
        # print(current_df[3].to_list()[0])
    df = df[~df[15].isin(greater_than_1_transcript)]
    final_df = pd.concat([df,max_abundance_df])
    
    final_df[30] = get_blocksize_length(final_df[11])
    final_df[31] = get_blocksize_length(final_df[23])
    final_df[32] = abs(final_df[30] - final_df[31])
#     final_df.to_csv(output_file + '/After.unique.known.bed',sep ='\t',index = None)

    greater_than_1_nagata = [k for k,v in Counter(final_df[3]).items() if v>1]
    filter_repeat_df = pd.DataFrame()
    for i in greater_than_1_nagata:
        current_df = final_df[final_df[3]==i]
        minimum_blockstarts = current_df[current_df[32]==current_df[32].min()]
        filter_repeat_df = pd.concat([filter_repeat_df,minimum_blockstarts])
    unique_NAGATA = final_df[~final_df[3].isin(greater_than_1_nagata)]
    output_overlaps_df = pd.concat([unique_NAGATA,filter_repeat_df])
    output_overlaps_df = output_overlaps_df.drop([24,25,26,28,29,27,30,31,32],axis = 1)
    
    
    
    overlap_dict = dict(zip(output_overlaps_df[3],output_overlaps_df[15]))
    nagata_overlaps = nagata_df[~nagata_df[3].isin(list(output_overlaps_df[3]))]

    annotation_overlaps = annotation_df[~annotation_df[3].isin(list(output_overlaps_df[15]))]
    nagata_overlaps[3] = nagata_overlaps[3] + '--' + nagata_overlaps[4].astype(str)
    nagata_overlaps.to_csv(output_file + '/NAGATA-specific.bed',sep ='\t',index = None,header = None)
    annotation_overlaps.to_csv(output_file + '/Annotation-specific.bed',sep ='\t',index = None,header = None)
    output_string = 'Validated\t' + str(output_overlaps_df.shape[0]) + '\n' + 'Not_annotated\t' + str(nagata_overlaps.shape[0]) + '\n' + 'Not_detected\t' + str(annotation_overlaps.shape[0]) 
#     print(output_string)
    text_file = open(output_file +"/overlap-performance.txt", "w")
    n = text_file.write(output_string)
    text_file.close()
    
    multiple_overlaps = [k for k,v in Counter(output_overlaps_df[3]).items() if v>1]

    tmp_df_single = output_overlaps_df[~output_overlaps_df[3].isin(multiple_overlaps)]
    tmp_df_multiple = output_overlaps_df[output_overlaps_df[3].isin(multiple_overlaps)]


    overlap_dict = {}
    for i in multiple_overlaps:
        current_df = tmp_df_multiple[tmp_df_multiple[3] == i]
        overlap_dict[i] = ','.join(current_df[15].tolist())
    tmp_df_multiple[15] = tmp_df_multiple[3].map(overlap_dict)
    final_output_overlap = pd.concat([tmp_df_multiple,tmp_df_single])
    final_output_overlap = final_output_overlap[~final_output_overlap[3].duplicated()]
    #final_output_overlap.sort_values(by=[1,2,9]).to_csv(output_file + '/NAGATA-Annotation.overlaps.bed',sep ='\t',index = None,header = None)
    return final_output_overlap, nagata_overlaps, annotation_overlaps
#     output_overlaps_df.to_csv(output_file + '/NAGATA-Annotation.overlaps.bed',sep ='\t',index = None,header = None)
    
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input intersect file (between NAGATA-isoforms and "some-known"-isoform) and returns best overlapping row')
    #bedtools intersect -wo -split -a Final_cluster.fwd.bed -b ../../Adenovirus/Annotation/Ad5_v10.0.forward.bed
    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-a",'--annotation', required=True, help="known annotation file (BED12 file)")
    requiredGrp.add_argument("-n",'--nagata_output', required=True, help="NAGATA output (BED12 file)")
    requiredGrp.add_argument("-o",'--output_directory', required=True, help="output file directory")
    requiredGrp.add_argument("-g",'--NAGATA_GFF3', required=False, help="GFF3 file")   
    requiredGrp.add_argument("-b",'--sum_blocksize_threshold', required=False, help="Overlap value in which to consider two transcripts similar or different (sum of blockSize)",default=50,type = float)
    requiredGrp.add_argument("-e",'--individual_exon_threshold', required=False, help="Overlap value in which to consider two transcripts similar or different (comparing individual exons)",default=20,type = float)
    args = vars(ap.parse_args())
    annotation_file = args['annotation']
    nagata_file = args['nagata_output']
    output_file = args['output_directory']
    gff_file = args['NAGATA_GFF3']
    overlap_parameter = args['sum_blocksize_threshold']
    ind_overlap_parameter = args['individual_exon_threshold']
    
    
    os.makedirs(output_file,exist_ok=True)
    fout_bedtools = open(output_file +"/Raw_BEDTOOLS_intersect_results.tmp",'wb')
    bedtools_p1 = f"bedtools intersect -wo -split -a {nagata_file} -b {annotation_file}".split(" ")
    bedtools_p2 = subprocess.run(bedtools_p1,stdout=fout_bedtools)
    df = pd.read_csv(output_file +"/Raw_BEDTOOLS_intersect_results.tmp",sep = '\t',header = None)
    nagata_df = pd.read_csv(nagata_file,sep = '\t',header = None)
    annotation_df = pd.read_csv(annotation_file,sep ='\t',header =None)

    df = df[df[9]==df[21]]  ## NAGATA-isoform and annotation-isoform need to have same number of exons
    df[25] = get_blocksize_length(df[10]) ## NAGATA-isoform blocksize sums
    df[26] = get_blocksize_length(df[22]) ## Annotation-isoform blocksize sums
    df[27] = abs(df[25]-df[24]) ## Absolute value of NAGATA-isoform - bedtools.intersect command
    df[28] = abs(df[26]-df[25]) ## Absolute value of Annotation-isoform - NAGATA-isoform)
    df.to_csv(output_file + '/Before-filtering.overlaps.bed',sep ='\t',index = None)

    full_df = df[df[28]< overlap_parameter]
    full_df = full_df[full_df[27]< overlap_parameter]  

    final_df = pd.DataFrame()
    for known_trans in set(full_df[15]):
        current_df = full_df[full_df[15]==known_trans]
        current_df = current_df[current_df[27]==current_df[27].min()] 
        current_df = current_df[current_df[28]==current_df[28].min()]
        final_df = pd.concat([final_df,current_df])

    keep_rows = []
    for _,rows in final_df.iterrows():
        current_indv = compare_ind_blocks(rows[10],rows[22],ind_overlap_parameter,rows[9])
        keep_rows.append(current_indv)
    final_df[29] = keep_rows
    final_df = final_df[final_df[29] == True]
    final_df = final_df.drop([24,25,26,29],axis = 1)
    
    overlap_dict = dict(zip(final_df[3],final_df[15]))
    nagata_overlaps = nagata_df[~nagata_df[3].isin(list(final_df[3]))]

    annotation_overlaps = annotation_df[~annotation_df[3].isin(list(final_df[15]))]

    nagata_overlaps.to_csv(output_file + '/NAGATA-specific.bed',sep ='\t',index = None)
    annotation_overlaps.to_csv(output_file + '/Annotation-specific.bed',sep ='\t',index = None)
    output_string = 'Overlapping\t' + str(final_df.shape[0]) + '\n' + 'NAGATA-specific\t' + str(nagata_overlaps.shape[0]) + '\n' + 'Annotation-specific\t' + str(annotation_overlaps.shape[0]) 
#     print(output_string)
    text_file = open(output_file +"/overlap-performance.txt", "w")
    n = text_file.write(output_string)
    text_file.close()
    final_df.to_csv(output_file + '/NAGATA-Annotation.overlaps.bed',sep ='\t',index = None)
    if gff_file != None:
        print(overlap_dict) 
        with open(gff_file, 'r') as file :
            filedata = file.read()
        for NAGATA, annotation in overlap_dict.items():
            filedata = filedata.replace(NAGATA, annotation)
        with open('NEW-TEST.gff3', 'w') as file:
            file.write(filedata)
