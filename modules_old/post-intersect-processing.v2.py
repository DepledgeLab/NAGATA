import pandas as pd
import subprocess

def get_blocksize_length(lst):
    """Uses bed file blockSize and returns overall length
    """
    return_lst = []
    for i in lst:
#         print(i,type(i))
        if len(str(i).split(',')) >1:
            ints = sum([int(item) for item in i.split(',') if len(item) != 0])
            return_lst.append(ints)
        else:
            return_lst.append(int(i))
    return return_lst    
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input intersect file (between NAGATA-isoforms and "some-known"-isoform) and returns best overlapping row')
    #bedtools intersect -wo -split -a Final_cluster.fwd.bed -b ../../Adenovirus/Annotation/Ad5_v10.0.forward.bed
    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-a",'--annotation', required=True, help="known annotation file (BED12 file)")
    requiredGrp.add_argument("-n",'--nagata_output', required=True, help="NAGATA output (BED12 file)")
    requiredGrp.add_argument("-o",'--output_directory', required=True, help="output file directory")
    requiredGrp.add_argument("-g",'--NAGATA_GFF3', required=False, help="GFF3 file")   
    requiredGrp.add_argument("-v",'--overlap_threshold', required=False, help="Overlap value in which to consider two transcripts similar or different",default=20,type = float)
    args = vars(ap.parse_args())
    annotation_file = args['annotation']
    nagata_file = args['nagata_output']
    output_file = args['output_directory']
    gff_file = args['NAGATA_GFF3']
    overlap_parameter = args['overlap_threshold']
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
    df[28] = abs(df[26]-df[25]) ## Absolute value of Annotation-isoform - NAGATA-isoform
    full_df = pd.DataFrame()
    for iso in set(df[3]):
    ### For every NAGATA-isoform, return the intersection that has the minimum difference between bedtool.intersect command and blocksize
        current_df = df[df[3]==iso]
        current_df = current_df[current_df[27]==current_df[27].min()] 
        current_df = current_df[current_df[28]==current_df[28].min()]
        full_df = pd.concat([full_df,current_df])

        
    full_df = full_df[full_df[28]< overlap_parameter]
    full_df = full_df[full_df[27]< overlap_parameter]
    
    final_df = pd.DataFrame()
    for known_trans in set(full_df[15]):
        current_df = full_df[full_df[15]==known_trans]
        current_df = current_df[current_df[27]==current_df[27].min()] 
        current_df = current_df[current_df[28]==current_df[28].min()]
        final_df = pd.concat([final_df,current_df])
    overlap_dict = dict(zip(final_df[3],final_df[15]))
    nagata_overlaps = nagata_df[~nagata_df[3].isin(list(overlap_dict.keys()))]
    annotation_overlaps = annotation_df[~annotation_df[3].isin(list(overlap_dict.values()))]
#     print('ANNOTATIONS',annotation_overlaps,annotation_overlaps.shape[0])
#     print('nagata_overlaps',nagata_overlaps,nagata_overlaps.shape[0])
    nagata_overlaps.to_csv(output_file + '/NAGATA-specific.bed',sep ='\t',index = None)
    annotation_overlaps.to_csv(output_file + '/Annotation-specific.bed',sep ='\t',index = None)
    output_string = 'Overlapping\t' + str(len(overlap_dict)) + '\n' + 'NAGATA-specific\t' + str(nagata_overlaps.shape[0]) + '\n' + 'Annotation-specific\t' + str(annotation_overlaps.shape[0]) 
    print(output_string)
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
