import pandas as pd
from collections import Counter, defaultdict
import numpy as np
# import modules.help_functions as help_functions
import warnings

warnings.filterwarnings("ignore")

def return_gene(TU_level_df):
    TU_level_df = TU_level_df[[0,1,2,4,5,12]]
    TU_level_df.columns = ['seq-id','feature-start','feature-end','score','strand','atr']
    TU_level_df['source'] = 'NAGATA'
    TU_level_df['score'] = '.'
    TU_level_df['feature-type'] = 'gene'
    TU_level_df['phase'] = '.'
    head_TU = TU_level_df.head(1)[['seq-id','source','feature-type','feature-start','feature-end','score','strand','phase','atr']]
    head_TU['atr'] = 'ID=gene.'+ head_TU['atr']
    return head_TU
    
def return_mRNA(TU_level_df):
    TU_level_df = TU_level_df[[0,1,2,4,5,3,12]]
    TU_level_df.columns = ['seq-id','feature-start','feature-end','score','strand','atr','gene']
    TU_level_df['source'] = 'NAGATA'
    TU_level_df['feature-type'] = 'mRNA'
    TU_level_df['phase'] = '.'
#     head_TU = TU_level_df[['seq-id','source','feature-type','feature-start','feature-end','score','strand','phase','atr']]
    TU_level_df['atr'] = 'ID=mRNA.'+ TU_level_df['atr'] + ';Parent=gene.' +TU_level_df['gene'] +';Name=' + TU_level_df['atr']
    TU_level_df = TU_level_df[['seq-id','source','feature-type','feature-start','feature-end','score','strand','phase','atr']]
    return TU_level_df

def return_exon(TU_level_df):
    TU_top_row = list(TU_level_df[[0,5,12]].iloc[0,:])
    seq_id,strand,gene = TU_top_row[0],TU_top_row[1],TU_top_row[2]
    sorted_exon_coord = get_sorted_exon(TU_level_df)
    parent_dict = get_parent_info(TU_level_df)
    final_df = pd.DataFrame()
    for ex_number,start_stops in sorted_exon_coord:
        exon_name = 'exon'+str(ex_number)
        current_exon_coordinate = '-'.join(list(map(str,start_stops)))
        parent_values = parent_dict[current_exon_coordinate]
        full_parent_name = ','.join(['mRNA.'+parent for parent in parent_values])
        atr_value = 'ID={}.{};Parent={}'.format(exon_name,gene,full_parent_name)
        current_exon_return = pd.DataFrame([seq_id,'NAGATA','exon',start_stops[0],start_stops[1],len(parent_values),strand,'.',atr_value]).T
        final_df = pd.concat([final_df,current_exon_return])
    final_df.columns = ['seq-id', 'source', 'feature-type', 'feature-start', 'feature-end', 'score', 'strand', 'phase', 'atr']
    final_df['score']  = '.'
    return final_df

def get_sorted_exon(test):
    full_list = []
    name = ''
    for _,rows in test.iterrows():
        exon_starts = list(map(lambda x:x+rows[1],map(int,rows[11].split(','))))
        exon_ends = [start + sizes for start, sizes in zip(exon_starts, map(int,rows[10].split(',')))]
        name = rows[12] 
        for start_ends in list(zip(exon_starts,exon_ends)):
            if start_ends not in full_list:
                full_list.append(start_ends)
    sorted_exons = list(enumerate(sorted(full_list),start =1))
    return sorted_exons

def get_parent_info(test):
    parent_dict = defaultdict(list)
    for _,rows in test.iterrows():
        exon_starts = list(map(lambda x:x+rows[1],map(int,rows[11].split(','))))
        exon_ends = [start + sizes for start, sizes in zip(exon_starts, map(int,rows[10].split(',')))]
        exon_starts = list(map(str,exon_starts))
        exon_ends = list(map(str,exon_ends))
        final_ = ['-'.join(exon_full) for exon_full in list(zip(exon_starts,exon_ends))]
        for exon_data in final_:
            parent_dict[exon_data].append(rows[3])
    return parent_dict
    
    
def run_BED2GFF3(df):
    df.columns = list(range(12))
    df[12] = [i.split('-')[0] for i in df[3]]
    final_GFF = pd.DataFrame()
    for TU in set(df[12]):
        current_TU = df[df[12] ==TU]
        get_gene = return_gene(current_TU)
        get_mrna =return_mRNA(current_TU)
        get_exon = return_exon(current_TU)
        final_GFF = pd.concat([final_GFF,get_gene,get_mrna,get_exon])
    return final_GFF

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input dataset')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="Bed12 file input location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']


    #### Handled in main function
    df = pd.read_csv(input_file,sep = '\t',header = None)
    df[12] = [i.split('-')[0] for i in df[3]]
    final_df = run_BED2GFF3(df)
    ####
    final_df.to_csv(output_file,sep = '\t',index = None,header = None)
