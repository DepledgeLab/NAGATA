import subprocess
import os
import pandas as pd
import modules.help_functions as help_functions

def get_blocksize_length(lst):
    """Uses bed file blockSize and returns overall length
    """
    return_lst = []
    for i in lst:
        if len(i.split(',')) >1:
            ints = sum([int(item) for item in i.split(',')])
            return_lst.append(ints)
        else:
            return_lst.append(int(i))
    return return_lst
    
def return_isoform_(TU_df,grouping_value):
    """Within TUs, uses blockCount to initially seperate isoforms, and then groups isoforms by sum of blockSizes
    """
    count = 0
    isoform_df = pd.DataFrame()
    for i in set(TU_df['blockCount']):
        current_TU = TU_df[TU_df['blockCount']==i]
        current_TU = current_TU.groupby("blockSize-sums").filter(lambda x: len(x) > 2)
        isoform_groups = list(help_functions.grouper(sorted(current_TU['blockSize-sums'].value_counts().index),grouping_value))
        for group,lengths in enumerate(isoform_groups):
            count += 1 
            current_TU_isos = current_TU[current_TU['blockSize-sums'].isin(lengths)]
            current_TU_isos['isoform'] = str(count)
            isoform_df = pd.concat([isoform_df,current_TU_isos])
    return isoform_df

def run_isoform(df,grouping_value):
    final_df = pd.DataFrame()
    all_TUs = set(df['Final_clusters'])
    for TU_names in all_TUs:
        individual_TUs = df[df['Final_clusters'] == TU_names]
        id_isoforms = return_isoform_(individual_TUs,grouping_value)
        final_df = pd.concat([final_df,id_isoforms])
    TU_map = dict(enumerate(set(final_df['Final_clusters']),start =1))
    TU_map = {v:k for k,v in TU_map.items()}
    final_df['TU.#-iso.#'] = 'TU.' + final_df['Final_clusters'].map(TU_map).astype(str) +'-iso.' + final_df['isoform']
    final_df['TU.#-iso.#-count'] = final_df['TU.#-iso.#'].map(dict(final_df['TU.#-iso.#'].value_counts()))
    final_df = final_df.drop(['Trans-unit','end-clusters','Final_clusters','isoform'],axis =1)
    return final_df

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input dataset and returns answer...')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="input file location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']
    final_df = run_isoform(input_file,filter_value)
#     print(final_df)
    final_df.to_csv(output_file,sep = '\t',index = None)