import pandas as pd
from collections import Counter
import numpy as np
import modules.help_functions as help_functions


    
def TSS_grouping(df,grouping_val):
    """Takes in dataframe returns each TSS with its grouped transcriptional unit
    """
    uniq_starts = sorted(df['start'].unique())
    results_dct = dict(enumerate(help_functions.grouper(uniq_starts,grouping_val)))
    swapped_starts = help_functions.swap_key_vals(results_dct)
    df['Trans-unit'] = df['start'].map(swapped_starts)
    return df
    
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input dataset')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="Bed12 file input location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    requiredGrp.add_argument("-g",'--grouping_val', required=True, help="Grouping value for clustering")
    requiredGrp.add_argument("-s",'--strand', required=True, help="Strand of bed file")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']
    grouping_val = int(args['grouping_val'])
    strand = args['strand']
    
    if strand == '+':
        names_list = ['chrom','start','end','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    else: 
        names_list = ['chrom','end','start','seq-name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    #### Handled in main function
    df = pd.read_csv(input_file,sep = '\t',header = None,names = names_list)
    df = df[df['strand'] == strand]
    df = df.drop(['score','strand','thickStart','thickEnd','itemRgb'],axis = 1)
    df = df.sort_values(by=['start','end']).reset_index(drop = True)
    ####
    df = TSS_grouping(df,grouping_val)
    df.to_csv(output_file,sep = '\t',index = None)

    
    
    
    
    
    
    
    