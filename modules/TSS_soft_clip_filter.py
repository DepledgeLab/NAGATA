import pandas as pd
import re

import warnings
warnings.filterwarnings("ignore")

def get_softclipping_TES(s):
    res = [re.findall(r'(\w+?)(\d+)', s)[0] ][0]
    if res[0] == 'S':
        return int(res[-1][::-1])
    else:
        return 0
        
def get_cigar_vals(df,strand):
    """Returns sequence-cigar pair
    """
    if strand == '+':
        return_cigars_vals = [int(cig.split('S')[0]) if cig.split('S')[0].isnumeric() else 0 for cig in df['cigar']]
    else:
        df['cigar'] = df['cigar'].apply(lambda x: x[::-1])
        return_cigars_vals = [get_softclipping_TES(cig) for cig in df['cigar']]
    df['soft_clip_values'] = return_cigars_vals
#     print(df.head())
    return df
    
def filter_sequences(seq_cigar_file_path,filt_val,strands):
    """ Uses sequence-cigar_string pair extracted from .sam file to extra length of 5' soft clipping
    and filter out values greater than filt_val.
    Returns: Filtered list of sequences
    """

    df_seq_cig = pd.read_csv(seq_cigar_file_path,sep ='\t',names = ['sequence','strandness','cigar'])
    new_strand = {16:'-',0:'+'}
    df_seq_cig['strandness'] = df_seq_cig['strandness'].map(new_strand)
    
    final_df = pd.DataFrame()
    df_seq_cig_current = df_seq_cig[df_seq_cig['strandness'] == strands]

    df_cigar_calc = get_cigar_vals(df_seq_cig_current,strands)
    df_cigar_calc_filt = df_cigar_calc[df_cigar_calc['soft_clip_values'] < filt_val]
#     final_df = pd.concat([final_df,df_cigar_calc_filt])
    return df_cigar_calc_filt,df_cigar_calc
    
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input dataset and returns answer...')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="input file location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    requiredGrp.add_argument("-v",'--filter_value', required=True, help="Value to filter soft-clipping TSS")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']
    filter_value = float(args['filter_value'])
    final_df = filter_sequences(input_file,filter_value,0)
#     print(final_df)
    final_df.to_csv(output_file,sep = '\t',index = None)
    