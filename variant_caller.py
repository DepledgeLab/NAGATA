import pandas as pd
import operator
import argparse


def parsing_indels(nucleotide_counts:'dict') -> 'DataFrame':
    """Returns the top3 indels sorted by frequency
    """
    NA_vals = ['NaN','NaN','NaN']
    potential_indels = {key:value for key,value in nucleotide_counts.items() if len(key) > 1 }
    nucleotide_counts['N'] = sum(map(int,potential_indels.values()))
    indel_top3 = [':'.join(map(str,i)) for i in sorted(potential_indels.items(), key=lambda x: int(x[1]), reverse=True)[:3]]
    raw_nucleotides_df = pd.DataFrame.from_dict(nucleotide_counts,orient ='index').T[['A','C','G','T','N']]
    indel_top3 += NA_vals
    top3_df = pd.DataFrame(indel_top3[:3], index = ['N#1','N#2','N#3']).T
    return pd.concat([raw_nucleotides_df,top3_df],axis =1)

def getting_variant_info(raw_dict:'dict',depth_value,ref) -> 'list':
    """Grabs the indel information and does minor formatting for downstream processes
    """
    max_variant = max(raw_dict.items(), key=operator.itemgetter(1))[0]
    max_variant_count = raw_dict[max_variant]
    variant_freq = round(max_variant_count/depth_value,3)
    if max_variant_count == 0:
        max_variant = '-'
        max_variant_count = '-'
        variant_freq = '-'
    if type(variant_freq) == float and variant_freq > .5:
        consensus = max_variant
    else:
        consensus = ref
    variant_lst = [max_variant,variant_freq,consensus]
#     print(variant_lst)
    return variant_lst

def initial_parsing(path)-> 'DF':
    """Takes in the raw bamreadcounts, seperates nucleotides and indels, returns an initial DF
    """
    full_df = pd.DataFrame()
    with open(path) as input_data:
        for lines in input_data:
            initial_split = lines.strip().split('\t')
            first_4_dict = dict(zip(['Chromosome','Position','Reference','Depth'],initial_split[:4]))
            indels = { i[0]:i[1] for i in list(map(lambda x: x.split(':'),initial_split[4:])) if i[0] != '='}
            df_start = pd.DataFrame.from_dict(first_4_dict,orient ='index').T
            parsed_events_df = parsing_indels(indels)
            combined_df = pd.concat([df_start,parsed_events_df],axis =1)
            full_df = pd.concat([full_df,combined_df])
    full_df['Depth'] = full_df[['A','C','G','T','N']].astype(int).sum(axis =1)
    return full_df

def indel_formatting(indel_counts:'list') -> 'dict':
    """Takes in the list of indels, and returns a dictionary with only deletions and subsitutions
    """
    current_dict = { j.split(':')[0]:int(j.split(':')[1]) for j in indel_counts if len(j.split(':')) > 1 }
    return current_dict

def return_variant_lst(full_df):
    """Formatting to determine the variant, variant_freq, and consensus columns
    """
    llst = []
    for k,v in full_df.iloc[:,2:].iterrows():
        reference = v['Reference']
        depth = v['Depth']
        nucleotide_dict = {k:int(v) for k,v in dict(v[2:-4]).items()}
        indels_dict = indel_formatting(v[-3:].values)
        nucleotide_dict.update(indels_dict)
        del nucleotide_dict[reference]
        current_lst = getting_variant_info(nucleotide_dict,depth,reference)
        llst.append(current_lst)
    return llst

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description = 'Takes in input dataset and returns answer...')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="input file location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']
    
    df = initial_parsing(input_file)
    adjusted_variant_lst = return_variant_lst(df)
    df[['Variant','Variant_Freq','Consensus']] = adjusted_variant_lst
    ordered_df = df[['Chromosome','Position','Reference','Variant','Variant_Freq','Consensus','Depth','A','C','G','T','N','N#1','N#2','N#3']]
    ordered_df.to_csv(output_file,sep = '\t',index = None)
    
    
    
    