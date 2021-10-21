import pandas as pd
from collections import Counter
import numpy as np
import modules.help_functions as help_functions
import warnings
warnings.filterwarnings("ignore")
def get_top_occurence(g):
    return g['end'].value_counts().idxmax()  

def TES_ends(df,grouping_val):
    """Takes in dataframe returns each TES with its grouped transcriptional unit
    """
    includ_end_df = pd.DataFrame()
    for start_clust in set(df['Trans-unit']):
        current_clust = df[df['Trans-unit'] ==start_clust]
        end_values = [exon_end for exon_end,counts in Counter(current_clust['end']).items() if counts > 2]
        end_values = sorted(end_values)
        grouped_ends = dict(enumerate(help_functions.grouper(end_values,50)))
        current_clust = current_clust[current_clust['end'].isin(end_values)]
        current_clust['end-clusters'] = current_clust['end'].map(help_functions.swap_key_vals(grouped_ends))
        current_clust = current_clust[current_clust['end-clusters'].notna()]
        includ_end_df = pd.concat([includ_end_df,current_clust])
    includ_end_df['end-clusters'] = includ_end_df['end-clusters'].astype(int)
    includ_end_df['Final_clusters'] = 'TU-' +includ_end_df['Trans-unit'].astype(str) + '-' + includ_end_df['end-clusters'].astype(str)
    return includ_end_df

def return_TU_counts(grouped_TES_df,lst_of_clusters):
    topdf = grouped_TES_df.groupby('Final_clusters').apply(get_top_occurence)
#     lst_of_clusters = sorted(set(grouped_TES_df['Final_clusters']))
    TU_counts = pd.DataFrame()
    for TU in lst_of_clusters:
        top_end = topdf[[TU == i for i in topdf.index]][0] 
        current_TU = grouped_TES_df[grouped_TES_df['Final_clusters']==TU]
        count_df = Counter(current_TU['end'])
        current_TU['counts'] = current_TU['end'].map(count_df)
        end_counts_df = current_TU[current_TU['end'].duplicated() == False][['Final_clusters','end','counts']]
        end_counts_df = end_counts_df.sort_values(by = 'counts',ascending = False).reset_index(drop = True)
        end_counts_df.at[0, 'chosen'] = '*'
        TU_counts = pd.concat([TU_counts,end_counts_df])
    return TU_counts

def TES_grouping(grouped_TSS_df,grouping_val):
    clustered_TES = TES_ends(grouped_TSS_df,grouping_val)
    lst_of_clusters = sorted(set(clustered_TES['Final_clusters']))
    TES_counts = return_TU_counts(clustered_TES,lst_of_clusters)
    final_TU_cluster_df = pd.DataFrame()
    for i in lst_of_clusters:
        current_TU = clustered_TES[clustered_TES['Final_clusters'] == i]
        chosen_end = TES_counts[TES_counts['Final_clusters'] == i]['end'][0] 
        current_TU = current_TU[current_TU['end'] == chosen_end]
        final_TU_cluster_df = pd.concat([final_TU_cluster_df,current_TU])
    return final_TU_cluster_df,TES_counts
    
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input dataset')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="Bed12 file input location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    requiredGrp.add_argument("-g",'--grouping_val', required=True, help="Grouping value for clustering")
#     requiredGrp.add_argument("-s",'--strand', required=True, help="Strand of bed file")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']
    grouping_val = int(args['grouping_val'])

    #### Handled in main function
    df = pd.read_csv(input_file,sep = '\t')

    ####
    df,counts_df = TES_grouping(df,grouping_val)
    df.to_csv(output_file,sep = '\t',index = None)
    counts_df.to_csv('TES.count.txt',sep = '\t',index = None)
