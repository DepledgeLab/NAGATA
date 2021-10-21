import pandas as pd

def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))
    
def return_nanopolish_dedup(df_bed_2,nano,strand):
    """Takes in Bed12 and nanopolish files filters 
    """
    df_bed_2_strand = df_bed_2[df_bed_2['strand'] == strand]
#     print(df_bed_2_strand['strand'].value_counts())
    combined_df = df_bed_2_strand.merge(nano,left_on = 'readname',right_on = 'readname')
    retain_reads = combined_df[combined_df['qc_tag'] == 'PASS']
    initial_filter = retain_reads[~retain_reads['readname'].duplicated()]
    return list(initial_filter['readname'])