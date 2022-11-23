import modules.help_scripts as helpers
import pandas as pd
def parse_TSS_within_TU(df:'DataFrame',TSS_daisy_chain_cal:int,padding_TSS:int)-> 'DataFrame':
    """Takes in dataframe with defined CPAS based TUs and returns 
    """
    unique_TUs = sorted(set(df['Transcriptional-unit']))
    TU_with_TSS_df = pd.DataFrame()
    for tu in unique_TUs:
        current_df = df[df['Transcriptional-unit'] == tu]
        unique_starts = sorted(set(current_df['start']))
        daisy_chain_TSS = helpers.identify_daisy_chain_groups(unique_starts,TSS_daisy_chain_cal)
#         print(current_df['start'].value_counts().head(50))
        
        new_TSS_groups = []
        for i in daisy_chain_TSS:
            df_identified_groups = current_df[current_df['start'].isin(i)]
#             print(i)
            most_abundant_tss = df_identified_groups['start'].value_counts().head(1).index[0]
#             print(most_abundant_tss)
            tmp_df = current_df.copy()
            tmp_df['absolute_start'] = abs(df['start']-most_abundant_tss)

            tmp_df = tmp_df[tmp_df['absolute_start'] < padding_TSS/2]
            new_TSS_groups.append(sorted(set(tmp_df['start'])))
#         print(daisy_chain_TSS)
#         print(new_TSS_groups,'\n\n\n')
        daisy_chain_TSS_id = helpers.swap_key_vals(dict(enumerate(daisy_chain_TSS,1)))
        current_df['TSS-group'] = current_df['start'].map(daisy_chain_TSS_id)
        current_df['TSS-group-count'] =current_df['TSS-group'].map(current_df['TSS-group'].value_counts())
        #print(dict(current_df.groupby('TSS-group')['start'].agg(lambda x: pd.Series.mode(x).iat[0])))
        most_abund_TSS_in_group = dict(current_df.groupby('TSS-group')['start'].agg(lambda x: pd.Series.mode(x).iat[0]))
        current_df['most_abund_TSS_in_TU'] = current_df['TSS-group'].map(most_abund_TSS_in_group)
        TU_with_TSS_df = pd.concat([TU_with_TSS_df,current_df])
    return TU_with_TSS_df