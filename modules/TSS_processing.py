import modules.help_scripts as helpers
import pandas as pd
def parse_TSS_within_TU(df:'DataFrame',TSS_daisy_chain_cal:int)-> 'DataFrame':
	"""Takes in dataframe with defined CPAS based TUs and returns 
	"""
	unique_TUs = sorted(set(df['Transcriptional-unit']))
	TU_with_TSS_df = pd.DataFrame()
	for tu in unique_TUs:
		current_df = df[df['Transcriptional-unit'] == tu]
		unique_starts = sorted(set(current_df['start']))
		daisy_chain_TSS = helpers.identify_daisy_chain_groups(unique_starts,TSS_daisy_chain_cal)
		daisy_chain_TSS_id = helpers.swap_key_vals(dict(enumerate(daisy_chain_TSS,1)))
		current_df['TSS-group'] = current_df['start'].map(daisy_chain_TSS_id)
		current_df['TSS-group-count'] =current_df['TSS-group'].map(current_df['TSS-group'].value_counts())
		#print(dict(current_df.groupby('TSS-group')['start'].agg(lambda x: pd.Series.mode(x).iat[0])))
		most_abund_TSS_in_group = dict(current_df.groupby('TSS-group')['start'].agg(lambda x: pd.Series.mode(x).iat[0]))
		current_df['most_abund_TSS_in_TU'] = current_df['TSS-group'].map(most_abund_TSS_in_group)
		TU_with_TSS_df = pd.concat([TU_with_TSS_df,current_df])
	return TU_with_TSS_df