import pandas as pd
import modules.help_scripts as helpers
def iso_deconv(df:'DataFrame',iso_group_vals,blocksizesum_noise_filter)-> 'DataFrame':
	"""For loop through TSS-CPAS groupings -> within groupings separate by number of exons ->
	Use sum of blocksizes to separate isoforms -> adds grouping information to create new name
	"""
	Final_df = pd.DataFrame()
	for TU_TSS in sorted(set(df['new-name'])):
		current_df = df[df['new-name'] == TU_TSS]
		for exon_count in sorted(set(current_df['blockcount'])):
			current_df_blocks = current_df[current_df['blockcount']==exon_count]
			current_df_blocks['blocksize-count'] = current_df_blocks['blocksize-sum'].map(current_df_blocks['blocksize-sum'].value_counts())
			current_df_blocks['splicing-count'] = current_df_blocks['blocksizes.new'].map(current_df_blocks['blocksizes.new'].value_counts())
			current_df_blocks = current_df_blocks[current_df_blocks['blocksize-count'] > blocksizesum_noise_filter]
			current_df_blocks = current_df_blocks[current_df_blocks['splicing-count'] > blocksizesum_noise_filter]

			if current_df_blocks.shape[0] > blocksizesum_noise_filter:  #
				blocksize_groups = helpers.identify_daisy_chain_groups(sorted(set(current_df_blocks['blocksize-sum'])),iso_group_vals)
				swapped_keys = helpers.swap_key_vals(dict(enumerate(blocksize_groups,1)))
				current_df_blocks['exon.blocksize.group'] = current_df_blocks['blocksize-sum'].map(swapped_keys)
				current_df_blocks['exon.blocksize.group'] = 'exon.'+current_df_blocks['blockcount'].astype(str) +'_blocksum.'+ current_df_blocks['exon.blocksize.group'].astype('str')
				Final_df = pd.concat([Final_df,current_df_blocks])
	Final_df['full-id'] = Final_df['new-name'] +'-'+ Final_df['exon.blocksize.group']
	return Final_df
