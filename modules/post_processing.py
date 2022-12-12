import pandas as pd

def dataframe_editing(df) -> 'BED':
	""" Takes in processed dataframe and returns an edited BED file
	"""
# 	df = df.drop(['new-name','exon.blocksize.group','TSS-diff','CPAS-diff','most_abund_CPAS','TSS-count','TSS.unique','TSS-group','most_abund_TSS_in_TU'],axis =1)
	
	df['TSS-unique-isoform-count-per-TU'] = df['full-id'].map(df['full-id'].value_counts())
	df['most_abund_CPAS'] = df['end'].map(df['end'].value_counts())
	df['TSS-abundance-per-TU']= round((df['TSS-unique-isoform-count-per-TU']/ df['most_abund_CPAS'])*100,2)
	df['blocksizes'] = df['blocksizes.new']
	
	df['score'] = df['TSS-unique-isoform-count-per-TU']
	df['thickstart'] = df['start']
	df['thickend'] = df['end']
	df['name'] = df['full-id']
	return df