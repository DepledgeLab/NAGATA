import pandas as pd

def correct_last(df:'DataFrame') -> 'DataFrame':
	"""Takes in dataframe returns a corrected blocksizes column with respect to 
	corresponding CPAS-diff value . CPAS-diff is the difference between most abundant 
	CPAS in TU and this read's CPAS value
	"""
	new_blocksizes = []
	for correction,blocksizes in list(zip(df['CPAS-diff'],df['blocksizes'])):
		parse_blocks = list(map(int,blocksizes.split(',')))
		parse_blocks[-1]+= correction
		parse_blocks = ','.join(list(map(str,parse_blocks)))
		new_blocksizes.append(parse_blocks)
	df['blocksizes.new'] = new_blocksizes
	return df

def correct_first(df:'DataFrame')-> 'DataFrame':
	"""Takes in dataframe returns a corrected blocksizes column with respect to 
	corresponding TSS-diff value . CPAS-diff is the difference between most abundant 
	TSS in TU and this read's TSS value
	"""
	new_blocksizes = []
	for correction,blocksizes in list(zip(df['TSS-diff'],df['blocksizes.new'])):
		parse_blocks = list(map(int,blocksizes.split(',')))
		parse_blocks[0]-= correction
		parse_blocks = ','.join(list(map(str,parse_blocks)))
		new_blocksizes.append(parse_blocks)
	df['blocksizes.new'] = new_blocksizes
	return df

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
    
    
    