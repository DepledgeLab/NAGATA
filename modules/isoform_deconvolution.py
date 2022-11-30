import pandas as pd
import modules.help_scripts as helpers
from collections import defaultdict,Counter

def group_by_exons(unique_blocksizes,iso_group_vals):
    some_count = 0
    some_dict = defaultdict(list)
    for i in unique_blocksizes:
        fixed_i = list(map(int,i.split(',')))
        if len(some_dict) == 0:
            some_dict[some_count].append(fixed_i)
            some_count +=1
        else:
            results = []
            for k,v in some_dict.items():
                compare_with_values = any([compare_ind_blocks(fixed_i,ind_v,iso_group_vals) for ind_v in v])
                if compare_with_values == True:
                    some_dict[k].append(fixed_i)
                    results.append(True)
                    break
                else:
                    results.append(False)
            if any(results) == False:
                some_dict[some_count].append(fixed_i)
                some_count +=1
    return some_dict
    
    
def minor_formatting(blocksize_groupings):
    final_blocksize_grouping = {}
    for k,v in blocksize_groupings.items():
        for block_sz in v:
            convert_blocksz = ','.join(map(str,block_sz))
            final_blocksize_grouping[convert_blocksz] = k
    return final_blocksize_grouping

# def iso_deconv(df:'DataFrame',iso_group_vals,blocksizesum_noise_filter)-> 'DataFrame':
#     """For loop through TSS-CPAS groupings -> within groupings separate by number of exons ->
#     Use sum of blocksizes to separate isoforms -> adds grouping information to create new name
#     """
#     Final_df = pd.DataFrame()
#     for TU_TSS in sorted(set(df['new-name'])):
#         current_df = df[df['new-name'] == TU_TSS]
#         for exon_count in sorted(set(current_df['blockcount'])):
#             current_df_blocks = current_df[current_df['blockcount']==exon_count]
#             current_df_blocks['blocksize-count'] = current_df_blocks['blocksize-sum'].map(current_df_blocks['blocksize-sum'].value_counts())
#             current_df_blocks['splicing-count'] = current_df_blocks['blocksizes.new'].map(current_df_blocks['blocksizes.new'].value_counts())
#             current_df_blocks = current_df_blocks[current_df_blocks['blocksize-count'] > blocksizesum_noise_filter]
#             current_df_blocks = current_df_blocks[current_df_blocks['splicing-count'] > blocksizesum_noise_filter]
# 
#             if current_df_blocks.shape[0] > blocksizesum_noise_filter:  #
#                 blocksize_groups = helpers.identify_daisy_chain_groups(sorted(set(current_df_blocks['blocksize-sum'])),iso_group_vals)
#                 swapped_keys = helpers.swap_key_vals(dict(enumerate(blocksize_groups,1)))
#                 current_df_blocks['exon.blocksize.group'] = current_df_blocks['blocksize-sum'].map(swapped_keys)
#                 current_df_blocks['exon.blocksize.group'] = 'exon.'+current_df_blocks['blockcount'].astype(str) +'_blocksum.'+ current_df_blocks['exon.blocksize.group'].astype('str')
#                 Final_df = pd.concat([Final_df,current_df_blocks])
#     Final_df['full-id'] = Final_df['new-name'] +'-'+ Final_df['exon.blocksize.group']
#     return Final_df
def compare_ind_blocks(nagata_blockSizes,annotate_blockSizes,some_threshold):
    pass_threshold = []
    for nagata_block,anno_block in zip(nagata_blockSizes,annotate_blockSizes):
        if abs(int(nagata_block)-int(anno_block)) < some_threshold:
            pass_threshold.append(True)
        else:
            pass_threshold.append(False)
    return all(pass_threshold)
def get_groups(seedgroups,current_df,iso_group_vals):
    final_blocksize_groups = defaultdict(list)
    for i in list(set(current_df['blocksizes.new'])):
        for k,v in seedgroups.items():
            seed_value = list(map(int,v[0].split(',')))
            compare_value = list(map(int,i.split(',')))
            if compare_ind_blocks(compare_value,seed_value,iso_group_vals) == True:
                final_blocksize_groups[k].append(compare_value)
    return dict(final_blocksize_groups)
    
def iso_deconv(df:'DataFrame',iso_group_vals,blocksizesum_noise_filter)-> 'DataFrame':
    Final_df = pd.DataFrame()
    for blockcounts in set(df['new-name.ex']):
        current_df = df[df['new-name.ex'] == blockcounts]
        unique_blocksizes = [k for k,v in Counter(current_df['blocksizes.new']).items() if v>blocksizesum_noise_filter]
        grouped_blocksizes = group_by_exons(unique_blocksizes,iso_group_vals)
        seed_groupings = defaultdict(list)
        for k,v in grouped_blocksizes.items():
            edit_vv = [','.join(list(map(str,vv))) for vv in v ]
            df_identified_groups = current_df[current_df['blocksizes.new'].isin(edit_vv)]
            most_abundant_blocksize = df_identified_groups['blocksizes.new'].value_counts().head(1).index[0]
            seed_groupings[k].append(most_abundant_blocksize)
        blocksize_groups = get_groups(seed_groupings,current_df,iso_group_vals)
        swap_key_values_pairs = minor_formatting(blocksize_groups)
        current_df['blocksize_group'] = current_df['blocksizes.new'].map(swap_key_values_pairs)
        current_df = current_df[current_df['blocksize_group'].notna()]
        current_df['blocksize_group'] = current_df['blocksize_group'].astype(int).astype(str)
        Final_df = pd.concat([Final_df,current_df])
    Final_df['full-id'] = Final_df['new-name.ex'] +'-'+ Final_df['blocksize_group']
#     output_df = pd.DataFrame()
#     for unique_groups in Final_df['full-id'].unique():
#         current_df = Final_df[Final_df['full-id'] == unique_groups]
#         blockstarts_group = group_by_exons(current_df['blockstarts'].unique(),iso_group_vals)
#         blocksize_groupings = minor_formatting(blockstarts_group)
#         current_df['blockstart_group'] = current_df['blockstarts'].map(blocksize_groupings)
#         current_df['full-id'] = current_df['full-id'] +'.' +current_df['blockstart_group'].astype(str)
#         output_df = pd.concat([output_df,current_df])
    return Final_df



  
#2 def iso_deconv(df:'DataFrame',iso_group_vals,blocksizesum_noise_filter)-> 'DataFrame':
#     Final_df = pd.DataFrame()
#     for TU_TSS in sorted(set(df['new-name'])):
#         current_df = df[df['new-name'] == TU_TSS]
#         for exon_count in sorted(set(current_df['blockcount'])):
#             current_df_blocks = current_df[current_df['blockcount']==exon_count]
#             unique_blocksizes = [k for k,v in Counter(current_df_blocks['blocksizes.new']).items() if v>blocksizesum_noise_filter]
# #             print(unique_blocksizes)
# #             unique_blocksizes = current_df_blocks['blocksizes.new'].unique()
#             grouped_blocksizes = group_by_exons(unique_blocksizes,iso_group_vals)
# #             print(grouped_blocksizes)
#             blocksize_groupings = minor_formatting(grouped_blocksizes)
#             current_df_blocks['exon.blocksize.group'] = current_df_blocks['blocksizes.new'].map(blocksize_groupings)
#             current_df_blocks = current_df_blocks.dropna(subset=['exon.blocksize.group'])
#             current_df_blocks['exon.blocksize.group'] = current_df_blocks['blockcount'].astype(str) + '.blocksizes.' + current_df_blocks['exon.blocksize.group'].astype(str)
# #             intermediate_df = pd.concat([intermediate_df,current_df_blocks])
#             Final_df = pd.concat([Final_df,current_df_blocks])
#     Final_df['full-id'] = Final_df['new-name'] +'-'+ Final_df['exon.blocksize.group']
#     output_df = pd.DataFrame()
#     for unique_groups in Final_df['full-id'].unique():
#         current_df = Final_df[Final_df['full-id'] == unique_groups]
#         blockstarts_group = group_by_exons(current_df['blockstarts'].unique(),iso_group_vals)
#         blocksize_groupings = minor_formatting(blockstarts_group)
#         current_df['blockstart_group'] = current_df['blockstarts'].map(blocksize_groupings)
#         current_df['full-id'] = current_df['full-id'] +'.' +current_df['blockstart_group'].astype(str)
#         output_df = pd.concat([output_df,current_df])
#     return output_df
#     
    
    
    
    