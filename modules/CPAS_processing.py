from collections import Counter

def filter_CPAS_noise(cpas_list,min_cpas_count)->dict:
    """Takes in CPAS list returns the CPAS sites that have a count with a minimum threshold
    """
    cpas_count_dict = Counter(cpas_list)
    cpas_count = sum(cpas_count_dict.values())
    cpas_unique_before = len(cpas_count_dict.keys())
    pass_cpas_filter = { k:v for k,v in cpas_count_dict.items() if v>= min_cpas_count}
    pass_count = sum(pass_cpas_filter.values())
    cpas_unique_after = pass_cpas_filter.keys()
    #print(f'Total reads before: {cpas_count}')
    #print(f'Total reads after CPAS noise filtering: {pass_count}')
    #print(f'Total unique CPAS before: {cpas_unique_before}')
    #print(f'Total unique CPAS after CPAS noise filtering: {len(cpas_unique_after)}')
    return pass_cpas_filter
    
