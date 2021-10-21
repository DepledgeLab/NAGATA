
def grouper(iterable,grouping):
    """Takes in an 'iterable' (list of TSSs or TESs) and returns the grouping of the given
    list seperated by the the value of 'grouping'
    """
    prev = None
    group = []
    for item in iterable:
        if not prev or item - prev <= grouping:
            group.append(item)
        else:
            yield group
            group = [item]
        prev = item
    if group:
        yield group
    
def swap_key_vals(initial_dict):
    """Reverses the key-value pair so actual TSS or TES sites become key and their grouping
    becomes the value
    """
    new_dict = {}
    for k,v in initial_dict.items():
        for ind_val in v:
            new_dict[ind_val] = k
    return new_dict