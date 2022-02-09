import os 
import sys

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
def which(program):
#Adapted from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def check_samtools():
#Adapted from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    cmd_path=which('samtools')
    if cmd_path:
        print('SAMTOOLS in path')
        return True
    else:
        sys.stdout.write('\nERROR: NAGATA requires Samtools')
        sys.stdout.write('\n\nPlease install samtools and add it to your path following the instructions at: http://www.htslib.org/download/\n')
        return False
def check_bedtools():
#Adapted from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    cmd_path=which('bedtools')
    if cmd_path:
        print('Bedtools in path\n')
        return True
    else:
        sys.stdout.write('\nERROR: NAGATA requires Bedtools\n')
        sys.stdout.write('\n\nPlease install Bedtools and add it to your path following the instructions at: https://bedtools.readthedocs.io/en/latest/content/installation.html\n')
        return False