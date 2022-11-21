import os
import subprocess
def swap_key_vals(initial_dict):
    """Reverses the key-value pair so actual TSS or CPAS sites become key and their grouping
    becomes the value
    """
    new_dict = {}
    for k,v in initial_dict.items():
        for ind_val in v:
            new_dict[int(ind_val)] = k
    return new_dict

def identify_daisy_chain_groups(cpas_list,daisy_chain_val)->list:
    """Takes in list of CPAS peaks returns grouped CPAS sites (with considerations to daisy chaining)
    """
    initial_list = [[cpas_list.pop(0)]]
    for cpas_val in cpas_list:
        if cpas_val - initial_list[-1][-1] < daisy_chain_val and abs(cpas_val - initial_list[-1][-1]) < daisy_chain_val:
            initial_list[-1].append(cpas_val)
        else:
            initial_list.append([cpas_val])
    return initial_list

def create_tmp_files(input_bam:'string',output_location:'string'):
    """Takes in an imput file and creates a temp file used for actual processing
    """
#     try:
#         print(override, '!!!!')
#         os.makedirs(output_location,exist_ok=override)
#     except:
#         print("Variable x is not defined")
    os.makedirs(output_location,exist_ok=True)
    ## Samtools to get read-ids, cigar strings, orientation into a temporary file
    samtools_p1 = subprocess.Popen(f"samtools view {input_bam}".split(' '),stdout=subprocess.PIPE)
    fout_samtools = open(output_location +'/seq-cigar-orient.tmp', 'wb')
    samtools_p2 = subprocess.run(['cut','-f1,2,6'], stdin=samtools_p1.stdout, stdout=fout_samtools)
    ## Bedtools to convert BAM file to BED
    fout_bedtools = open(output_location +'/1.raw-alignment.bed','wb')
    bedtools_p1 = f"bedtools bamtobed -bed12 -i {input_bam}".split(" ")
    bedtools_p2 = subprocess.run(bedtools_p1,stdout=fout_bedtools)
    print('temp files created')
