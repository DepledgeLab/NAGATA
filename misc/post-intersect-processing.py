import pandas as pd

def get_blocksize_length_(lst):
    """Uses bed file blockSize and returns overall length
    """
    return_lst = []
    for i in lst:
        if type(i) == int:
            return_lst.append(int(i))
        else:
            ints = sum([int(item) for item in i.split(',')])
            return_lst.append(ints)
    return return_lst
    
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description = 'Takes in input intersect file (between NAGATA-isoforms and "some-known"-isoform) and returns best overlapping row')

    requiredGrp = ap.add_argument_group('required arguments')
    requiredGrp.add_argument("-i",'--input_file', required=True, help="input file location")
    requiredGrp.add_argument("-o",'--output_location', required=True, help="output file location")
    args = vars(ap.parse_args())
    input_file = args['input_file']
    output_file = args['output_location']
    df = pd.read_csv(input_file,sep = '\t',header = None)
    df = df[df[9]==df[21]]
    df[25] = get_blocksize_length_(df[10])
    df[26] = get_blocksize_length_(df[22])
    df[27] = abs(df[25]-df[24])
    df[28] = abs(df[26]-df[25])
    df = df[df[27]==df[27].min()]
    df = df[df[28]==df[28].min()]
    
    df.to_csv(output_file,sep ='\t',header = None,index = None)