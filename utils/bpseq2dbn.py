import numpy as np
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--inputs', default='inputs', type=str, help='Path to input file in fasta format, accept multiple sequences as well in fasta format; default = ''inputs/2zzm-1-B.fasta''\n', metavar='')
parser.add_argument('--outputs',default='inputs', type=str, help='Path to output files; SPOT-RNA outputs at least three files .ct, .bpseq, and .prob files; default = ''inputs/\n', metavar='')
parser.add_argument('--rna_id', default='sample_seq', type=str, help='Name of the input sequence file\n')

args = parser.parse_args()

with open(os.path.join(args.inputs, args.rna_id + ".bpseq.unknotted")) as f:
    temp = pd.read_csv(f,comment='#', delim_whitespace=True, header=None, usecols=[0,1,2]).values
seq = temp[:,1]

pairs = [[i,j] for i,j in zip(temp[:,0], temp[:,2]) if i!=0 and j!=0 and i<j]


dbn = ['.']*temp.shape[0]
for pair in pairs:
	dbn[pair[0]-1] = '('
	dbn[pair[1]-1] = ')'

row1 = seq
row2 = np.array(dbn)
temp = np.vstack((row1, row2))

np.savetxt(os.path.join(args.outputs, args.rna_id + '.dbn'), temp, delimiter='', fmt="%s", header='>' + 'single_seq', comments='')

