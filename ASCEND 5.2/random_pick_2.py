# -*- coding: utf-8 -*-
"""
### VERSION v5.1

Created on Thu Feb 21 09:13:51 2019
@author: Windows
# randomly pick n random samples from a given population, in/output in EIGENSTRAT
# 15 jan 2020: added an option to select multiple target populations
"""

import numpy as np
import argparse
from shutil import copyfile
import sys

parser = argparse.ArgumentParser(description='randomly pick n random samples from a given population, in/output in EIGENSTRAT format.')
parser.add_argument('-f', '--filePrefix', type=str, required=True, help='Prefix of the files to analyze.')
parser.add_argument('-p', '--targetPopulation', type=str, nargs='+', required=True, help='Name of the target population to analyze.')
parser.add_argument('-o', '--outfilePrefix', type=str, required=True, help='Prefix of the output file.')
parser.add_argument('-n', '--sampleSize', type=int, nargs='+', required=True, help='Sample size.')

args = parser.parse_args()
input_prefix = args.filePrefix
target_popname = args.targetPopulation
N = args.sampleSize
output_prefix = args.outfilePrefix

print(input_prefix)
print(output_prefix)
print('Target population(s): '+' '.join(target_popname))
print('Final sample size(s): '+' '.join([str(x) for x in N]))
print('Output file prefix: '+str(output_prefix))

if input_prefix==output_prefix:
    sys.exit('input prefix must differ from output prefix')

###
#####
###

IND = np.genfromtxt(input_prefix+'.ind', dtype=str)

OK = np.array([])
for i in range(len(target_popname)):
    pop = target_popname[i]
    n = N[i]
    ok = np.where(IND[:,2]==pop)[0]
    if len(ok) > n:
        ok = np.random.choice(ok, n, replace=False)
    ok = np.sort(ok)
    OK = np.concatenate((OK, ok))

OK = [int(x) for x in OK]
print(' '.join([str(x) for x in OK]))

ind = IND[OK,:]
with open(output_prefix+'.ind', 'w') as fout:
    for i in range(ind.shape[0]):
        fout.write(' '.join([str(x) for x in ind[i,:]])+'\n')
fout.close()

G = np.genfromtxt(input_prefix+'.geno', delimiter=[1]*IND.shape[0], dtype=int)
del ind, IND

G = G[:,OK]
with open(output_prefix+'.geno', 'w') as fout:
    for i in range(G.shape[0]):
        fout.write(''.join([str(x) for x in G[i,:]])+'\n')
fout.close()

copyfile(input_prefix+'.snp', output_prefix+'.snp')

