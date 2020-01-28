# -*- coding: utf-8 -*-
"""
@author: Remi Tournebize
# v1.0 27 jan 2020: randomly pick individuals from the dataset as outgroup
# v1.1 27 jan 20: debugged seed error if none was provided
"""

import numpy as np
import sys, argparse, random
from shutil import copyfile

###
#####
###

parser = argparse.ArgumentParser(description='Randomly pick outgroup individuals')
parser.add_argument('-p', '--parfile', type=str, required=True, help='Name of the parameter file.')
args = parser.parse_args()
parfile = args.parfile

###
#####
###

print('Reading from '+parfile)

PAR = open(parfile, 'r')
params = []
for line in PAR:
    line = line.strip()
    if line != '':
        params = params+[line]
PAR.close()

options = {
    'genotypename' : '', 
    'snpname' : '', 
    'indivname' : '', 
    'genooutfilename' : '', 
    'snpoutfilename' : '', 
    'indoutfilename' : '', 
    'outgroupsize' : '', 
    'targetpop' : '',
    'seed' : None
}
    
for param in params:
    name, value = param.split(':')
    value = value.strip()
    options[name] = value

if options['genotypename']=='':
    sys.exit('genotypename argument is missing or is empty')
else:
    genotypename = options['genotypename']

if options['snpname']=='':
    sys.exit('snpname argument is missing or is empty')
else:
    snpname = options['snpname']

if options['indivname']=='':
    sys.exit('indivname argument is missing or is empty')
else:
    indivname = options['indivname']
    
if options['genooutfilename']=='':
    sys.exit('genooutfilename argument is missing or is empty')
else:
    genooutfilename = options['genooutfilename']

if options['snpoutfilename']=='':
    sys.exit('snpoutfilename argument is missing or is empty')
else:
    snpoutfilename = options['snpoutfilename']

if options['indoutfilename']=='':
    sys.exit('indoutfilename argument is missing or is empty')
else:
    indoutfilename = options['indoutfilename']

if options['targetpop']=='':
    sys.exit('targetpop argument is missing or is empty')
else:
    targetpop = options['targetpop']

if options['outgroupsize']=='':
    sys.exit('outgroupsize argument is missing or is empty')
else:
    x = int(options['outgroupsize'])
    if x<=0:
        sys.exit('outgroupsize must not be strictly greater than 0')
    outgroupsize = x
    
if options['seed']=='':
    sys.exit('seed argument is missing or is empty')
else:
    if options['seed'] != None:
        seed = int(options['seed'])

##############################################
##############################################
##############################################

IND = np.genfromtxt(indivname+'', dtype=str)

IND = IND.astype('<U1000')

if (targetpop in IND[:,2]) == False:
    sys.exit('targetpop is not present in the ind file')

rep = np.where((IND[:,2]!=targetpop) & (IND[:,2]!='Ignore'))[0]

if outgroupsize > len(list(rep)):
    sys.exit('outgroupsize exceeds the number of available individuals')

if seed != None:
    print('Using seed for the random sampling: '+str(seed))
    random.seed(seed)
out = random.sample(list(rep), outgroupsize)

non_out = rep[ np.array([x not in out for x in rep]) ]

IND[out,2] = 'OUTGROUP'
IND[non_out,2] = 'Ignore'

keep_inds = np.where(IND[:,2] != 'Ignore')[0]

#print(keep_inds)

copyfile(snpname, snpoutfilename)

G = np.genfromtxt(genotypename, delimiter=[1]*IND.shape[0], dtype=np.int8, usecols = (keep_inds))

IND2 = IND[keep_inds,:]
del IND

with open(indoutfilename, 'w') as fout:
    for i in range(IND2.shape[0]):
        fout.write(IND2[i,0]+' '+IND2[i,1]+' '+IND2[i,2]+'\n')
fout.close()

with open(genooutfilename, 'w') as fout:
    for i in range(G.shape[0]):
        fout.write(''.join([str(x) for x in G[i,:]])+'\n')
fout.close()


