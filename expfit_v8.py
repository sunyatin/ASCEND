# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 15:12:07 2019
@author: Windows
#220319 15:16 added model comparison
#250319 15:01 maj jackknife procedure (now weighted) and changed output format, verified by comparing results with those obtained using dowtjack given by priya
#2803201 11:15 added rearrangement for jacquelin's method in case of negative values
#29032019 11:10 constrained to handle exp1 fit for bg LD substraction
#15052019 added try in case of errors, now always output the pop
#10062019: replaced 1/A and 1/(A+c) computation by -Alog(abs(c)) and 1/(-Alog(abs(c))) respectively
#02082019: simplified code, checked with R its working OK
#11082019: added a constraint on the chromosomes present in the chrom (-n) file that should be analyzed
"""

version = 'v8'

import numpy as np
import io, os, sys, argparse, warnings
from scipy.optimize import curve_fit, OptimizeWarning
np.seterr(all='ignore')
warnings.simplefilter("error", OptimizeWarning)


parser = argparse.ArgumentParser(description='Calcule allele sharing statistics.')
parser.add_argument('-f', '--inputFile', type=str, required=True, help='Input file (with the .out extension).')
parser.add_argument('-p', '--targetPopulation', type=str, required=True, help='Name of the target population to analyze.')
parser.add_argument('-o', '--outfilePrefix', type=str, required=True, help='Prefix of the output file.')
parser.add_argument('-n', '--blockSizeFile', type=str, required=True, help='File containing two tab-separated columns: chromosome number, number of SNPs. Used for jackknife.')
parser.add_argument('-minD', '--minD', type=float, required=True, help='Minimum genetic distance.')
parser.add_argument('-maxD', '--maxD', type=float, default = 20.0, help='Maximum genetic distance.')
parser.add_argument('--noBgLDSubstraction', action='store_true', default = False, help='Add this option if you do not want to substract covariance by background LD.')

args = parser.parse_args()

inputFile = args.inputFile
target_popname = args.targetPopulation
minD = args.minD
maxD = args.maxD
use_substracted = not args.noBgLDSubstraction
out = args.outfilePrefix
blockSizeFile = args.blockSizeFile

print('================================================================================================')
print('expfit '+version+'   |   '+'Input file: '+inputFile+'   |   Target population:  '+target_popname+'   |   Blocksize:  '+blockSizeFile+'   |   Output:  '+out)
print('Distance: '+str(minD)+' to '+str(maxD)+' cM   |   Do background LD substraction: '+str(use_substracted))

##
#####
##


def jacquelin_exp1d(XY, min_D_cM, max_D_cM):

    XY = XY[XY[:,0]>=min_D_cM,:]
    XY = XY[XY[:,0]<=max_D_cM,:]
    
    x = XY[:,0]
    y = XY[:,1]
    n = XY.shape[0]

    S = [0]
    for k in np.arange(1, n):
        S += [ S[-1] + 1/2 * (y[k]+y[k-1])*(x[k]-x[k-1]) ]
              
    S = np.asarray(S)
    x = np.asarray(x)
    y = np.asarray(y)

    M1 = [
        [
            sum( (x-x[0])**2 ),
            sum( (x-x[0])*S )
            
        ],
        [
            sum( (x-x[0])*S ),
            sum( S**2 )
        ]]
        
    M2 = [
          sum( (y-y[0])*(x-x[0]) ),
          sum( (y-y[0])*S )
          ]
        
    M1 = np.array(M1)
    M2 = np.array(M2)
    
    K = np.dot(np.linalg.inv(M1), M2)
    c = K[1]
    
    N1 = [
        [
            n,
            sum( np.exp(c*x) )
            
        ],
        [
            sum( np.exp(c*x) ),
            sum( np.exp(2*c*x) )
        ]]
        
    N2 = [
          sum( y ),
          sum( y*np.exp(c*x) )
          ]
  
    AB = np.dot(np.linalg.inv(N1), N2)
          
    a = AB[0]
    b = AB[1]
  
    def exp(x, a, b, c):
        return a + b*np.exp(c*x)
        
    if np.isnan(a):
        return None

    try:
        popt, pcov = curve_fit(exp, x, y, maxfev = 5000, p0 = [a, b, c])
        ret = np.array([popt[1], popt[2]/-2, popt[0], 0, 0])
        ret = np.append(ret, [1/ret[0], -1.0*ret[0]*np.log(abs(ret[2])), 1/(-1.0*ret[0]*np.log(abs(ret[2])))])
    except:
        ret = None

    return ret


def expfit_1D_jackknife(input, blocksizes, min_D_cM = 0, max_D_cM = 30, use_substracted = True):
    
    issue = False
    
    F = np.genfromtxt(input, delimiter='\t', dtype='float', skip_header=1)
    F_q = np.asarray([int(x) for x in F[:,0]])
    chrom = [int(x) for x in np.unique(F[:,0]).tolist()]

    intersection = np.intersect1d(np.array(chrom), blocksizes[:,0]).tolist()
    
    print('\tChromosomes in input file:\n'+' '.join([str(int(x)) for x in np.unique(F[:,0]).tolist()]))
    print('\tChromosomes at the intersection of the two files:\n'+' '.join([str(x) for x in intersection]))
    
    if len(intersection) == 0:
        sys.exit('\n*** Error. The intersection of chromosomes between the two files is empty.')
    
    chrom = intersection
    ok = [blocksizes[i,0] in intersection for i in range(blocksizes.shape[0])]
    blocksizes = blocksizes[ok,:]
    
    # general mean
    XY, first = [], True
    for q in chrom:
        xy = F[F_q==q,:]
        if first:
            XY = xy
            first = False
        else:
            XY[:,3:7] += xy[:,3:7]
    
    if use_substracted:
        XY = np.column_stack((XY[:,1], np.divide(XY[:,5], XY[:,6])))
    else:
        XY = np.column_stack((XY[:,1], np.divide(XY[:,3], XY[:,6])))
    
    params = jacquelin_exp1d(XY, min_D_cM, max_D_cM)

    if params is None:
        general_mean = np.array([None]*8)
    else:
        general_mean = params

    # per chromosome estimate
    PARAMS = []
    for i in range(len(chrom)):
        use = chrom.copy()
        if len(chrom) > 1:
            use.pop(i)
        
        XY, first = [], True
        for q in use:
            xy = F[F_q==q,:]
            if first:
                XY = xy
                first = False
            else:
                XY[:,3:7] += xy[:,3:7]

        if use_substracted:
            XY = np.column_stack((XY[:,1], np.divide(XY[:,5], XY[:,6])))
        else:
            XY = np.column_stack((XY[:,1], np.divide(XY[:,3], XY[:,6])))
        
        params = jacquelin_exp1d(XY, min_D_cM, max_D_cM)
        
        if params is None:
            PARAMS.append(np.array([None]*8))
            issue = True
        else:
            PARAMS.append(params.tolist())
    
    PARAMS = np.asmatrix(PARAMS)
    
    N = []
    for q in chrom:
        N += [int(blocksizes[blocksizes[:,0]==int(q),1])]
    N = np.asmatrix(N).T
        
    PARAMS = np.hstack((PARAMS, N))
    N = PARAMS[:,-1].A1
    PARAMS = PARAMS[:,0:-1]
    
    if issue is False:
        # delete-m jackknife
        n_tot = np.sum(N)
        g = PARAMS.shape[0] 
        MEAN = []
        SD = []
        for j in range(PARAMS.shape[1]):
            x = PARAMS[:,j].A1
            mean = g*general_mean[j] - np.sum( (1-N/n_tot)*x )
            H = n_tot/N
            term = np.sum( (1-N/n_tot)*x )
            term2 = H*general_mean[j] - (H-1)*x - g*general_mean[j] + term
            var = 1/g * np.sum( 1/(H-1) * term2**2   )
            sd = np.sqrt(var)     
            MEAN += [mean]
            SD += [sd]
        
        PARAMS = np.vstack((PARAMS, general_mean, MEAN, SD))
        PARAMS = np.hstack((PARAMS, np.asmatrix(np.append(N, (None, None, None))).T ))
    else:
        PARAMS = np.vstack((PARAMS, general_mean, [None]*8, [None]*8))
        PARAMS = np.hstack((PARAMS, np.asmatrix(np.append(N, (None, None, None))).T ))
    
    return PARAMS


###
######
###
    
with open(out+'.fit', 'w') as fout:
    fout.write('method\tpop\tchromosome\tA\tt\tc\tk2\tt2\t1/A\t-Alog(abs(c))\t-1/Alog(abs(c))\tsize\n')

    if not os.path.exists(blockSizeFile):
        fout.close()
        sys.exit('\n*** Error. File does not exist: '+blockSizeFile+'\n================================================================================================')
        
    if not os.path.exists(inputFile):
        fout.close()
        sys.exit('\n*** Error. File does not exist: '+inputFile+'\n================================================================================================')
    
    ###
    ######
    ###
    
    s = io.StringIO(open(blockSizeFile).read().replace('\t', ' '))
    bs = np.genfromtxt(s, dtype=int, delimiter = ' ')
    if bs.ndim == 0:
        fout.close()
        sys.exit('\n*** Error. The block size file is empty.')
    elif bs.ndim == 1:
        bs = np.reshape(bs, (-1, 2))
        
    R = expfit_1D_jackknife(input = inputFile,
            blocksizes = bs,
            min_D_cM = minD,
            max_D_cM = maxD,
            use_substracted = use_substracted)
    
    for i in range(R.shape[0]-3):
        fout.write('Jacquelin\t'+target_popname+'\t'+str(i+1)+'\t'+'\t'.join([str(x) for x in R[i,:].A1])+'\n')
    fout.write('Jacquelin\t'+target_popname+'\t'+'MEAN'+'\t'+'\t'.join([str(x) for x in R[R.shape[0]-3,:].A1])+'\n')
    fout.write('Jacquelin\t'+target_popname+'\t'+'JK.MEAN'+'\t'+'\t'.join([str(x) for x in R[R.shape[0]-2,:].A1])+'\n')
    fout.write('Jacquelin\t'+target_popname+'\t'+'JK.SD'+'\t'+'\t'.join([str(x) for x in R[R.shape[0]-1,:].A1])+'\n')

fout.close()

print('================================================================================================')

