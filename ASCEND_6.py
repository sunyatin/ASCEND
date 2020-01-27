# -*- coding: utf-8 -*-
"""
### VERSION v6
Created on Thu Feb 21 09:13:51 2019
@author: Remi Tournebize
# 24 jun 2019: created a mean to compute the correlation instead of the covariance in allele sharing
# 25 jul 2019: debugged (before we added the SNP in the count if its cor value was Nan or Inf)
# 30 jul 2019: compute MAF across the populations
# 09 aug 2019: handle missing data
# 10 aug 2019: added the quick mode in case there is no missing data
# 10 aug 2019: the script has been carefully checked with test_profile_script.R with various maxPropNA and pNA: it works fine
# 11 aug 2019: significantly improved the speed of execution up to 3-fold and also better memory management than in v5, it has also been checked with test_profile_script.R and validated the test
# 11 sep 2019: added --haploid and --pseudodiploid modes
# 14 oct 2019: added a --chrom option to run the analysis on a particular chromosome
# 26 jan 2020: use an external file to get the parameters
# 27 jan 2020: added the expfit directly to the script and added the random outgroup selection too
"""

import numpy as np
import time, sys, warnings, io, os, argparse
from scipy.optimize import curve_fit, OptimizeWarning
np.seterr(all='ignore')
warnings.simplefilter('error', OptimizeWarning)
np.seterr(divide='ignore', invalid='ignore')

###
#####
###

parser = argparse.ArgumentParser(description='ASCEND v6')
parser.add_argument('-p', '--parfile', type=str, required=True, help='Name of the parameter file.')
args = parser.parse_args()
parfile = args.parfile

###
#####
###

def number_shared_alleles(pop1, pop2):
    dif = 2 - abs(pop1 - pop2)
    prod = pop1 * pop2
    dif[prod==1] = 1
    dif[prod==81] = -9
    # the missing data will be encoded as negative values (-5, -6, -7 or -9)
    return dif.astype(np.int8)

def compute_correlation(ASD_MATRIX, i, ASD_MATRIX_noNA, snpIndices, max_proportion_NA):
    SNPs = ASD_MATRIX[(i+1):ASD_MATRIX.shape[0],:]
    SNPs = SNPs[snpIndices,:]
    SNP1 =  np.tile(ASD_MATRIX[i,:], (SNPs.shape[0], 1))
    del ASD_MATRIX
    # use broadcasting here instead of tile
    
    SNPs_noNA = ASD_MATRIX_noNA[(i+1):ASD_MATRIX_noNA.shape[0],:]
    SNPs_noNA = SNPs_noNA[snpIndices,:]
    SNP1_noNA = np.tile(ASD_MATRIX_noNA[i,:], (SNPs_noNA.shape[0], 1))
    del ASD_MATRIX_noNA

    PROD = SNPs * SNP1
    nMaxPairs = PROD.shape[1]
    PROD_noNA = SNPs_noNA * SNP1_noNA
    del SNP1_noNA, SNPs_noNA

    nPairs = PROD_noNA.sum(1)
    
    SNPs *= PROD_noNA
    SNP1 *= PROD_noNA
    PROD *= PROD_noNA

    mean_SNP1 = SNP1.sum(1) / nPairs
    mean_SNPs = SNPs.sum(1) / nPairs
    
    r_ = PROD.sum(1)/nPairs - (mean_SNP1*mean_SNPs)
    
    PROD_noNA = PROD_noNA.transpose()
    DIF = (SNP1.transpose() - mean_SNP1) * PROD_noNA
    std_SNP1 = np.sqrt( (DIF**2).sum(0) / (nPairs - 1*0) )
    DIF = (SNPs.transpose() - mean_SNPs) * PROD_noNA
    std_SNPs = np.sqrt( (DIF**2).sum(0) / (nPairs - 1*0) )
    del PROD_noNA, mean_SNP1, SNP1, mean_SNPs, SNPs
    
    r_ = r_ / (std_SNP1 * std_SNPs)
    del std_SNP1, std_SNPs
    
    bads1 = np.where(nPairs < nMaxPairs*(1-max_proportion_NA))
    bads2 = np.where(np.isnan(r_))
    bads = np.unique(np.concatenate((bads1, bads2), axis=None))
    del bads1, bads2
    
    return((r_, bads))


def calculate_allele_sharing(input_prefix, 
                             output_prefix,
                             target_popname, 
                             input_geno_is_diploid, 
                             out_popname = None, 
                             minMAF = 0.05,
                             stepD_cM = 0.5, 
                             minD_cM = 0, 
                             maxD_cM = None,
                             input_distance_in_cM = False,
                             chrom_to_analyze = None):
    start = time.time()
    if out_popname is not None:
        print('Analyzing '+target_popname+' and using an outgroup population: '+out_popname)
    else:
        print('Analyzing only '+target_popname+', without outgroup.')
    corr_decimals = 1e-10
    
    D_FULL = np.genfromtxt(input_prefix[1]+'', dtype=float, usecols = (1,2))
    CHR = D_FULL[:,0].astype(int)
    D_FULL = D_FULL[:,1]
    if input_distance_in_cM is False:
        print('Converting Morgans input into centiMorgans.')
        D_FULL = 100.0 * D_FULL
    POP = np.genfromtxt(input_prefix[2]+'', dtype=str, usecols = 2)
    
    if out_popname == None:
        print('There are '+str(np.sum(POP==target_popname))+' target samples\n')
    else:
        print('There are '+str(np.sum(POP==target_popname))+' target and '+str(np.sum(POP==out_popname))+' outgroup samples\n')
    
    if chrom_to_analyze is None:
        chr_values = np.unique(CHR).tolist()
    else:
        chr_values = np.array([chrom_to_analyze])
    print('Chromosomes: '+' '.join([str(x) for x in chr_values])+'\n')
    
    RR, first = [], True
    for chrom in chr_values:
        foc = np.where(CHR == chrom)[0]

        r0, r1 = foc[0], foc[-1]  
        print('\n\n>> Chrom:    '+str(chrom)+'   -~-   Range: '+str(r0)+' => '+str(r1))
        D = D_FULL[foc]
        del foc
        G = np.genfromtxt(input_prefix[0]+'', delimiter=[1]*len(POP), dtype=int,
                               skip_header = r0, max_rows = r1-r0+1)
        
        if input_geno_is_diploid is False:
            G[G==1] = 2
            
        if pseudodiploidize is True:
            G[G==1] = np.random.choice([0,2], size = len(G[G==1]), replace = True)
    
        if target_popname not in POP:
            sys.exit("ERROR! The target population you specified is not present in the ind file")

        there_is_NA_in_dataset = False

        if out_popname is not None:
            G_out = G[:,POP==out_popname]
            if np.any(G_out==9):
                there_is_NA_in_dataset = True
            nh_2 = G_out.shape[1]
            if out_popname not in POP:
                sys.exit("ERROR! The outgroup population you specified is not present in the ind file")
    
        G = G[:,POP==target_popname]
        if np.any(G==9):
                there_is_NA_in_dataset = True
                
        if there_is_NA_in_dataset:
            print('Slow Mode')
        
            # filter on the MAF
            nh_1 = G.shape[1]*2
            G_no9 = G.copy()
            G_no9[G_no9<=2] = 10
            G_no9 -= 9
            NH1 = G_no9.sum(1) * 2
            G0 = G.copy()
            G0 = G0 * G_no9
            del G_no9
            
            if out_popname is not None:
                nh_2 = G_out.shape[1]*2
                G_out_no9 = G_out.copy()
                G_out_no9[G_out_no9<=2] = 10
                G_out_no9 -= 9
                NH2 = G_out_no9.sum(1) * 2
                G0_out = G_out.copy()
                G0_out = G0_out * G_out_no9
                del G_out_no9
                # maf
                maf = 1.0 * (G0.sum(1)+G0_out.sum(1)) / (NH1+NH2)
                del G0_out, NH2
            else:
                maf = 1.0 * G0.sum(1) / NH1
            del G0, NH1
            maf[maf>0.5] = 1 - maf[maf>0.5]
            goods = np.where(maf>minMAF)
            del maf
            nSNP_raw = G.shape[0]
            D = D[goods]
            G = G[goods,:][0]

            if out_popname is not None:
                G_out = G_out[goods,:][0]     
        
            if maxD_cM is None:
                maxD_cM = D[len(D)]
            
            nSNP = G.shape[0]
            print('Raw nSNP: '+str(nSNP_raw)+'   MAF-filtered nSNP: '+str(nSNP))
                    
            # allele sharing calculation   
            ASD1 = np.zeros((nSNP, int( (nh_1/2)*(nh_1/2-1)/2) ))
            ASD1 = ASD1.astype(np.int8)
            i = 0
            for j in range(0, G.shape[1]):
                for k in range(j+1, G.shape[1]):
                    ASD1[:,i] = number_shared_alleles(G[:,j], G[:,k])
                    i += 1
            ASD1_noNA = np.clip(ASD1, -1, 0) + 1

            if out_popname is not None:
                ASD2 = np.zeros((nSNP, int((nh_2/2)*(nh_1/2)) ))
                ASD2 = ASD2.astype(np.int8)
                i = 0
                for j in range(0, G.shape[1]):
                    for k in range(0, G_out.shape[1]):
                        ASD2[:,i] = number_shared_alleles(G[:,j], G_out[:,k])
                        i += 1
                ASD2_noNA = np.clip(ASD2, -1, 0) + 1
    
            del G
            if out_popname is not None:
                del G_out
                   
            # populate the bins
            bins_left_bound = np.arange(minD_cM, maxD_cM, step=stepD_cM)
            n_bins = len(bins_left_bound)
            print('There will be '+str(n_bins)+' bins.')
   
            BINS_R = np.asarray([0]*n_bins).astype(np.float)
            BINS_Rbg = np.asarray([0]*n_bins).astype(np.float)
            BINS_Rcor = np.asarray([0]*n_bins).astype(np.float)
            BINS_N = np.asarray([0]*n_bins).astype(np.int64)
            v_last = 0
            for i in range(nSNP):
            
                dd = D[(i+1):nSNP]
                d = D[i]
                d_ = abs(dd - d)
                indices = ( (d_ - minD_cM + corr_decimals) // stepD_cM).astype(np.int32)
                goods_Dref = np.where((0<=indices) & (indices<=(n_bins-1)))[0]
                indices = indices[goods_Dref]
                
                (r_, inf1) = compute_correlation(ASD1, i, ASD1_noNA, goods_Dref, max_proportion_NA)

                if out_popname is not None:
                    (r_bg_, inf2) = compute_correlation(ASD2, i, ASD2_noNA, goods_Dref, max_proportion_NA)
                    inf = np.unique(np.concatenate((inf1, inf2), axis=None))
                    r_[inf] = np.array([0]*len(inf))
                    r_bg_[inf] = np.array([0]*len(inf))
                    np.add.at(BINS_R, indices, r_)
                    np.add.at(BINS_Rbg, indices, r_bg_)
                    r_sub_ = (r_ - r_bg_)
                    np.add.at(BINS_Rcor, indices, r_sub_)
                    indices = np.delete(indices, inf)
                    del inf1, inf2, inf, r_, r_bg_
                    np.add.at(BINS_N, indices, 1)
                else:
                    r_[inf1] = np.array([0]*len(inf1))
                    np.add.at(BINS_R, indices, r_)
                    indices = np.delete(indices, inf1)
                    del inf1, r_
                    np.add.at(BINS_N, indices, 1)
            
                v_i = int(i*100/nSNP)
                if v_i > v_last:
                    print(str(v_i)+'%', end=' ', flush=True)
                    v_last = v_i
                    if out_popname is None:
                        BINS_Rbg = ['NA']*len(BINS_R)
                        BINS_Rcor = ['NA']*len(BINS_R)
                    R = np.column_stack(([chrom]*len(BINS_R),
                                   bins_left_bound,
                                   bins_left_bound + stepD_cM/float(2),
                                   BINS_R,
                                   BINS_Rbg,
                                   BINS_Rcor,
                                   BINS_N))
            if first:
                RR = R
                first = False
            else:
                RR = np.concatenate((RR, R), axis=0)
            
            with open(output_prefix+'', 'w') as fout:
                fout.write('chrom\tbin.left.bound\tbin.center\tcor.pop\tcor.bg\tcor.substracted\tn.pairs\n')
                for i in range(RR.shape[0]):
                    fout.write(str(int(RR[i,0]))+'\t')
                    fout.write('{:.5f}'.format(RR[i,1])+'\t')
                    fout.write('{:.5f}'.format(RR[i,2])+'\t')
                    fout.write('\t'.join([str('{:.5f}'.format(x)) for x in RR[i,3:]])+'\n')
            fout.close()
        
        else:
            print('Quick Mode')
        
            nh_1 = G.shape[1]
            nh_1 *= 2
            if out_popname is not None:
                nh_2 *= 2
            if out_popname is not None:
                maf = 1.0 * (G.sum(1)+G_out.sum(1)) / (nh_1+nh_2)
            else:
                maf = 1.0 * G.sum(1) / nh_1
            maf[maf>0.5] = 1-maf[maf>0.5]
            goods = np.where(maf>minMAF)
            del maf
            nSNP_raw = G.shape[0]
            D = D[goods]
            G = G[goods,:][0]
                
            if out_popname is not None:
                G_out = G_out[goods,:][0]     
        
            if maxD_cM is None:
                maxD_cM = D[len(D)]
            
            nSNP = G.shape[0]
            print('Raw nSNP: '+str(nSNP_raw)+'   MAF-filtered nSNP: '+str(nSNP))
                    
            # allele sharing calculation
            ASD1 = np.zeros((nSNP, int( (nh_1/2)*(nh_1/2-1)/2) ))
            ASD1 = ASD1.astype(np.int8)
            i = 0
            for j in range(0, G.shape[1]):
                for k in range(j+1, G.shape[1]):
                    ASD1[:,i] = number_shared_alleles(G[:,j], G[:,k])
                    i += 1
        
            if out_popname is not None:
                ASD2 = np.zeros((nSNP, int((nh_2/2)*(nh_1/2)) ))
                ASD2 = ASD2.astype(np.int8)
                i = 0
                for j in range(0, G.shape[1]):
                    for k in range(0, G_out.shape[1]):
                        ASD2[:,i] = number_shared_alleles(G[:,j], G_out[:,k])
                        i += 1
                    
            del G
            if out_popname is not None:
                del G_out
                    
            # populate the bins
            bins_left_bound = np.arange(minD_cM, maxD_cM, step=stepD_cM)
            n_bins = len(bins_left_bound)
            print('There will be '+str(n_bins)+' bins.')
            
            BINS_R = np.asarray([0]*n_bins).astype(np.float)
            BINS_Rbg = np.asarray([0]*n_bins).astype(np.float)
            BINS_Rcor = np.asarray([0]*n_bins).astype(np.float)
            BINS_N = np.asarray([0]*n_bins).astype(np.int64)
            v_last = 0
            for i in range(nSNP):
                dd = D[(i+1):nSNP]
                d = D[i]
                d_ = abs(dd - d)
                indices = ( (d_ - minD_cM + corr_decimals) // stepD_cM).astype(np.int32)
                goods_Dref = np.where((0<=indices) & (indices<=(n_bins-1)))[0]
                indices = indices[goods_Dref]

                aadd = ASD1[(i+1):nSNP,:]
                aadd = aadd[goods_Dref,:]
                ad = ASD1[i,:]
                r_ = (aadd*ad).mean(1) - aadd.mean(1)*ad.mean()
                r_ = r_ / (aadd.std(1) * ad.std())
                inf1 = np.where(np.isnan(r_))
                
                if out_popname is not None:
                    aadd_bg = ASD2[(i+1):nSNP,:]
                    aadd_bg = aadd_bg[goods_Dref,:]
                    ad_bg = ASD2[i,:]
                    r_bg_ = (aadd_bg*ad_bg).mean(1) - aadd_bg.mean(1)*ad_bg.mean()
                    r_bg_ = r_bg_ / (aadd_bg.std(1) * ad_bg.std())
                    inf2 = np.where(np.isnan(r_bg_))

                if out_popname is not None:
                    inf = np.unique(np.append(inf1, inf2))
                    r_[inf] = np.array([0]*len(inf))
                    r_bg_[inf] = np.array([0]*len(inf))
                    np.add.at(BINS_R, indices, r_)
                    np.add.at(BINS_Rbg, indices, r_bg_)
                    r_sub_ = (r_ - r_bg_)
                    np.add.at(BINS_Rcor, indices, r_sub_)
                    del aadd, ad, aadd_bg, ad_bg, r_, r_bg_, r_sub_
                    indices = np.delete(indices, inf)
                    np.add.at(BINS_N, indices, 1)
                
                else:
                    r_[inf1] = np.array([0]*len(inf1))
                    np.add.at(BINS_R, indices, r_)
                    del aadd, ad, r_
                    indices = np.delete(indices, inf1)
                    np.add.at(BINS_N, indices, 1)
                
                v_i = int(i*100/nSNP)
                if v_i > v_last:
                    print(str(v_i)+'%', end=' ', flush=True)
                    v_last = v_i
                    
                    if out_popname is None:
                        BINS_Rbg = ['NA']*len(BINS_R)
                        BINS_Rcor = ['NA']*len(BINS_R)
                        
                    R = np.column_stack(([chrom]*len(BINS_R),
                                   bins_left_bound,
                                   bins_left_bound + stepD_cM/float(2),
                                   BINS_R,
                                   BINS_Rbg,
                                   BINS_Rcor,
                                   BINS_N))

            if first:
                RR = R
                first = False
            else:
                RR = np.concatenate((RR, R), axis=0)
            
            with open(output_prefix+'', 'w') as fout:
                fout.write('chrom\tbin.left.bound\tbin.center\tcor.pop\tcor.bg\tcor.substracted\tn.pairs\n')
                for i in range(RR.shape[0]):
                    fout.write(str(int(RR[i,0]))+'\t')
                    fout.write('{:.5f}'.format(RR[i,1])+'\t')
                    fout.write('{:.5f}'.format(RR[i,2])+'\t')
                    fout.write('\t'.join([str('{:.5f}'.format(x)) for x in RR[i,3:]])+'\n')
            fout.close()
        
    print('\nIt took:   '+str(round((time.time()-start)/60,2))+' min\n')
    
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
        yfit = exp(x, popt[0], popt[1], popt[2])
        NRMSD = 1./(np.max(yfit)-np.min(yfit)) * np.sqrt(np.mean( (yfit-y)**2 ))
        ret = np.array([popt[1], popt[2]/-2, popt[0], NRMSD])
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
        general_mean = np.array([None]*4)
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
            PARAMS.append(np.array([None]*4))
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
        PARAMS = np.vstack((PARAMS, general_mean, [None]*4, [None]*4))
        PARAMS = np.hstack((PARAMS, np.asmatrix(np.append(N, (None, None, None))).T ))
    
    return PARAMS

##
### script
##

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
    'targetpop' : '', 
    'outpop' : None, 
    'outputprefix' : '', 
    'minmaf' : 0, 
    'haploid' : 'NO', 
    'dopseudodiploid' : 'NO', 
    'binsize' : 0, 
    'mindis' : 0.1, 
    'maxdis' : 30, 
    'morgans' : 'NO', 
    'maxpropmissing' : 1,
    'chrom' : None,
    'blocksizename': '',
    'onlyfit': 'NO'
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
    
if options['targetpop']=='':
    sys.exit('targetpop argument is missing or is empty')
else:
    target_popname = options['targetpop']
    
if options['outpop']=='':
    sys.exit('outpop argument is missing or is empty')
else:
    out_popname = options['outpop']
    
if options['outputprefix']=='':
    sys.exit('outputprefix argument is missing or is empty')
else:
    output_prefix = options['outputprefix']
    
if options['minmaf']=='':
    sys.exit('minmaf argument is missing or is empty')
else:
    x = float(options['minmaf'])
    if x>1:
        sys.exit('minmaf must not be strictly greater than 1')
    minMAF = x
    
if options['haploid']=='':
    sys.exit('haploid argument is missing or is empty')
else:
    if options['haploid'] == 'YES':
        input_geno_is_diploid = False
    elif options['haploid'] == 'NO':
        input_geno_is_diploid = True
    else:
        sys.exit('haploid argument you provided is unrecognized, must be either YES or NO')
        
if options['onlyfit']=='':
    sys.exit('onlyfit argument is missing or is empty')
else:
    if options['onlyfit'] == 'YES':
        onlyfit = True
    elif options['onlyfit'] == 'NO':
        onlyfit = False
    else:
        sys.exit('onlyfit argument you provided is unrecognized, must be either YES or NO')

if options['dopseudodiploid']=='':
    sys.exit('dopseudodiploid argument is missing or is empty')
else:
    if options['dopseudodiploid'] == 'YES':
        pseudodiploidize = True
    elif options['dopseudodiploid'] == 'NO':
        pseudodiploidize = False
    else:
        sys.exit('dopseudodiploid argument you provided is unrecognized, must be either YES or NO')

if options['binsize']=='':
    sys.exit('binsize argument is missing or is empty')
else:
    stepD_cM = float(options['binsize'])
    
if options['mindis']=='':
    sys.exit('mindis argument is missing or is empty')
else:
    minD_cM = float(options['mindis'])
    
if options['maxdis']=='':
    sys.exit('maxdis argument is missing or is empty')
else:
    maxD_cM = float(options['maxdis'])
    
if options['morgans']=='':
    sys.exit('morgans argument is missing or is empty')
else:
    if options['morgans'] == 'YES':
        input_distance_in_cM = False
    elif options['morgans'] == 'NO':
        input_distance_in_cM = True
    else:
        sys.exit('morgans argument you provided is unrecognized, must be either YES or NO')
 
if options['maxpropmissing']=='':
    sys.exit('maxpropmissing argument is missing or is empty')
else:
    x = float(options['maxpropmissing'])
    if x>1 or x<0:
        sys.exit('maxpropmissing must be contained in the range 0-1')
    max_proportion_NA = x
     
if options['chrom']=='':
    sys.exit('chrom argument is missing or is empty')
else:
    if options['chrom'] != None:
        cta = int(options['chrom'])
    else:
        cta = None
        
if options['blocksizename']=='':
    sys.exit('blocksizename argument is missing or is empty')
else:
    blockSizeFile = options['blocksizename']
        
if not os.path.exists(blockSizeFile):
    sys.exit('\n*** Error. File does not exist: '+blockSizeFile+'\n')
    
####
    
print('ASCEND version 6\n')
print(genotypename+' '+snpname+' '+indivname)
print('Output file name: '+output_prefix)
print('Target population: '+target_popname)
print('Out population: '+str(out_popname))
print('minMAF: '+str(minMAF))
if input_geno_is_diploid:
    print('Assumes the genotypes are diploids')
else:
    print('Assumes the genotypes are haploids')
if pseudodiploidize:
    print('Genotypes are pseudodiploidized')
print('D:'+str(minD_cM)+' '+str(maxD_cM)+' by '+str(stepD_cM)+' cM')
print('Input distance in cM: '+str(input_distance_in_cM))
print('Maximum proportion of allele sharing values that can be missing: '+str(max_proportion_NA))
print('Block size file: '+blockSizeFile)
print('================================================================================')

if input_geno_is_diploid is False and pseudodiploidize is True:
    sys.exit('Error. You cannot do haploid and dopseudodiploid at the same time.')

# run the analysis
if onlyfit == False:
    calculate_allele_sharing(input_prefix = [genotypename, snpname, indivname], 
        output_prefix = output_prefix+'.out',
        target_popname = target_popname, 
        out_popname = out_popname, 
        minMAF = minMAF, 
        input_geno_is_diploid = input_geno_is_diploid, 
        stepD_cM = stepD_cM, 
        minD_cM = minD_cM, 
        maxD_cM = maxD_cM,
        input_distance_in_cM = input_distance_in_cM,
        chrom_to_analyze = cta)
else:
    print('Estimating the parameters only, with the weighted block jackknife procedure.')
    print('Input: '+output_prefix+'.out')

if out_popname == None:
    use_substracted = False
else:
    use_substracted = True

# run the expfit with jackknife procedure
print('================================================================================')
print('Running the exponential fitting with weighted jackknife')
    
with open(output_prefix+'.fit', 'w') as fout:
    fout.write('pop\tchromosome\tA\tt\tc\tNRMSD\tblocksize\n')
    
    s = io.StringIO(open(blockSizeFile).read().replace('\t', ' '))
    bs = np.genfromtxt(s, dtype=int, delimiter = ' ')
    if bs.ndim == 0:
        fout.close()
        sys.exit('\n*** Error. The block size file is empty.')
    elif bs.ndim == 1:
        bs = np.reshape(bs, (-1, 2))
        
    R = expfit_1D_jackknife(input = output_prefix+'.out',
            blocksizes = bs,
            min_D_cM = minD_cM,
            max_D_cM = maxD_cM,
            use_substracted = use_substracted)
    
    for i in range(R.shape[0]-3):
        fout.write(''+target_popname+'\t'+str(i+1)+'\t'+'\t'.join([str(x) for x in R[i,:].A1])+'\n')
    fout.write(''+target_popname+'\t'+'MEAN'+'\t'+'\t'.join([str(x) for x in R[R.shape[0]-3,:].A1])+'\n')
    fout.write(''+target_popname+'\t'+'JK.MEAN'+'\t'+'\t'.join([str(x) for x in R[R.shape[0]-2,:].A1])+'\n')
    fout.write(''+target_popname+'\t'+'JK.SD'+'\t'+'\t'.join([str(x) for x in R[R.shape[0]-1,:].A1])+'\n')

fout.close()

print('\nEnd\n')
