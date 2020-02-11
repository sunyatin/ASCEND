# -*- coding: utf-8 -*-
"""
### VERSION v5.1

Created on Thu Feb 21 09:13:51 2019
@author: Windows
# 24 jun 2019: created a mean to compute the correlation instead of the covariance in allele sharing
# 25 jul 2019: debugged (before we added the SNP in the count if its cor value was Nan or Inf)
# 30 jul 2019: compute MAF across the populations
# 09 aug 2019: handle missing data
# 10 aug 2019: added the quick mode in case there is no missing data
# 10 aug 2019: the script has been carefully checked with test_profile_script.R with various maxPropNA and pNA: it works fine
# 11 aug 2019: significantly improved the speed of execution up to 3-fold and also better memory management than in v5, it has also been checked with test_profile_script.R and validated the test
"""

import numpy as np
import time
import argparse
import sys
import warnings

np.seterr(divide='ignore', invalid='ignore')

parser = argparse.ArgumentParser(description='Calcule allele sharing statistics.')
parser.add_argument('-f', '--filePrefix', type=str, required=True, help='Prefix of the files to analyze.')
parser.add_argument('-p', '--targetPopulation', type=str, required=True, help='Name of the target population to analyze.')
parser.add_argument('-o', '--outfilePrefix', type=str, required=True, help='Prefix of the output file.')
parser.add_argument('-r', '--referencePopulation', type=str, default=None, help='Name of reference population that will allow substracting background LD.')
parser.add_argument('-maxProportionNA', '--maxProportionNA', type=float, required=True, help='Maximum proportion of NA to accept.')
parser.add_argument('-minMAF', '--minMAF', type=float, default=-1, help='Minimum allele frequency.')
parser.add_argument('-minD', '--minD', type=float, default=0.01, help='Minimum genetic distance.')
parser.add_argument('-maxD', '--maxD', type=float, default=20.0, help='Maximum genetic distance.')
parser.add_argument('-stepD', '--stepD', type=float, default=0.1, help='Bin size.')
parser.add_argument('--haploidGeno', action='store_true', default=False, help='Add this option if geno is haploid (by default assumes diploid).')
parser.add_argument('--Morgans', action='store_true', default=False, help='Add this option if genetic distance are in Morgans (by default assumes centiMorgans).')

args = parser.parse_args()

"""
import os

input_prefix = 'for_Vysya'
target_popname = 'Vysya'
out_popname = 'ANI47'
output_prefix = 'handleNA_noNA'
minMAF = 0
input_geno_is_diploid = True
stepD_cM = 0.1
minD_cM = 0
maxD_cM = 30
input_distance_in_cM = False
max_proportion_NA = 0

"""
input_prefix = args.filePrefix
target_popname = args.targetPopulation
out_popname = args.referencePopulation
output_prefix = args.outfilePrefix
minMAF = args.minMAF
input_geno_is_diploid = not args.haploidGeno
stepD_cM = args.stepD
minD_cM = args.minD
maxD_cM = args.maxD
input_distance_in_cM = not args.Morgans
max_proportion_NA = args.maxProportionNA

print('\n\n')
print(input_prefix)
print(output_prefix)
print('Target pop: '+target_popname)
print('Out pop: '+str(out_popname))
print('minMAF: '+str(minMAF))
print('Assumes the genotypes are diploids')
print('D:'+str(minD_cM)+' '+str(maxD_cM)+' by '+str(stepD_cM))
print('Input distance in cM: '+str(input_distance_in_cM))
print('Maximum proportion of allele sharing values that can be missing: '+str(max_proportion_NA))
print('===============================================================')

###
#####
###



""" 
add exception if number of snps not number of rows in genotype matrix
"""


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


def calculate_allele_sharing(input_prefix, output_prefix,
                             target_popname, 
                             input_geno_is_diploid, 
                             out_popname = None, 
                             minMAF = 0.05,
                             stepD_cM = 0.5, minD_cM = 0, maxD_cM = None,
                             input_distance_in_cM = False):
    start = time.time()
    print('>   '+input_prefix)
    if out_popname is not None:
        print('Analyzing '+target_popname+' and using an outgroup population: '+out_popname)
    else:
        print('Analyzing only '+target_popname+', without outgroup.')
    corr_decimals = 1e-10
    
    D_FULL = np.genfromtxt(input_prefix+'.snp', dtype=float, usecols = (1,2))
    CHR = D_FULL[:,0].astype(int)
    D_FULL = D_FULL[:,1]
    if input_distance_in_cM is False:
        print('Converting Morgans input into centiMorgans.')
        D_FULL = 100.0 * D_FULL
    POP = np.genfromtxt(input_prefix+'.ind', dtype=str, usecols = 2)
    
    chr_values = np.unique(CHR).tolist()
    print('Chromosomes: '+' '.join([str(x) for x in chr_values]))
    
    RR, first = [], True
    for chrom in chr_values:
        foc = np.where(CHR == chrom)[0]

        r0, r1 = foc[0], foc[-1]  
        print('>> chrom:    '+str(chrom)+'   -~-   Range: '+str(r0)+' => '+str(r1))
        D = D_FULL[foc]
        del foc
        G = np.genfromtxt(input_prefix+'.geno', delimiter=[1]*len(POP), dtype=int,
                               skip_header = r0, max_rows = r1-r0+1)
    
        if target_popname not in POP:
            sys.exit("ERROR! Pop1 not in list!")

        there_is_NA_in_dataset = False

        if out_popname is not None:
            G_out = G[:,POP==out_popname]
            if np.any(G_out==9):
                there_is_NA_in_dataset = True
            nh_2 = G_out.shape[1]
            if out_popname not in POP:
                sys.exit("ERROR! PopOut not in list!")
    
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
                    print(str(v_i)+'%')
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
            
            with open(output_prefix+'.out', 'w') as fout:
                fout.write('chrom\tbin.left.bound\tbin.center\tcov.pop\tcov.bg\tcov.substracted\tn.pairs\n')
                for i in range(RR.shape[0]):
                    fout.write(str(int(RR[i,0]))+'\t'+'\t'.join([str(x) for x in RR[i,1:]])+'\n')
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
            if input_geno_is_diploid is False:
                G = G[:,0:int(nh_1/2)] + G[:,int(nh_1/2):nh_1]
                if out_popname is not None:
                    G_out = G_out[:,0:int(nh_2/2)] + G_out[:,int(nh_2/2):nh_2]
        
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
                    print(str(v_i)+'%')
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
            
            with open(output_prefix+'.out', 'w') as fout:
                fout.write('chrom\tbin.left.bound\tbin.center\tcov.pop\tcov.bg\tcov.substracted\tn.pairs\n')
                for i in range(RR.shape[0]):
                    fout.write(str(int(RR[i,0]))+'\t'+'\t'.join([str(x) for x in RR[i,1:]])+'\n')
            fout.close()
        
    print('It took:   '+str(round((time.time()-start)/60,2))+' min\n')

##
### script
##

calculate_allele_sharing(input_prefix = input_prefix, 
    output_prefix = output_prefix,
    target_popname = target_popname, 
    out_popname = out_popname, 
    minMAF = minMAF, 
    input_geno_is_diploid = input_geno_is_diploid, 
    stepD_cM = stepD_cM, 
    minD_cM = minD_cM, 
    maxD_cM = maxD_cM,
    input_distance_in_cM = input_distance_in_cM)
    

