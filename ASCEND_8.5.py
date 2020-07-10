# -*- coding: utf-8 -*-
"""
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
# 26 jan 2020: use an external file to get the parameters, major change compared to previous versions (which was reading params from command line)
# 27 jan 2020: added the expfit directly to the script and added the random outgroup selection too
# 29 jan 2020: v7.0: now export 6 different output files: .log .out .perchrom.outs .perjk.outs .perjk.fits .fit // changed their format also
# 29 jan 2020: v7.0: now take as chromosome weights for the jk: the number of SNPs, automatically, if option not provided
# 29 jan 2020: v7.0: output the correct NRMSD based on the jackknife mean estimates of A, t, c
# 29 jan 2020: v7.0: output the plot of NRMSD
# 29 jan 2020: v7.0: changed the binning of naive
# 30 jan 2020: v7.0: changed input of distances as Morgans by default, made all subsequent changes
# 30 jan 2020: v7.0: changed the name of the max prop of missing allele sharing values
# 30 jan 2020: v7.0: changed all "is" or "is not" to == or !=
# 30 jan 2020: v7.0: checked consistency with ASCEND 6.0 on one example: OK
# 30 jan 2020: v8.0: implemented the FFT and added two options: usefft, qbins
# 31 jan 2020: v8.0: implemented FFT with cross pop correlation
# 31 jan 2020: v8.1: checked that objects are copied / checked pb with references
# 03 feb 2020: v8.1.1: speeded up the assignation of the SNPs to the subbins (@fft_core) + added possibility to specify multiple chr to analyze + user can provide None as option argument and this will be converted to np.None
# 04 feb 2020: v8.1.2: added sys.exit in case where chr provided does not exist in the snp file
# 08 feb 2020: v8.2: added the option "randomhet" to randomly pick up 0 or 2 shared alleles in case of hetero ind pairs (instead of setting 1, by default)
#              v8.3: ignore this version which introduced a variance inflation correction factor (removed for v8.4)
# 22 apr 2020: v8.4: introduced a check on the geno file (check that values are only 0, 1, 2 or 9) + introduced "main" syntax + very minor reformatting in the print2 log
# 29 apr 2020: v8.5: implemented the random outgroup picking within ASCEND for easier use, added param 'outpopsize' and 'seed'
                     added also a seed for dopseudodiploid and for randomhet
                     added also the possibility to create the output dir if it does not exist
"""

import numpy as np
import time, sys, warnings, io, os, argparse, random
from scipy.optimize import curve_fit, OptimizeWarning
from scipy import signal
import matplotlib.pyplot as plt
from copy import deepcopy as cp

####################################################################################################################
####################################################################################################################
####################################################################################################################

def main():
    
    version = '8.5'

    # np.seterr(divide='ignore', invalid='ignore', all='ignore')
    warnings.simplefilter('error', OptimizeWarning)

    parser = argparse.ArgumentParser(description='ASCEND v'+str(version))
    parser.add_argument('-p', '--parfile', type=str, required=True, help='Name of the parameter file.')
    args = parser.parse_args()
    parfile = args.parfile
    
    precision_cM = 6
    
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
        'targetpop' : '', 
        'outpop' : None, 
        'outputprefix' : '', 
        'minmaf' : 0, 
        'haploid' : 'NO', 
        'dopseudodiploid' : 'NO', 
        'binsize' : 0.001, # in Morgans
        'mindis' : 0.001, # in Morgans
        'maxdis' : 0.3,  # in Morgans
        'morgans' : 'YES', 
        'maxpropsharingmissing' : 1,
        'chrom' : None,
        'blocksizename': None,
        'onlyfit' : 'NO',
        'usefft': 'YES',
        'qbins' : 100,
        'randomhet' : 'NO',
        'outpopsize' : None,
        'seed' : None
    }
    
    for param in params:
        if param.strip() == '\n' or param.startswith('#') == True:
            continue
        name, value = param.split(':')
        if name not in options.keys():
            sys.exit(name+' is not recognized as a valid option, please check the documentation')
        value = value.strip()
        if value == 'None':
            value = None
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
        stepD_cM = 100.0 * float(options['binsize'])
        
    if options['mindis']=='':
        sys.exit('mindis argument is missing or is empty')
    else:
        minD_cM = 100.0 * float(options['mindis'])
        
    if options['maxdis']=='':
        sys.exit('maxdis argument is missing or is empty')
    else:
        maxD_cM = 100.0 * float(options['maxdis'])
        
    if options['morgans']=='':
        sys.exit('morgans argument is missing or is empty')
    else:
        if options['morgans'] == 'YES':
            input_distance_in_cM = False
        elif options['morgans'] == 'NO':
            input_distance_in_cM = True
        else:
            sys.exit('morgans argument you provided is unrecognized, must be either YES or NO')
     
    if options['maxpropsharingmissing']=='':
        sys.exit('maxpropsharingmissing argument is missing or is empty')
    else:
        x = float(options['maxpropsharingmissing'])
        if x>1 or x<0:
            sys.exit('maxpropsharingmissing must be contained in the range 0-1')
        max_proportion_NA = x
         
    if options['chrom']=='':
        sys.exit('chrom argument is missing or is empty')
    else:
        if options['chrom'] != None:
            chrar = options['chrom']
            chrar = chrar.split(',')
            chrar = [int(x) for x in chrar]
            cta = chrar
        else:
            cta = None
            
    if options['blocksizename']=='':
        sys.exit('blocksizename argument is missing or is empty')
    else:
        if options['blocksizename'] != None:
            blockSizeFile = options['blocksizename']
            if not os.path.exists(blockSizeFile):
                sys.exit('\n*** Error. File does not exist: '+blockSizeFile+'\n')
        else:
            blockSizeFile = None
            
    if options['usefft']=='':
        sys.exit('usefft argument is missing or is empty')
    else:
        if options['usefft'] == 'YES':
            usefft = True
        elif options['usefft'] == 'NO':
            usefft = False
            
    if options['qbins']=='':
        sys.exit('qbins argument is missing or is empty')
    else:
        qbins = int(options['qbins'])
        
    if options['randomhet']=='':
        sys.exit('randomhet argument is missing or is empty')
    else:
        if options['randomhet'] == 'YES':
            randomhet = True
        elif options['randomhet'] == 'NO':
            randomhet = False
        else:
            sys.exit('randomhet argument you provided is unrecognized, must be either YES or NO')
            
    if options['outpopsize']=='':
        sys.exit('outpopsize argument is missing or is empty')
    else:
        outpopsize = options['outpopsize']
        if outpopsize != None:
            outpopsize = int(outpopsize)
        
    if options['seed']=='':
        sys.exit('seed argument is missing or is empty')
    else:
        seed = options['seed']
        if seed != None:
            seed = int(seed)
        
    ########################
    # preliminary checks
    
    if out_popname == 'RANDOM' and outpopsize == None:
        sys.exit('You defined a random sampling of the outgroup population with "outpop: RANDOM" but did not provide the size of this population with the "outpopsize" argument')
    
    if out_popname == 'RANDOM':
        indivname2 = output_prefix+'.RandomOutpop.ind'
    
    ########################
    
    outdir = os.path.dirname(os.path.normpath(output_prefix))
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    flog = open(output_prefix+'.log', 'w')
        
    print2(flog, '\n\n                ᗩ ᔕ ᑕ Eᑎ ᗪ                \n\n')
    print2(flog, 'Version '+str(version)+'\n')
    print2(flog, genotypename+' '+snpname+' '+indivname)
    print2(flog, 'Output file name: '+output_prefix)
    print2(flog, 'Target population: '+target_popname)
    if out_popname == 'RANDOM':
        print2(flog, 'Outgroup population: defined by picking '+str(outpopsize)+' random individuals, the new ind file is stored at: '+indivname2)
    else:
        print2(flog, 'Outgroup population: '+str(out_popname))
    if out_popname != 'RANDOM' and out_popname != None and outpopsize != None:
        print2(flog, 'Note that the "outpopsize" argument has no effect if you specify an outgroup population which is not "RANDOM"')
    print2(flog, 'minMAF: '+str(minMAF))
    if input_geno_is_diploid:
        print2(flog, 'Assumes the genotypes are diploids')
    else:
        print2(flog, 'Assumes the genotypes are haploids')
    if pseudodiploidize:
        print2(flog, 'Converting genotypes into pseudohomozygous')
    if seed != None:
        print2(flog, 'Using seed: '+str(seed))
    print2(flog, 'Distance bins: '+str(minD_cM)+' cM to '+str(maxD_cM)+' cM by steps of '+str(stepD_cM)+' cM')
    print2(flog, 'Input distance in cM: '+str(input_distance_in_cM))
    print2(flog, 'Maximum proportion of allele sharing values that can be missing: '+str(max_proportion_NA))
    if blockSizeFile == None:
        print2(flog, 'No block size file provided, will use the number of SNPs per chromosome as weights')
    else:
        print2(flog, 'Block size file: '+blockSizeFile)
    
    if usefft:
        print2(flog, '\nUsing FFT algorithm')
        print2(flog, 'qbins: '+str(qbins))
    else:
        print2(flog, '\nUsing naive algorithm')
        
    if randomhet:
        print2(flog, 'CAUTION! You set "randomhet: YES" which means that in case of a pair of individuals being heterozygous at a SNP, the number of alleles shared will be randomly set as 0 or 2 (instead of 1).')
    
    print2(flog, '============================================================================================')
    
    if input_geno_is_diploid == False and pseudodiploidize == True:
        sys.exit('Error. You cannot do haploid and dopseudodiploid at the same time.')
        
    if out_popname == 'RANDOM':
        print2(flog, 'Random outgroup definition:')
        print2(flog, 'Number of outgroup individuals: '+str(outpopsize))
        random_outgroup(indivname, indivname2, target_popname, outpopsize, seed)
        indivname = indivname2
        out_popname = 'OUTGROUP'
        print2(flog, '============================================================================================')
    
    # run the analysis
    if onlyfit == False:
        calculate_allele_sharing(input_prefix = [genotypename, snpname, indivname],
                                        usefft = usefft,
                                        qbins = qbins,
                                        output_prefix = output_prefix,
                                        target_popname = target_popname, 
                                        out_popname = out_popname, 
                                        minMAF = minMAF, 
                                        input_geno_is_diploid = input_geno_is_diploid, 
                                        pseudodiploidize = pseudodiploidize,
                                        flog = flog,
                                        max_proportion_NA = max_proportion_NA,
                                        precision_cM = precision_cM,
                                        stepD_cM = stepD_cM, 
                                        minD_cM = minD_cM, 
                                        maxD_cM = maxD_cM,
                                        input_distance_in_cM = input_distance_in_cM,
                                        chrom_to_analyze = cta,
                                        randomhet = randomhet,
                                        seed = seed)
    else:
        print2(flog, 'Estimating the parameters only, with the weighted block jackknife procedure.')
        print2(flog, 'Input: '+output_prefix+'.out')
    
    if out_popname == None:
        use_substracted = False
        print2(flog, 'Not substracting by cross-population allele sharing correlation')
    else:
        use_substracted = True
        print2(flog, 'Substracting by cross-population allele sharing correlation')
    
    # run the expfit with jackknife procedure
    print2(flog, '============================================================================================')
    print2(flog, 'Running the exponential fitting with weighted jackknife')
    
    if blockSizeFile == None:
        blockSizeFile = output_prefix+'.qweights'
    
    s = io.StringIO(open(blockSizeFile).read().replace('\t', ' '))
    bs = np.genfromtxt(s, dtype=int, delimiter = ' ')
    if bs.ndim == 0:
        sys.exit('\n*** Error. The block size file is empty.')
    elif bs.ndim == 1:
        bs = np.reshape(bs, (-1, 2))    
    
    R = expfit_1D_jackknife(input = output_prefix+'.perchrom.outs',
            blocksizes = bs,
            output_prefix = output_prefix,
            flog = flog,
            min_D_cM = minD_cM,
            max_D_cM = maxD_cM,
            use_substracted = use_substracted)
    
    NRMSD = R[1]
    R = R[0]
    
    with open(output_prefix+'.perjk.fits', 'w') as fout:
        fout.write('run\tA\tt\tc\tblockweights\n')
        for i in range(R.shape[0]-3):
            fout.write(str(i+1)+'\t'+'\t'.join([str(x) for x in R[i,0:3].A1])+'\t'+str(R[i,4])+'\n')
    
    with open(output_prefix+'.fit', 'w') as fout:
        fout.write('param\tmean\tjk.mean\tjk.se\n')
        fout.write('A\t'+'\t'.join([str(x) for x in R[(R.shape[0]-3):(R.shape[0]-0),0].A1])+'\n')
        fout.write('t\t'+'\t'.join([str(x) for x in R[(R.shape[0]-3):(R.shape[0]-0),1].A1])+'\n')
        fout.write('c\t'+'\t'.join([str(x) for x in R[(R.shape[0]-3):(R.shape[0]-0),2].A1])+'\n')
        fout.write('NRMSD\t'+str(NRMSD[0])+'\t'+str(NRMSD[1])+'\tNA\n')
        
    # plot function    
    print2(flog, '\nEnd\n')
    
    flog.close()
    
    # plot the curve with associated fit
    
    R = np.genfromtxt(output_prefix+'.out', delimiter='\t', dtype='float', skip_header=1)
    FIT = np.genfromtxt(output_prefix+'.fit', delimiter='\t', dtype='float', skip_header=1)
    FIT = FIT[:,1:]
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    if use_substracted == True:
        plt.plot(R[:,0], R[:,3], 'o')
    else:
        plt.plot(R[:,0], R[:,1], 'o')
        
    yfit = FIT[0,1] * np.exp(-2 * FIT[1,1] * R[:,0]) + FIT[2,1]
    
    plt.plot(R[:,0], yfit, 'r-')
    
    title = 'target: '+target_popname+' - outgroup: '+str(out_popname)
    if use_substracted == True:
        title += '\nsubtracted by cross-pop correlation'
    else:
        title += '\nnot subtracted by cross-pop correlation'
    if pseudodiploidize == True:
        title += '\ngenotypes converted to pseudohomozygous'
    
    plt.title(title)
    plt.xlabel('Genetic distance (cM)')
    plt.ylabel('Allele sharing correlation')
    
    e1 = np.exp(1)*100
    
    textstr = '\n'.join((
        r'$T_f=%.0f\ [%.0f;%.0f]\ gBP$' % (FIT[1,1]*100, 100*(FIT[1,1]-1.96*FIT[1,2]), 100*(FIT[1,1]+1.96*FIT[1,2])),
        r'$I_f=%.1f\%%\ [%.1f;%.1f]$' % (FIT[0,1]*e1, e1*(FIT[0,1]-1.96*FIT[0,2]), e1*(FIT[0,1]+1.96*FIT[0,2])),
        r'$NRMSD=%.3f$' % (FIT[3,1] )))
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.6)
    
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            horizontalalignment='right', verticalalignment='top', bbox=props)
    
    fig.savefig(output_prefix+'.png', dpi=110)
    


####################################################################################################################
####################################################################################################################
####################################################################################################################

def print2(flog, text):
    print(text)
    flog.write(text+'\n')

####################################################################################################################
def random_outgroup(indivname, indivname2, targetpop, outgroupsize, seed):
    
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
    
    with open(indivname2, 'w') as fout:
        for i in range(IND.shape[0]):
            fout.write(IND[i,0]+' '+IND[i,1]+' '+IND[i,2]+'\n')
    fout.close()
    
    del IND

####################################################################################################################
def number_shared_alleles(pop1, pop2, randomhet = False, seed = None):
    dif = 2 - abs(pop1 - pop2)
    prod = pop1 * pop2
    if randomhet == True:
        if seed != None:
            random.seed(seed+1)
        dif[prod==1] = np.random.choice([0,2], size = sum(prod==1), replace = True)
    else:
        dif[prod==1] = 1
    dif[prod==81] = -9
    # the missing data will be encoded as negative values (-5, -6, -7 or -9)
    return dif.astype(np.int8)

####################################################################################################################
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

####################################################################################################################
# for FFT
def standardize(M, isDefined, nDefined):
    M = cp(M)
    isDefined = cp(isDefined)
    nDefined = cp(nDefined)
    
    M[isDefined==0] = 0
    means = np.sum(M, axis = 1) / nDefined
    M = np.transpose(np.transpose(M) - means).astype('float')
    M[isDefined==0] = 0
    SD = np.sqrt((M**2).sum(axis = 1) / nDefined)

    isDefined[SD==0] = 0
    nDefined = np.sum(isDefined, axis = 1)
    
    SD[SD==0] = 1 # this replacement will not bias anything because the values were excluded
    M = np.transpose(np.transpose(M) / SD)
    M[isDefined==0] = 0
    
    return [M, isDefined, nDefined]

####################################################################################################################
# for FFT
def correlation(M, INDEX, isDefined, bins_out, sub_bins):
    
    # sum over individual pairs
    for i in range(M.shape[1]):
        
        x = cp(M[INDEX,i])
        x[INDEX==-1] = 0
        n = cp(isDefined[INDEX,i])
        n[INDEX==-1] = 0
        
        x = x.astype(np.float)
        n = n.astype(np.float)
        
        c = signal.correlate(x, x, mode = 'full', method = 'fft')
        npairs = signal.correlate(n, n, mode = 'full', method = 'fft')
        
        c = c[c.size//2:]
        npairs = npairs[npairs.size//2:]
        npairs = npairs.astype(np.uint64)
        
        W = np.array([0] * len(bins_out)).astype(float)
        N = np.array([0] * len(bins_out))
        for g in range(len(bins_out)):
            if g == (len(bins_out) - 1):
                dstep = bins_out[1] - bins_out[0]
                idx = np.where( (sub_bins>=bins_out[g]) & (sub_bins<(bins_out[g]+dstep)) )[0]
            else:
                idx = np.where( (sub_bins>=bins_out[g]) & (sub_bins<bins_out[g+1]) )[0]
            W[g] = np.sum( c[idx] )
            N[g] = np.sum( npairs[idx] )
            #W[g] = np.sum( c[(g*qbins):((g+1)*qbins)] )
            #N[g] = np.sum( npairs[(g*qbins):((g+1)*qbins)] )
        
        if i == 0:
            C = W
            NN = N
        else:
            C = C + W
            NN = NN + N

    C = C / NN
    return [C, NN]

####################################################################################################################
# for FFT
def fft_core(ASD, ASD_noNA, D, max_proportion_NA, stepD_cM, qbins, bins_left_bound):
    isDefined = cp(ASD_noNA)
    nDefined = np.sum(isDefined, axis = 1)
    
    fft_goods = 1. - (1.*nDefined/ASD_noNA.shape[1]) <  max_proportion_NA
    isDefined = isDefined[fft_goods,:]
    nDefined = nDefined[fft_goods]
    
    subbins_max = max(D[fft_goods])
    
    subbins = np.arange(0, subbins_max, step = stepD_cM / qbins)
    pos_bin = np.digitize(D[fft_goods], subbins, right=False) - 1

    # assign unique SNPs to subbins
    uniq = np.unique(pos_bin, return_index = True)
    pos_bin = uniq[0]
    pos_bin_idx = uniq[1]
    INDEX = np.full(len(subbins), -1, dtype=np.int32)
    np.put(INDEX, pos_bin, pos_bin_idx)

    asd = standardize(ASD[fft_goods,:], isDefined, nDefined)
    
    isDefined = asd[1]
    nDefined = asd[2]
    asd = asd[0]
    
    BINS_R = correlation(asd, INDEX, isDefined, bins_left_bound, subbins)
                
    return(BINS_R)

####################################################################################################################
def calculate_allele_sharing(input_prefix, 
                             usefft,
                             qbins,
                             output_prefix,
                             target_popname, 
                             input_geno_is_diploid, 
                             pseudodiploidize,
                             flog,
                             max_proportion_NA = 1.0,
                             precision_cM = 6,
                             out_popname = None, 
                             minMAF = 0,
                             stepD_cM = 0.1, 
                             minD_cM = 0, 
                             maxD_cM = 30.0,
                             input_distance_in_cM = False,
                             chrom_to_analyze = None,
                             randomhet = False,
                             seed = False):
    
    start = time.time()
    if out_popname != None:
        print2(flog, 'Analyzing '+target_popname+' and using an outgroup population: '+out_popname)
    else:
        print2(flog, 'Analyzing only '+target_popname+', without outgroup')
    
    D_FULL = np.genfromtxt(input_prefix[1], dtype=float, usecols = (1,2))
    CHR = D_FULL[:,0].astype(int)
    D_FULL = D_FULL[:,1]
    if input_distance_in_cM == False:
        print2(flog, 'Converting Morgans input into centiMorgans')
        D_FULL = 100.0 * D_FULL
    POP = np.genfromtxt(input_prefix[2], dtype=str, usecols = 2)
    
    if out_popname == None:
        print2(flog, 'There are '+str(np.sum(POP==target_popname))+' target samples\n')
    else:
        print2(flog, 'There are '+str(np.sum(POP==target_popname))+' target and '+str(np.sum(POP==out_popname))+' outgroup samples\n')
    
    if np.sum(POP==target_popname) == 0:
        sys.exit('We found 0 target samples, there must a problem! Have you provided the right .ind file?')
    
    if chrom_to_analyze == None:
        chr_values = np.unique(CHR).tolist()
    else:
        chr_values = np.array(chrom_to_analyze)
    print2(flog, 'Chromosomes: '+' '.join([str(x) for x in chr_values])+'\n')
    
    # v8.4: check that the geno file contains only values in [0, 1, 2, 9]
    G = np.genfromtxt(input_prefix[0], delimiter=[1]*len(POP), dtype=np.int8, max_rows = 1000)
    G = np.unique(G)
    if np.all(np.isin(G, np.array([0,1,2,9]))) == False:
        sys.exit('Error. The geno file contains invalid genotype values (i.e. neither 0, 1, 2 nor 9)')
    del G
    
    RR, first, nSNP_per_chrom = [], True, []
    for chrom in chr_values:
        foc = np.where(CHR == chrom)[0]
        if len(foc) == 0:
            sys.exit('\nError. Chromosome '+str(chrom)+' does not seem to exist in your .snp file!')
        nSNP_per_chrom += [[chrom, len(foc)]]

        r0, r1 = foc[0], foc[-1]  
        print2(flog, '\n\n>> Chrom:    '+str(chrom)+'   -~-   Range: ['+str(r0+1)+'; '+str(r1+1)+']')
        D = D_FULL[foc]
        del foc
        G = np.genfromtxt(input_prefix[0], delimiter=[1]*len(POP), dtype=np.int8,
                               skip_header = r0, max_rows = r1-r0+1)
        
        if input_geno_is_diploid == False:
            G[G==1] = 2
            
        if pseudodiploidize == True:
            if seed != None:
                random.seed(seed+2)
            G[G==1] = np.random.choice([0,2], size = len(G[G==1]), replace = True)
    
        if target_popname not in POP:
            sys.exit("ERROR! The target population you specified is not present in the ind file")

        if out_popname != None:
            G_out = cp(G[:,POP==out_popname])
            nh_2 = G_out.shape[1]
            if out_popname not in POP:
                sys.exit("ERROR! The outgroup population you specified is not present in the ind file")
    
        G = G[:,POP==target_popname]
    
        # filter on the MAF
        nh_1 = G.shape[1]*2
        G_no9 = cp(G)
        G_no9[G_no9<=2] = 10
        G_no9 -= 9
        NH1 = G_no9.sum(1) * 2
        G0 = cp(G)
        G0 = G0 * G_no9
        del G_no9
        
        if out_popname != None:
            nh_2 = G_out.shape[1]*2
            
            if usefft == True:
                # maf
                bads = np.where(NH1==0)
                NH1[NH1==0] = 1
                maf = 1.0 * G0.sum(1) / NH1
                maf[maf>0.5] = 1 - maf[maf>0.5]
                maf[bads] = -1. * np.inf
                goods = np.where(maf>minMAF)
                del maf, bads
            else:
                G_out_no9 = cp(G_out)
                G_out_no9[G_out_no9<=2] = 10
                G_out_no9 -= 9
                NH2 = G_out_no9.sum(1) * 2
                G0_out = cp(G_out)
                G0_out = G0_out * G_out_no9
                # maf
                NH = NH1 + NH2
                bads = np.where(NH==0)
                NH[NH==0] = 1
                maf = 1.0 * (G0.sum(1)+G0_out.sum(1)) / NH
                maf[maf>0.5] = 1 - maf[maf>0.5]
                maf[bads] = -1. * np.inf
                goods = np.where(maf>minMAF)
                del maf, G0_out, G_out_no9, NH2, NH, bads
            
        else:
            # maf
            bads = np.where(NH1==0)
            NH1[NH1==0] = 1
            maf = 1.0 * G0.sum(1) / NH1
            maf[maf>0.5] = 1 - maf[maf>0.5]
            maf[bads] = -1. * np.inf
            goods = np.where(maf>minMAF)
            del maf, bads

        del G0, NH1

        nSNP_raw = G.shape[0]
        D = D[goods]
        G = G[goods,:][0]

        if out_popname != None:
            G_out = G_out[goods,:][0]     
    
        if maxD_cM == None:
            maxD_cM = D[len(D)]
        
        nSNP = G.shape[0]
        print2(flog, 'Raw nSNP: '+str(nSNP_raw)+'   MAF-filtered nSNP: '+str(nSNP))
                
        # allele sharing calculation   
        ASD1 = np.zeros((nSNP, int( (nh_1/2)*(nh_1/2-1)/2) ))
        ASD1 = ASD1.astype(np.int8)
        i = 0
        for j in range(0, G.shape[1]):
            for k in range(j+1, G.shape[1]):
                ASD1[:,i] = number_shared_alleles(G[:,j], G[:,k], randomhet, seed)
                i += 1
        ASD1_noNA = np.clip(ASD1, -1, 0) + 1

        if out_popname != None:
            ASD2 = np.zeros((nSNP, int((nh_2/2)*(nh_1/2)) ))
            ASD2 = ASD2.astype(np.int8)
            i = 0
            for j in range(0, G.shape[1]):
                for k in range(0, G_out.shape[1]):
                    ASD2[:,i] = number_shared_alleles(G[:,j], G_out[:,k], randomhet, seed)
                    i += 1
            ASD2_noNA = np.clip(ASD2, -1, 0) + 1

        del G
        if out_popname != None:
            del G_out
            
        bins_left_bound = np.arange(minD_cM, maxD_cM, step=stepD_cM)
        n_bins = len(bins_left_bound)
        #print2(flog, 'There will be '+str(n_bins)+' bins.')
        
        ########################################################
        ########################################################
        ############          NAIVE APPROACH        ############
        ########################################################
        ########################################################
        
        if usefft == False:
            
            BINS_R = np.asarray([0]*n_bins).astype(np.float)
            BINS_Rbg = np.asarray([0]*n_bins).astype(np.float)
            BINS_Rcor = np.asarray([0]*n_bins).astype(np.float)
            BINS_N = np.asarray([0]*n_bins).astype(np.int64)
            v_last = 0
            for i in range(nSNP):

                dd = D[(i+1):nSNP]
                d = D[i]
                
                d_ = abs(dd - d)
                d_ = np.round_(d_, precision_cM)
                indices = (d_ - minD_cM) / stepD_cM
                indices = np.round_(indices, precision_cM)
                indices = (np.floor(indices)).astype(np.int32)
                del d, dd, d_
                
                goods_Dref = np.where((0<=indices) & (indices<=(n_bins-1)))[0]
                indices = indices[goods_Dref]
                
                (r_, inf1) = compute_correlation(ASD1, i, ASD1_noNA, goods_Dref, max_proportion_NA)
    
                if out_popname != None:
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
                    if out_popname == None:
                        BINS_Rbg = [np.nan]*len(BINS_R)
                        BINS_Rcor = [np.nan]*len(BINS_R)
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
                
                
        ########################################################
        ########################################################
        ############           FFT APPROACH         ############
        ########################################################
        ########################################################
        
        else:
        
            X = fft_core(ASD1, ASD1_noNA, D, max_proportion_NA, stepD_cM, qbins, bins_left_bound)
            BINS_R, BINS_N = X[0], X[1]

            if out_popname != None:
                #print('.....cross-pop correlation')
                X = fft_core(ASD2, ASD2_noNA, D, max_proportion_NA, stepD_cM, qbins, bins_left_bound)
                BINS_Rbg = X[0]
                BINS_Rcor = BINS_R - BINS_Rbg
            else:
                BINS_Rbg = BINS_R.copy() * np.nan
                BINS_Rcor = BINS_R.copy() * np.nan

            R = np.column_stack(([chrom]*len(BINS_R),
                                   bins_left_bound,
                                   bins_left_bound + stepD_cM/float(2),
                                   BINS_R * BINS_N,
                                   BINS_Rbg * BINS_N,
                                   BINS_Rcor * BINS_N,
                                   BINS_N))
            
            if first:
                RR, first = R, False
            else:
                RR = np.concatenate((RR, R), axis=0) 

        ########################################################
        ####################         export to file      #######
        ########################################################
        with open(output_prefix+'.perchrom.outs', 'w') as fout:
            fout.write('chrom\tbin.left.bound\tsum.cor.pop\tsum.cor.bg\tsum.cor.subtracted\tn.pairs\n')
            for i in range(RR.shape[0]):
                fout.write(str(int(RR[i,0]))+'\t')
                fout.write('{:.5f}'.format(float(RR[i,1]))+'\t')
                fout.write('\t'.join([str('{:.5f}'.format(float(x))) for x in RR[i,3:6]])+'\t')
                fout.write(str(int(RR[i,6]))+'\n')

    print2(flog, '\n\nIt took:   '+str(round((time.time()-start)/60,2))+' min\n')

    with open(output_prefix+'.qweights', 'w') as fout:
        for j in range(len(nSNP_per_chrom)):
            fout.write(str(nSNP_per_chrom[j][0])+'\t'+str(nSNP_per_chrom[j][1])+'\n')
            
            
####################################################################################################################
def exp(x, a, b, c):
        return a + b*np.exp(c*x)
    
####################################################################################################################
def exp_final(x, A, t, c):
        return A*np.exp(-2*t*x)+c

####################################################################################################################
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


####################################################################################################################
def expfit_1D_jackknife(input, blocksizes, output_prefix, flog, min_D_cM = 0, max_D_cM = 30.0, use_substracted = False):
    
    issue = False
    
    F = np.genfromtxt(input, delimiter='\t', dtype='float', skip_header=1)
    F_q = np.asarray([int(x) for x in F[:,0]])
    chrom = [int(x) for x in np.unique(F[:,0]).tolist()]

    intersection = np.intersect1d(np.array(chrom), blocksizes[:,0]).tolist()
    
    print2(flog, '\tChromosomes in input file:\n'+' '.join([str(int(x)) for x in np.unique(F[:,0]).tolist()]))
    print2(flog, '\tChromosomes at the intersection of the two files:\n'+' '.join([str(x) for x in intersection]))
    
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
            XY[:,2:6] += xy[:,2:6]

    with open(output_prefix+'.out', 'w') as fout:
        fout.write('bin.left.bound\tmean.cor.pop\tmean.cor.bg\tmean.cor.subtracted\tn.pairs\n')
        for i in range(XY.shape[0]):
            fout.write('{:.5f}'.format(XY[i,1])+'\t')
            fout.write('{:.5f}'.format(np.divide(XY[i,2], XY[i,5]))+'\t')
            fout.write('{:.5f}'.format(np.divide(XY[i,3], XY[i,5]))+'\t')
            fout.write('{:.5f}'.format(np.divide(XY[i,4], XY[i,5]))+'\t')
            fout.write(str(int(XY[i,5]))+'\n')

    if use_substracted:
        XY = np.column_stack((XY[:,1], np.divide(XY[:,4], XY[:,5])))
    else:
        XY = np.column_stack((XY[:,1], np.divide(XY[:,2], XY[:,5])))
        
    GLOBAL_CURVE = XY.copy()

    params = jacquelin_exp1d(XY, min_D_cM, max_D_cM)

    if type(params) == type(None) and params == None:
        general_mean = np.array([None]*4)
    else:
        general_mean = params

    # per chromosome estimate
    PARAMS = []
    with open(output_prefix+'.perjk.outs', 'w') as fout:
        
        fout.write('run\tbin.left.bound\tmean.cor.pop\tmean.cor.bg\tmean.cor.subtracted\tn.pairs\n')
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
                    XY[:,2:6] += xy[:,2:6]

            for dd in range(XY.shape[0]):
                fout.write(str(i+1)+'\t')
                fout.write('{:.5f}'.format(XY[dd,1])+'\t')
                fout.write('{:.5f}'.format(np.divide(XY[dd,2], XY[dd,5]))+'\t')
                fout.write('{:.5f}'.format(np.divide(XY[dd,3], XY[dd,5]))+'\t')
                fout.write('{:.5f}'.format(np.divide(XY[dd,4], XY[dd,5]))+'\t')
                fout.write(str(int(XY[dd,5]))+'\n')
    
            if use_substracted:
                XY = np.column_stack((XY[:,1], np.divide(XY[:,4], XY[:,5])))
            else:
                XY = np.column_stack((XY[:,1], np.divide(XY[:,2], XY[:,5])))
            
            params = jacquelin_exp1d(XY, min_D_cM, max_D_cM)
            
            if type(params) == type(None) and params == None:
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
    
    if issue == False:
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
        
        # compute the global NRMSD
        x = GLOBAL_CURVE[:,0]
        y = GLOBAL_CURVE[:,1]
        yfit = exp_final(x, general_mean[0], general_mean[1], general_mean[2])
        NRMSD = 1./(np.max(yfit)-np.min(yfit)) * np.sqrt(np.mean( (yfit-y)**2 ))
        yfit = exp_final(x, MEAN[0], MEAN[1], MEAN[2])
        NRMSD2 = 1./(np.max(yfit)-np.min(yfit)) * np.sqrt(np.mean( (yfit-y)**2 ))
        NRMSD = [NRMSD, NRMSD2]
        
        # append
        PARAMS = np.vstack((PARAMS, general_mean, MEAN, SD))
        PARAMS = np.hstack((PARAMS, np.asmatrix(np.append(N, (None, None, None))).T ))
    else:
        PARAMS = np.vstack((PARAMS, general_mean, [None]*4, [None]*4))
        PARAMS = np.hstack((PARAMS, np.asmatrix(np.append(N, (None, None, None))).T ))
        NRMSD = [None, None]
    
    return [PARAMS, NRMSD]


####################################################################################################################
if __name__=="__main__":
    main()

