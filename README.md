# ASCEND
ASCEND (Allele Sharing Correlation for the Estimation of Nonequilibrium Demography) is a method to estimate the age and intensity of founder events (or bottlenecks) using population genotype data.

# Installation

ASCEND is a Python script and does not require prior installation apart specific modules (numpy, scipy) that you can install using the command:

`pip3 install numpy`

# Requirements

In its current implementation, ASCEND requires that your population is comprised of at least 5 samples.

Since ASCEND relies on the recombination map, make sure your SNPs have the most accurate genetic position (see https://github.com/sunyatin/itara/blob/master/liftover.py to lift over a new genetic map using physical positions).

# Input

ASCEND requires that the input data is in EIGENSTRAT format (see https://reich.hms.harvard.edu/software/InputFileFormats). The EIGENSTRAT format is comprised of three files:

- `*.ind` A 3-column file with each individual on one row, and columns are (1) individual name, (2) individual sex (note that sex is not used for the ASCEND analysis), (3) population label
- `*.snp` A 6-column file with each SNP on one row, and columns are (1) SNP name, (2) chromosome, (3) genetic position (in Morgans or centiMorgans), (4) physical position (in base pairs), 5th and 6th columns are the ancestral/derived or reference/alternate alleles but these 2 columns are not taken into account for the ASCEND analysis
- `*.geno` The genotype matrix with no delimiter between genotypes, each row is one SNP and each column is one individual, genotypes are encoded as 0 (= 1/1), 1 (=0/1) or 0 (=0/0). Missing genotypes are encoded as 9.

You can convert your file into EIGENSTRAT using the CONVERTF program (see https://github.com/argriffing/eigensoft/tree/master/CONVERTF).

Note that the .geno file must not be compressed/binary.

# Command line

To run an ASCEND analysis:

`python3 ASCEND.py -p [NameOfTheParameterFile].par`

Note that by default, ASCEND assumes that the genetic positions are in Morgans and that the samples are diploid.

For reliable estimation, use a minimum of 5 diploids in the target and outgroup populations.

For best performance, we advise to use at most 20-30 individuals for the target and ~15 individuals for the outgroup population.

# Full list parameters

Note that you can comment any line and option using "#" (the software will ignore those lines). Also, the options can be written in any order.

*Mandatory options*

- `genotypename: STRING` name of the input .geno file
- `snpname: STRING` name of the input .snp file
- `indivname: STRING` name of the input .ind file
- `outputprefix: STRING` prefix of the output file, ASCEND will automatically append the appropriate extensions of the output files
- `targetpop: STRING` name of the target population to analyze

*Optional options (if not provided, ASCEND will take the default values)*

- `outpop: STRING` name of the outgroup population (if not provided, ASCEND will not compute the cross-population correlation and will output NaN in the corresponding column)

Related to genetic data

- `chrom: INT` add this option to restrict the analysis to a specific chromosome
- `haploid: NO` ASCEND assumes genotypes are diploid but if you set this option to YES it will interpret your genotypes as haploid (default: NO)
- `dopseudodiploid: YES` set YES if your genotypes have to be pseudodiploidized (i.e. for heterozygous genotypes, one allele will be randomly picked and set in two copies) (default: NO)

Related to SNP filtering

- `maxpropsharingmissing: 1.0` maximum proportion of missing allele sharing items allowed, above this threshold the SNP will be discarded (default: 1.0)
- `minmaf: 0` minimum Minor Allele Frequency that is allowed for a SNP to be taken into account, i.e. if a SNP has MAF<minMAF, it will be excluded (default: 0.0)

Related to the decay curve parameters

- `mindis: 0` minimum genetic distance in Morgans (default: 0.0 M)
- `maxdis: 0.3` maximum genetic distance in Morgans (default: 0.3 M)
- `binsize: 0.001` size of the genetic distance bin in Morgans (default: 0.001 M)
- `morgans: YES` set NO if your input genetic distances are in centiMorgans (by default ASCEND assumes Morgans) (default: YES)

Related to the algorithm

- `usefft: YES` whether to use the Mesh + Fast Fourier Transforms (FFT) algorithm which speeds up the calculation by up to 8,000 times with only marginal approximations, note that if you have less than 10,000 SNPs per chromosome, we would advice using the naive algorithm i.e. setting `usefft: NO` (default: YES)
- `qbins: 100` number of mesh points within each bins of the decay curve to consider for the mesh-FFT approximation (a higher number increases the mesh resolution and hence the accuracy of the decay curve) (default: 100)

Related to the fitting

- `onlyfit: NO` set YES if you want to do the estimation of the parameters directly, using files that have been already output by the script (default: NO)
- `blocksizename: STRING` add this option to indicate the name of a file containing the per-chromosome weights to use for the weighted jackknife analysis; the file has two tab-separated columns: (i) the chromosome label (should be the same as in the .snp file) and (ii) the number of SNPs on the chromosome or the chromosome length in bp; if this option is not provided, ASCEND will automatically calculate the chromosome weights as the number of SNPs in the input .snp file.

# Output

For each analysis (if `onlyfit: NO`), ASCEND will output two files with extensions `.out` (the decay curves) and `.fit` (the fits) respectively. If `onlyfit: YES`, ASCEND will output only the `.fit` file.

### `.out`

The `.out` file contains the calculated decay curves over each chromosome and has 7 columns:
- `chrom` the chromosome number
- `bin.left.bound` the left boundary of the genetic distance bins
- `bin.center` the center of the genetic distance bins
- `cor.pop` the within-population allele sharing correlation for the bin
- `cor.bg` the cross-population allele sharing correlation for the bin
- `cor.substracted` the within-population allele sharing correlation subtracted by the cross-population for the bin
- `n.pairs` the number of SNP pairs used in the calculation of the allele sharing correlation for the bin

Note that in case where you do not provide an outgroup population, the `cor.bg` and `cor.substracted` will be output with `nan` as values.

### `.fit`

The `.fit` file provides the estimates of the exponential model (that you can then use to estimate the founder age and the founder intensity) with their associated standard errors. To compute standard errors, the script performs a weighted block jackknife where blocks are the chromosomes and weights are their sizes or their number of SNPs (that are provided in the `blocksizename` file). We fit an exponential function of the form `z(d) = A exp(-2dt) + c` where `z(d)` is the allele sharing correlation at the genetic distance bin `d`, `A` is the amplitude and `t` is the rate of exponential decay. 

The `.fit` file has 7 columns:
- `pop` the target population label
- `chromosome` the chromosome number, or else `MEAN` for the parameter estimates when averaging the decay curves over all chromosomes, `JK.MEAN` for the jackknife mean estimate and `JK.SE` for the jackknife standard error
- `A` the amplitude of the exponential model
- `t` the exponential decay rate
- `c` the affine term of the exponential model
- `NRMSD` the root mean squared error of the exponential fit, normalized by the range of the fitted allele sharing correlation values (range = max - min)
- `blocksize` the chromosome weights

# Full example

An example run is provided in the repository `example`. You can re-run it using the command:

`python3 ASCEND_6.py -p example.par`

To plot the decay curves with their associated fits, you can use the R script:

`Rscript plot_ASCEND.R example.out example.fit example.png TRUE 0.2`

The parameters of the `plot_ASCEND.R` script are: [decay_curve.out] [fits.fit] [output_fig_name.png] [TRUE or FALSE if you want to subtract the within-correlation by the cross-correlation or not] [NRMSD threshold]

This script will output a plot of the allele sharing correlation decay curve (black points) along with the fitted exponential model (red line). In the top-right corner it provides the estimates of founder age (Tf) and intensity (If) with their associated 95% confidence intervals within brackets as well as the NRMSD. Note that the plot will be grayed if the NRMSD is below the user-provided threshold or if the fitting failed or if any of the CI95% includes the value 0.

The example provided is a simulation with 3 chromosomes of a founder event occurring 50 generations ago with intensity 10% so your estimates of Tf and If in the output plot should overlap with these numbers.

## Picking random samples as outgroups

If you want to pick `n` random samples (random sampling without replacement) as an outgroup population from your original dataset, first run the following script:

`python3 pickoutgroups.py -p NameOfTheParameterFile.par`

The parameter file takes 8 arguments and 1 optional (`seed`):

- `genotypename: STRING` the input .geno file
- `snpname: STRING` the input .snp file
- `indivname: STRING` the input .ind file
- `genooutfilename: STRING` the output .geno file
- `snpoutfilename: STRING` the output .snp file
- `indoutfilename: STRING` the output .ind file
- `outgroupsize: INT` number of outgroup individuals to sample (we recommand a size of 15 individuals)
- `targetpop: STRING` the label of the target population
- `seed: INT` a seed for the random sampling of outgroup individuals (if this option is not provided, will use the current timestamp as a seed)

The script will basically output the data subset to the target samples along with `outgroupsize` random individuals that have been set with the population label `OUTGROUP`.

Note that if you use ASCEND after `pickoutgroups.py`, make sure that the option `outpop` in ASCEND is set as:

`outpop: OUTGROUP`

### Full usage example

Example of a full run using this outgroup strategy:

`python3 pickoutgroups.py -p outgroup.par`

`python3 ASCEND_v6.py -p example_OUTGROUP.par`

`Rscript plot_ASCEND.R example_OUTGROUP.out example_OUTGROUP.fit example_OUTGROUP.png TRUE 0.2`

# Troubleshooting

### UnicodeDecodeError
If your input are in PACKED EIGENSTRAT format (i.e. the .geno file is binary), ASCEND will throw the error
`UnicodeDecodeError: 'utf-8' codec can't decode byte 0x86 in position 1936: invalid start byte`

To solve this problem, you will have to convert your input dataset to EIGENSTRAT using `convertf`, cf. https://github.com/DReichLab/AdmixTools/tree/master/convertf

# License

This software is licensed for academic and non-profit use only.

# Support
Send queries to Remi Tournebize (remi dot tournebize at gmail dot com) or Priya Moorjani (moorjani at berkeley dot edu).



