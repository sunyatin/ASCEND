# ASCEND
ASCEND (Allele Sharing Correlation for the Estimation of Nonequilibrium Demography) is a method to estimate the age and intensity of founder events (or bottlenecks) using population genotype data.

# Installation

ASCEND is a Python script and does not require prior installation apart specific modules (numpy) that you can install using the command:

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

# Command line

To run an ASCEND analysis:

`python3 ASCEND_v6.py -p NameOfTheParameterFile.par`

Note that by default, ASCEND assumes that the genetic positions are in centiMorgans and that the samples are diploid.

# Full list parameters

- `genotypename: STRING` name of the input geno file
- `snpname: STRING` name of the input snp file
- `indivname: STRING` name of the input ind file
- `blocksizename: STRING` the name of a file containing two tab-separated columns: (i) the chromosome label (should be the same as in the .snp file) and (ii) the number of SNPs on the chromosome or the chromosome length in bp
- `outputprefix: STRING` prefix of the output file, ASCEND will automatically append the extension `.out`
- `targetpop: STRING` name of the target population to analyze
- `outpop: STRING` name of the outgroup population (if not provided, ASCEND will not compute the cross-population correlation)
- `maxpropmissing: 1` maximum proportion of missing data allowed in the allele sharing vectors to be considered in the calculation of the decay curve (default: 1.0)
- `minmaf: 0` minimum allele frequency, if a SNP has MAF>minMAF it is excluded (default: 0.0)
- `mindis: 0.1` minimum genetic distance in centiMorgans (default: 0.0 cM)
- `maxdis: 30.0` maximum genetic distance in centiMorgans (default: 30.0 cM)
- `binsize: 0.1` size of the genetic distance bin (default: 0.1 cM)
- `haploid: NO` set YES if your genotypes are haploid (default: NO)
- `dopseudodiploid: YES` set YES if your genotypes have to be pseudodiploidized (i.e. for heterozygous genotypes, one allele will be randomly picked and set as two copies) (default: NO)
- `morgans: NO` set YES if your input genetic distances are in Morgans (by default ASCEND assumes centiMorgans) (default: NO)
- `chrom: INT` add this option to restrict the analysis to a specific chromosome
- `onlyfit: NO` set YES if you want to do the estimation of the parameters directly, using a file that has been already output by the script, with name `outputname.out` (default: NO)

# Output

For each analysis (if `onlyfit: NO`), ASCEND will output two files with extensions `.out` (the decay curves) and `.fit` (the fits) respectively. If `onlyfit: YES`, ASCEND will output only the `.fit` file.

### `.out`

The `.out` file contains the calculated decay curves over each chromosome and has 7 columns:
- `chrom` the chromosome numbers
- `bin.left.bound` the left boundary of the genetic distance bins
- `bin.center` the center of the genetic distance bins
- `cor.pop` the within-population allele sharing correlation for the bin
- `cor.bg` the cross-population allele sharing correlation for the bin
- `cor.substracted` the within-population allele sharing correlation substracted by the cross-population for the bin
- `n.pairs` the number of SNP pairs used in the calculation of the allele sharing correlation for the bin

Note that in case where you do not provide an outgroup population, the `cor.bg` and `cor.substracted` will be empty.

### `.fit`

The `.fit` file provides the estimates of the exponential model (that you can then use to estimate the founder age and the founder intensity) with their associated standard errors. To compute standard errors, the script performs a weighted block jackknife where blocks are the chromosomes and weights are their sizes or their number of SNPs (that are provided in the `blocksizename` file). We fit an exponential function of the form `z(d) = A exp(-2dt) + c` where `z(d)` is the allele sharing correlation at the genetic distance bin `d`, `A` is the amplitude and `t` is the rate of exponential decay. 

The `.fit` file has 7 columns:
- `pop` the target population label
- `chromosome` the chromosome number, or else `MEAN` for the parameter estimates when averaging the decay curves over all chromosomes, `JK.MEAN` for the jackknife mean estimate and `JK.SE` for the jackknife standard error.
- `A` the amplitude of the exponential model
- `t` the exponential decay rate
- `c` the affine term of the exponential model
- `NRMSD` the root mean squared error of the exponential fit, standardized by the range of the fitted correlation values
- `blocksize` for each chromosome, their corresponding weights

# Full example

An example run is provided in the repository `example`. You can run it using the command:

`python3 ASCEND_6.py example.par`

To plot the decay curves with their associated fits, you can use the RScript:

`Rscript plot_ASCEND.R example.out example.fit example.png TRUE 0.2`

The parameters of the `plot_ASCEND.R` script are: [decay_curve.out] [fits.fit] [output_fig_name.png] [TRUE or FALSE if you want to subtract the within-correlation by the cross-correlation or not] [NRMSD threshold]

## Picking random samples as outgroups

If you want to pick `n` random samples (random sampling without replacement) as an outgroup population from your original dataset, first run the following script:

`python3 pickoutgroups.py -p NameOfTheParameterFile.par`

The parameter file takes 8 arguments:

- `genotypename:` input_data.geno
- `snpname:` input_data.snp
- `indivname:` input_data.ind
- `genooutfilename:` output_data.geno
- `snpoutfilename:` output_data.snp
- `indoutfilename:` output_data.ind
- `outgroupsize:` number_of_outgroup_samples (we recommand to use a size of 15 individuals)
- `targetpop:` target_population

The script will basically output the data subset to the target samples along with `number_of_outgroup_samples` random individuals that have been set with the label `OUTGROUP`.

### Full usage example

Example of a full run using this outgroup strategy:

`python3 pickoutgroups.py -p outgroup.par`

`python3 ASCEND_v6.py -p example_OUTGROUP.par`

`Rscript plot_ASCEND.R example_OUTGROUP.out example_OUTGROUP.fit example_OUTGROUP.png TRUE 0.2`

# Troubleshooting

### UnicodeDecodeError
If your input are in PACKED EIGENSTRAT format (i.e. the geno file is compressed as a binmary), ASCEND will output an error:
UnicodeDecodeError: 'utf-8' codec can't decode byte 0x86 in position 1936: invalid start byte

To solve this problem, convert your input dataset to EIGENSTRAT using `convertf`: https://github.com/DReichLab/AdmixTools/tree/master/convertf

# Support
Send queries to Remi Tournebize (remi dot tournebize at gmail dot com) or Priya Moorjani (moorjani at berkeley dot edu).



