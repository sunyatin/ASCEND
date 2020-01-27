# ASCEND
ASCEND (Allele Sharing Correlation for the Estimation of Nonequilibrium Demography) is a method to estimate the age and intensity of founder events (or bottlenecks) using population genotype data.

## Installation

ASCEND is a Python script and does not require prior installation apart specific modules (numpy) that you can install using the command:

`pip3 install numpy`

## Requirements

In its current implementation, ASCEND requires that your population is comprised of at least 5 samples.

Since ASCEND relies on the recombination map, make sure your SNPs have the most accurate genetic position (see https://github.com/sunyatin/itara/blob/master/liftover.py to lift over a new genetic map using physical positions).

## Input

ASCEND requires that the input data is in EIGENSTRAT format (see https://reich.hms.harvard.edu/software/InputFileFormats). The EIGENSTRAT format is comprised of three files:

- `*.ind` A 3-column file with each individual on one row, and columns are (1) individual name, (2) individual sex (note that sex is not used for the ASCEND analysis), (3) population label
- `*.snp` A 6-column file with each SNP on one row, and columns are (1) SNP name, (2) chromosome, (3) genetic position (in Morgans or centiMorgans), (4) physical position (in base pairs), 5th and 6th columns are the ancestral/derived or reference/alternate alleles but these 2 columns are not taken into account for the ASCEND analysis
- `*.geno` The genotype matrix with no delimiter between genotypes, each row is one SNP and each column is one individual, genotypes are encoded as 0 (= 1/1), 1 (=0/1) or 0 (=0/0). Missing genotypes are encoded as 9.

You can convert your file into EIGENSTRAT using the CONVERTF program (see https://github.com/argriffing/eigensoft/tree/master/CONVERTF).

## Command line

To run an ASCEND analysis:

`python3 ASCEND.py [NameOfTheParameterFile]`

Note that by default, ASCEND assumes that the genetic positions are in centiMorgans and that the samples are diploid.

### Full list parameters

- `genotypename: STRING` name of the input geno file
- `snpname: STRING` name of the input snp file
- `indivname: STRING` name of the input ind file
- `outputname: STRING` name of the output file
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

## Output

For each analysis, ASCEND outputs a single file with 7 columns:
- `chrom` the chromosome numbers
- `bin.left.bound` the left boundary of the genetic distance bins
- `bin.center` the center of the genetic distance bins
- `cor.pop` the within-population allele sharing correlation for the bin
- `cor.bg` the cross-population allele sharing correlation for the bin
- `cor.substracted` the within-population allele sharing correlation substracted by the cross-population for the bin
- `n.pairs` the number of SNP pairs used in the calculation of the allele sharing correlation for the bin

Note that in case where you do not provide an outgroup population, the `cor.bg` and `cor.substracted` will be empty.

## Example

An example run is provided in the repository `example`. You can run it using the command:

`python3 ASCEND_5.3.py example.par`

## Support
Send queries to Remi Tournebize (remi dot tournebize at gmail dot com) or Priya Moorjani (moorjani at berkeley dot edu).



