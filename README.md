# ASCEND
ASCEND (Allele Sharing Correlation for the Estimation of Nonequilibrium Demography) is a method to estimate the age and intensity of founder events (or bottlenecks) using population genotype data.

## Installation

ASCEND is a Python script and does not require installation.

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

To run an ASCEND analysis (minimal usage):

`python ASCorrelation_v5.2.1.py -f prefix_of_the_EIGENSTRAT_files -p label_of_the_population_to_analyze -o output_file -maxProportionNA maximum_proportion_of_missing_data`

By default, ASCEND assumes that the genetic positions are in centiMorgans and that the samples are diploid.

### Mandatory parameters

- `-f`
- `-p`
- `-o`
- `-maxProportionNA`

### Optional parameters

- `-r`
- `-minMAF`
- `-minD`
- `-maxD`
- `-stepD`
- `--haploid`
- `--pseudodiploidize`
- `--Morgans`
- `--chrom`

## Output

For each analysis, ASCEND outputs a single file with the following format:


## Example

An example run is provided in the repository `example`.

## Support
Send queries to Remi Tournebize (remi dot tournebize at gmail dot com) or Priya Moorjani (moorjani at berkeley dot edu).



