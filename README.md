# gollum
Detection of rings on acrocentric chromosomes using short reads sequencing

## Installation

This tool is available through the pip package manager in python (>=3.8) and can be installed with:  
`pip install gollumpy`.

`gollum` depends on 2 additional tools:
- `samtools` which can be installed by various means including `conda` or `brew` for mac
- `blat` which can be installed using `conda` or a binary can directly be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/

## Usage
`gollum` is designed to function on whole genome sequencing aligned using T2T reference genome.

```
gollum -h
usage: gollum [-h] [--chr CHR] [-s SAMPLE_NAME] -f FASTA -o OUTPUT_DIR [--blacklist-bed BLACKLIST_BED] [--e-value E_VALUE] [--bit-score BIT_SCORE] [--dbscan-eps DBSCAN_EPS] [--dbscan-min DBSCAN_MIN] cram

positional arguments:
  cram                  cram file to search for ring

optional arguments:
  -h, --help            show this help message and exit
  --chr CHR             chromosome for which to search for a ring (need to be acrocentric)
  -s SAMPLE_NAME, --sample-name SAMPLE_NAME
                        sample name, used for output files (default is the same as the cram name)
  -f FASTA, --fasta FASTA
                        path to the fasta used for the cram alignment
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        path to the output directory
  --blacklist-bed BLACKLIST_BED
                        bed of blacklisted regions
  --e-value E_VALUE     e-value threshold for blat (keeping <=)
  --bit-score BIT_SCORE
                        bit score threshold for blat (keeping >=)
  --dbscan-eps DBSCAN_EPS
                        EPS value used by DBSCAN inside when running the clustering of the positions (default: 100)
  --dbscan-min DBSCAN_MIN
                        Min number of reads inside a DBSCAN cluster (default: 5)
```

## Example

We are providing a simplified cram file in the test directory in this repository. In addition to this file, to run gollum, you will also need the T2T reference fasta file, which can be downloaded from https://github.com/marbl/CHM13 ([direct link to the file](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz)).

You will then be able to run the exemple using:  
`gollum C0011LO.disc.cram -s C0011LO -f chm13v2.0.fa -o output`

Expected output:
```
sample  cluster breakpoint      supporting_reads        acrocentric_specificity
C0011LO cluster_0       chr22:47097797  11      9.704158145875937
```

additionaly several files are produced (and are provided in this repository in the test/output sub directory):
- `C0011LO.log` a log file of the execution of `gollum`
- `C0011LO.supporting_reads.tsv` which list all the reads supporting the presence of a ring in this individual.
- and 2 temporary files are left `C0011LO.read2_for_blat.fasta` and `C0011LO.read2.blat.tsv` which are the entry and exit for the blat command.

## Remarques

`gollum` was developped and tested for rings on chromosome 22. It should work for rings on other acrocentric chromosomes, but was not tested in those cases.

We are providing in the test directory a `blacklist.bed` file. This is the regions we excluded from our analysis on chromosome 22 because they contained recuring events that are most likely due to something different than a ring. For instance, the beginning of the 22q contains large regions homologous to other regions from the 22p ([Guarracino et al 2023](https://doi.org/10.1038/s41586-023-05976-y)), additionaly highly recurring inversion in low complexity regions can also match to the 22p and create false positives. This file will most likely need to be expanded depending on the cohort analysed.
