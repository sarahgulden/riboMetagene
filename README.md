
<!-- README.md is generated from README.Rmd. Please edit that file -->

# riboMetagene

riboMetagene extracts 5’ end counts from Ribo-seq BAM files, shows the
read length distribution, and determines the P-site offset per read
length (work in progress).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("sarahgulden/riboMetagene")
```

## Workflow

``` r
library(riboMetagene)

samples <- list_bam_samples("path/to/bams")
counts  <- extract_5prime(samples, out_dir = "5prime_tables")
qc      <- readlength_distribution(counts, samples)
```
