This is the repository of the paper entitled
"the social environment has little impact on inbreeding depression in a social mammal"
by King To Chan, Alexandre Courtiol, Leonie F. Walter, Oliver P. Honer, Larissa S. Arantes, Eve Davidian, Philemon Naman, Loeske E.B. Kruuk, Josephine M. Pemberton & Kasha Strickland.

## Preparing your R session

To run the script in this repository, you need to install several R packages:

```r
std_pkgs_needed <- c("tidyverse", "patchwork", "ggExtra",            # for wrangling and plotting
                     "kinship2", "ggroups", "nadiv", "sequoia",      # for handling relatedness
                     "spaMM", "doSNOW", "mgcv")                      # for modelling

install.packages(std_pkgs_needed)
```

## Reproducing the prepared dataset

To reproduce the prepared dataset, you need:

- the full hyena database (not provided, but not needed to reproduce the results, see below)

- the variant call format (VCF) file storing the SNP-information (provided at https://zenodo.org/records/19709944)


## Reproducing the results

To reproduce the results of this paper, run the script named `00_master_script.R`.

All other scripts are called from within this main script.

Note that you do not need the full database to be able to reproduce the results since those are derived from the prepared dataset provided in the folder data.