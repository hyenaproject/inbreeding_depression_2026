This is the repository of the paper entitled
***The social environment has little impact on inbreeding depression in a social mammal***
by King To Chan, Alexandre Courtiol, Leonie F. Walter, Oliver P. Honer, Larissa S. Arantes, Eve Davidian, Philemon Naman, Loeske E.B. Kruuk, Josephine M. Pemberton & Kasha Strickland.

### Workflow instruction

To reproduce the results of this paper, you need to **clone this repository** on your computer.
You can do this by downloading the repository as a zip file (green button above) and unzipping it on your computer.
Alternatively you can use the git interface of your choice (e.g. git, GitHub desktop, RStudio).

Once this is done, you can run the entire analyses on pre-stored data.
For this, simply **open the RStudio project in RStudio** (contained in the unzipped folder) in RStudio, **select the script `00_main_script.R`, and run it line by line**.
This script will call other scripts (in the correct order) to reproduce the results of the paper.
For details, do inspect those other scripts too.

If you are not using RStudio, just make sure to set your working directory to the folder where you unzipped the repository, and run the script `00_main_script.R`.

#### Notes

The extracted dataset derived from the full hyena dataset needed for all analyses is already contained in the project folder *data*. 

In case you want to reproduce the SNP filtering procedure, you need to download [plink](https://www.cog-genomics.org/plink/) (we used the version 1.9).
Depending on your system, simply storing this program in the project folder may work for R to find it, otherwise add `plink` to the system PATH (search 'add program to system path' or similar on the web if you don't know how to).

The data frames produced by the SNP filtering are already pre-stored in the folder *data*, in case you want to skip this step. 
 
All other scripts are called from within the main script `00_main_script.R`.

The main script and the auxiliary scripts called by it will make sure to download the R packages you need and the variant call format (VCF) file storing the SNP-information (archived at https://zenodo.org/records/19709944)

#### Devel info

```r
> sessioninfo::package_info()
 package      * version    date (UTC) lib source
 backports      1.5.1      2026-04-03 [2] CRAN (R 4.5.3)
 boot           1.3-32     2025-08-29 [2] CRAN (R 4.5.1)
 cachem         1.1.0      2024-05-16 [2] CRAN (R 4.5.1)
 checkmate      2.3.4      2026-02-03 [2] CRAN (R 4.5.2)
 cli            3.6.6      2026-04-09 [2] CRAN (R 4.5.3)
 codetools      0.2-20     2024-03-31 [2] CRAN (R 4.5.1)
 CoprManager    0.5.8      2026-01-29 [4] local
 devtools       2.5.1      2026-04-16 [2] CRAN (R 4.5.3)
 digest         0.6.39     2025-11-19 [2] CRAN (R 4.5.2)
 doSNOW       * 1.0.20     2022-02-04 [2] CRAN (R 4.5.1)
 dplyr        * 1.2.1      2026-04-03 [2] CRAN (R 4.5.3)
 ellipsis       0.3.3      2026-04-04 [2] CRAN (R 4.5.3)
 evaluate       1.0.5      2025-08-27 [2] CRAN (R 4.5.1)
 farver         2.1.2      2024-05-13 [2] CRAN (R 4.5.1)
 fastmap        1.2.0      2024-05-15 [2] CRAN (R 4.5.1)
 forcats      * 1.0.1      2025-09-25 [2] CRAN (R 4.5.1)
 foreach      * 1.5.2      2022-02-02 [2] CRAN (R 4.5.1)
 fs             2.1.0      2026-04-18 [2] CRAN (R 4.5.3)
 generics       0.1.4      2025-05-09 [2] CRAN (R 4.5.1)
 ggExtra      * 0.11.0     2025-09-01 [2] CRAN (R 4.5.1)
 ggplot2      * 4.0.3      2026-04-22 [2] CRAN (R 4.5.3)
 ggroups      * 2.1.2      2022-03-27 [2] CRAN (R 4.5.1)
 glue           1.8.1      2026-04-17 [2] CRAN (R 4.5.3)
 gtable         0.3.6      2024-10-25 [2] CRAN (R 4.5.1)
 hms            1.1.4      2025-10-17 [2] CRAN (R 4.5.1)
 htmltools      0.5.9      2025-12-04 [2] CRAN (R 4.5.2)
 httpuv         1.6.17     2026-03-18 [2] CRAN (R 4.5.2)
 iterators    * 1.0.14     2022-02-05 [2] CRAN (R 4.5.1)
 kinship2     * 1.9.6.2    2025-09-04 [2] CRAN (R 4.5.1)
 knitr          1.51       2025-12-20 [2] CRAN (R 4.5.2)
 later          1.4.8      2026-03-05 [2] CRAN (R 4.5.2)
 lattice        0.22-9     2026-02-09 [2] CRAN (R 4.5.2)
 lifecycle      1.0.5      2026-01-08 [2] CRAN (R 4.5.2)
 lubridate    * 1.9.5      2026-02-04 [2] CRAN (R 4.5.2)
 magrittr       2.0.5      2026-04-04 [2] CRAN (R 4.5.3)
 MASS           7.3-65     2025-02-28 [2] CRAN (R 4.5.1)
 Matrix       * 1.7-5      2026-03-21 [2] CRAN (R 4.5.3)
 memoise        2.0.1      2021-11-26 [2] CRAN (R 4.5.1)
 mgcv         * 1.9-4      2025-11-07 [2] CRAN (R 4.5.2)
 mime           0.13       2025-03-17 [2] CRAN (R 4.5.1)
 miniUI         0.1.2      2025-04-17 [2] CRAN (R 4.5.1)
 minqa          1.2.8      2024-08-17 [2] CRAN (R 4.5.1)
 nadiv        * 2.18.0     2024-05-23 [2] CRAN (R 4.5.1)
 nlme         * 3.1-169    2026-03-27 [2] CRAN (R 4.5.3)
 numDeriv       2016.8-1.1 2019-06-06 [2] CRAN (R 4.5.1)
 otel           0.2.0      2025-08-29 [2] CRAN (R 4.5.1)
 patchwork    * 1.3.2      2025-08-25 [2] CRAN (R 4.5.1)
 pbapply        1.7-4      2025-07-20 [2] CRAN (R 4.5.1)
 pillar         1.11.1     2025-09-17 [2] CRAN (R 4.5.1)
 pkgbuild       1.4.8      2025-05-26 [2] CRAN (R 4.5.1)
 pkgconfig      2.0.3      2019-09-22 [2] CRAN (R 4.5.1)
 pkgload        1.5.2      2026-04-22 [2] CRAN (R 4.5.3)
 plyr           1.8.9      2023-10-02 [2] CRAN (R 4.5.1)
 promises       1.5.0      2025-11-01 [2] CRAN (R 4.5.2)
 proxy          0.4-29     2025-12-29 [2] CRAN (R 4.5.2)
 purrr        * 1.2.2      2026-04-10 [2] CRAN (R 4.5.3)
 quadprog     * 1.5-8      2019-11-20 [2] CRAN (R 4.5.1)
 R6             2.6.1      2025-02-15 [2] CRAN (R 4.5.1)
 rbibutils      2.4.1      2026-01-21 [2] CRAN (R 4.5.2)
 RColorBrewer   1.1-3      2022-04-03 [2] CRAN (R 4.5.1)
 Rcpp           1.1.1-1    2026-04-23 [1] Custom
 Rdpack         2.6.6      2026-02-08 [2] CRAN (R 4.5.2)
 readr        * 2.2.0      2026-02-19 [2] CRAN (R 4.5.2)
 reformulas     0.4.4      2026-02-02 [2] CRAN (R 4.5.2)
 registry       0.5-1      2019-03-05 [2] CRAN (R 4.5.1)
 remotes        2.5.0      2024-03-17 [2] CRAN (R 4.5.1)
 rlang          1.2.0      2026-04-06 [2] CRAN (R 4.5.3)
 ROI            1.0-2      2026-01-12 [2] CRAN (R 4.5.2)
 rstudioapi     0.18.0     2026-01-16 [2] CRAN (R 4.5.2)
 S7             0.2.2      2026-04-22 [2] CRAN (R 4.5.3)
 scales         1.4.0      2025-04-24 [2] CRAN (R 4.5.1)
 sequoia      * 3.2.0      2026-01-09 [2] CRAN (R 4.5.2)
 sessioninfo    1.2.3      2025-02-05 [2] CRAN (R 4.5.1)
 shiny          1.13.0     2026-02-20 [2] CRAN (R 4.5.2)
 slam           0.1-55     2024-11-13 [2] CRAN (R 4.5.1)
 snow         * 0.4-4      2021-10-27 [2] CRAN (R 4.5.1)
 spaMM        * 4.6.65     2026-04-06 [2] CRAN (R 4.5.3)
 stringi        1.8.7      2025-03-27 [2] CRAN (R 4.5.1)
 stringr      * 1.6.0      2025-11-04 [2] CRAN (R 4.5.2)
 tibble       * 3.3.1      2026-01-11 [2] CRAN (R 4.5.2)
 tidyr        * 1.3.2      2025-12-19 [2] CRAN (R 4.5.2)
 tidyselect     1.2.1      2024-03-11 [2] CRAN (R 4.5.1)
 tidyverse    * 2.0.0      2023-02-22 [2] CRAN (R 4.5.1)
 timechange     0.4.0      2026-01-29 [2] CRAN (R 4.5.2)
 tzdb           0.5.0      2025-03-15 [2] CRAN (R 4.5.1)
 usethis        3.2.1      2025-09-06 [2] CRAN (R 4.5.1)
 vctrs          0.7.3      2026-04-11 [2] CRAN (R 4.5.3)
 withr          3.0.2      2024-10-28 [2] CRAN (R 4.5.1)
 xfun           0.57       2026-03-20 [2] CRAN (R 4.5.3)
 xtable         1.8-8      2026-02-22 [2] CRAN (R 4.5.2)

 [1] /home/courtiol/R/x86_64-redhat-linux-gnu-library/4.5
 [2] /usr/local/lib/R/library
 [3] /usr/lib64/R/library
 [4] /usr/share/R/library
 * ── Packages attached to the search path.
```