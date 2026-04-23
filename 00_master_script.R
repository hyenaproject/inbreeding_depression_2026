# Main script

# This script will guide you through all the analyses.
# For it to work, please make sure you copy the entire repository and
# have the root of the repository set as working directory.

cleanup <- TRUE  ## switch to TRUE to delete files produced by this script

if (cleanup) {
  unlink("data/complete_table_Fhat3.csv")
  unlink("data/ID_Fhat3.csv")
  unlink("data/fitted_models")
  dir.create("data/fitted_models")
  unlink("figures_and_tables")
  dir.create("figures_and_tables")
}

# Download the VCF file
# if the download fails for some reasons, you can also download the file manually
# at https://zenodo.org/records/19709944
# If doing so, you must place the downloaded file within the folder data and 
# rename it as hyenas.vcf

opt <- options(timeout = 1200) # increase download timeout to 20min for slow connections
download.file("https://zenodo.org/records/19709944/files/hyenas_deIDd_Apr26.vcf?download=1",
              destfile = "data/hyenas.vcf")
options(opt) # restore default setting




# Load packages and local functions --------------------------------------------

source("01_run_always.R")


# SNP data extraction ----------------------------------------------------------

## Run only if you want to recreate the file complete_table_Fhat3.csv: (15min/4 cores)
source("02_extract_SNP_data.R")   

complete_table_Fhat3


# Stats on sample sizes mentioned in the Result section ------------------------

source("03_samplesize_stats.R")


# Plots part 1 -----------------------------------------------------------------

source("04_plotting_fig_1_S1-2.R")
Figure1
FigureS1
FigureS2


# Modelling and tables ---------------------------------------------------------
nb_cores <- 15      # adjust depending on system (we used 50)
nb_boot  <- 1000    # use 1000 for final use!

source("05_modelling.R")

Table1  <- read_table("figures_and_tables/Table1.txt") 
Table2  <- read_table("figures_and_tables/Table2.txt")
Table3  <- read_table("figures_and_tables/Table3.txt")
TableS1 <- read_table("figures_and_tables/TableS1.txt") 
TableS2 <- read_table("figures_and_tables/TableS2.txt")
TableS3 <- read_table("figures_and_tables/TableS3.txt")
TableS4 <- read_table("figures_and_tables/TableS4.txt")
TableS5 <- read_table("figures_and_tables/TableS5.txt")


# Plots part 2 -----------------------------------------------------------------

source("06_plotting_fig_2_S3-4_S7-11.R")

Figure2
FigureS3
FigureS4
FigureS7
FigureS8
FigureS9
FigureS10
FigureS11


# Plots part 3 -----------------------------------------------------------------

source("07_plotting_fig_3.R")

Figure3


# Plots part 4 -----------------------------------------------------------------

source("08_plotting_fig_S5-6.R")

FigureS5
FigureS6


# Additional results -----------------------------------------------------------

## Number of highly inbred individuals ----
complete_table_Fhat3 <- read_csv("data/complete_table_Fhat3.csv")
sum(complete_table_Fhat3$Fhat3 > 0.1, na.rm = TRUE) # 17
sum(complete_table_Fhat3$Fhat3 > 0.05, na.rm = TRUE) # 48

## Correlation Fgrm - Fped ----
signif(with(complete_table_Fhat3,
            cor(Fped, Fhat3, use = "pairwise.complete.obs")), digits = 3) # 0.588
signif(with(complete_table_Fhat3 |>
              filter(keep_sextets_and_septets),
            cor(Fped, Fhat3, use = "pairwise.complete.obs")), digits = 3) # 0.625

## Sample sizes for data with both Fped and Fgrm ----
complete_table_Fhat3 |>
  filter(!is.na(Fped), !is.na(Fhat3)) |>
  nrow() # 1119
complete_table_Fhat3 |>
  filter(keep_sextets_and_septets) |>
  filter(!is.na(Fped), !is.na(Fhat3)) |>
  nrow() # 704

## Variances Fgrm - Fped ----
signif(var(complete_table_Fhat3$Fhat3, na.rm = TRUE), digits = 3) # 0.00105
signif(var(complete_table_Fhat3$Fped, na.rm = TRUE), digits = 3) # 0.000399
signif(var(complete_table_Fhat3 |>
             filter(keep_sextets_and_septets) |>
             pull(Fped), na.rm = TRUE), digits = 3) # 0.000579


# Additional information -------------------------------------------------------

# We do not provide the full hyena database but we provide the script 
# XX_hyenaR_derived_data.R to document how we used hyenaR to create the 
# prepared dataset as well as to compute some additional statistics provided
# in method.
