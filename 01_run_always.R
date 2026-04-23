# Libraries ####
std_pkgs_needed <- c("tidyverse", "kinship2", "purrr", "spaMM", "mgcv", "doSNOW", "patchwork",
                     "sequoia", "ggroups", "ggExtra", "nadiv") 

# install.packages(std_pkgs_needed) ## shortcut to install all standard packages at once, if needed. Uncomment and run this line if you need to install all the packages at once.

sink <- lapply(std_pkgs_needed, \(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressMessages(library(pkg, character.only = TRUE, quietly = TRUE))
    } else {
      stop(paste0("Please rerun the script after installing the R package ", pkg,
                  ":\n `install.packages('", pkg, "')`"))
    }
  })
rm(sink) # delete variable storing loaded pkg names

# sourcing functions
message("Sourcing helper functions...")
source("local_functions.R")

# checking that the VCF file is at the right place
message("Checking presence of VCF file...")
if (!file.exists("data/hyenas.vcf")) warning("No file called hyenas.vcf was found in /data")

message("Done!")
