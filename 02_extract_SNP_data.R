# Extraction/filtering of SNP data. 
# Please call '01_run_always.R' first, in case you want to inspect this script separately:
source("01_run_always.R") 

# To run the below code, you need: 
# 1. The file hyenas.vcf placed in the folder data (see below for download)
# 2. 'plink' to be installed on your machine (https://www.cog-genomics.org/plink/).
# (for Linux, after extracting the program, run `sudo cp plink /usr/bin/` in a terminal so that calls below work)
# (for Windows, store plink.exe in the working directory of this project)

if (system("plink --version") != 0) {
  stop("plink is not accessible. Please check installation notes.")
}

# Download the VCF file --------------------------------------------------------
# if the download fails for some reasons, you can also download the file manually
# at https://zenodo.org/records/19709944
# If doing so, you must place the downloaded file within the folder data and 
# rename it as hyenas.vcf

opt <- options(timeout = 1200) # increase download timeout to 20min for slow connections
download.file("https://zenodo.org/records/19709944/files/hyenas_deIDd_Apr26.vcf?download=1",
              destfile = "data/hyenas.vcf")
md5 <- tools::md5sum("data/hyenas.vcf")
if (md5 != "3def661153f84cce687a025dec9a262f") {
  file.remove("data/hyenas.vcf")
  message("MD5 mismatch! Please download the file manually, rename it to hyenas.vcf and store it in the folder `data`.")
} else {
  message("MD5 check passed! The file is correct.")
}
options(opt) # restore default setting

# for maintainers only:
#path_vcf <- "private/Hyenas_2ndRun_mincov10MaxCov110.vcf"
#complete_table <- read_csv("private/complete_table_original.csv")

complete_table <- read_csv("data/complete_table.csv")
path_vcf <- "data/hyenas.vcf"

message("Adding Fgrm to complete_table and investigating SNP/IDs selection... (check script for details)")

# Data preparation --------------------------------------------------------
# Create a scratch directory for temporary output 
dir.create("./scratch")

# Extract IDs born in the crater
craterIDs <- complete_table |> 
  filter(keep_mainclan_born) |> 
  select(ID)  

# Store the craterIDs temporarily as txt-file (for plink)
craterIDs_df <- data.frame(FID = fetch_FID_plink(craterIDs$ID), IID = fetch_IID_plink(craterIDs$ID))
write.table(craterIDs_df, "scratch/craterIDs.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Build pedigree (for sequoia)
base_pedigree <- complete_table |> 
  arrange(birth_date) |> 
  select(ID, father, mother)


# Fgrm/Fhat3 value extraction ---------------------------------------------

## 1. Individual genotype missingness (MIND) - MAX 0.5 ---- 

system(paste("plink --vcf", path_vcf, "--allow-extra-chr --missing --out scratch/missing")) 
# 69816 variants and 1181 IDs loaded

base_genomic_mind <- read.table("scratch/missing.imiss", header = TRUE)
hist(base_genomic_mind$F_MISS, xlab = "Frequency of variants not genotyped")
all_genotyped_IDs <- get_ID_plink(base_genomic_mind)

### Removing samples with missing call frequencies greater than a threshold of 0.5 ----
system(paste("plink --vcf", path_vcf, " --allow-extra-chr --mind 0.5 --make-bed --out scratch/MIND_filtered")) 
# 69816 variants, 1167 IDs remaining

#### Checking basic characteristics of 14 dropped individuals
read.table("scratch/MIND_filtered.fam") |> 
  rename(FID = V1, IID = V2) |> 
  get_ID_plink() -> IDs_kept_good_genotyped
IDs_dropped_missing_genotyped <- setdiff(all_genotyped_IDs, IDs_kept_good_genotyped)
complete_table |> 
  filter(ID %in% IDs_dropped_missing_genotyped) |> 
  select(ID, sex, birth_clan, birth_year, Fped, mother_rank_std_0_1) |> 
  print(n = Inf)

### Removing individuals born outside of the crater ----
system("plink --bfile scratch/MIND_filtered --allow-extra-chr --keep scratch/craterIDs.txt --make-bed --out scratch/MIND_crater_filtered")
# 69816 variants, 1119 IDs remaining

#### Checking basic characteristics of 41 dropped individuals
read.table("scratch/MIND_crater_filtered.fam") |> 
  rename(FID = V1, IID = V2) |> 
  get_ID_plink() -> IDs_kept_within_crater
IDs_dropped_outside_crater <- setdiff(IDs_kept_good_genotyped, IDs_kept_within_crater)
complete_table |> 
  filter(ID %in% IDs_dropped_outside_crater) |> 
  select(ID, sex, birth_clan, birth_year, Fped, mother_rank_std_0_1) |> 
  print(n = Inf)


## 2. Mendelian error (ME) - filter SNPs with 1% ME rate or more -----------

system(paste("plink --vcf", path_vcf, "--allow-extra-chr --recode A --out scratch/snp")) 
# 69816 variants, 1181 IDs loaded

### Convert genotype data in various formats (Slowish: 1-3min) ----
system.time(temp_GenoM <- GenoConvert(InFile = "scratch/snp.raw", InFormat = "raw"))
# 1181 IDs, 60247  SNPs

### Fix IDs ----
rownames(temp_GenoM) <- complete_IID_CroCro(rownames(temp_GenoM))

### Estimate allele frequency (AF), missingness and Mendelian errors per SNP ----
# (Slow: 2-8min)
base_pedigree <- as.data.frame(base_pedigree)
system.time(base_genomic_ME <- SnpStats(Pedigree = base_pedigree, # Entire pedigree of all hyenas. Non_genotyped IDs are ignored 
                                        GenoM = temp_GenoM, 
                                        Plot = FALSE,
                                        ErrFlavour = NULL)) # ErrFlavour = DEPRECATED AND IGNORED but needed.

### Mendelian error rate (as a proportion) per parent–offspring pair ----
base_genomic_ME$MEpercent <- base_genomic_ME$MEpair/(base_genomic_ME$n.dam + base_genomic_ME$n.sire)

### Select ME rates greater than 1% ----
base_genomic_ME_0.01 <- row.names(base_genomic_ME[base_genomic_ME$MEpercent > 0.01, ])

### Store the ME rates greater than 1% temporarily as txt-file (for plink)
write.table(base_genomic_ME_0.01, 
            file = "scratch/ME_0.01.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

### Filtering out all variants with ME rates greater than 1% ----
system("plink --bfile scratch/MIND_crater_filtered --allow-extra-chr --exclude scratch/ME_0.01.txt --make-bed --out scratch/ME_filtered_0.01") 
# 66577 variants and 1119 IDs remains


## 3. Local genotype missingness (GENO) & MAF filtering --------------------

### Missingness @ 0.2 and MAF @ 1% -  ----
system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.2 --maf 0.01 --ibc --out scratch/maf_0.01_miss_0.2") 
# 21929 variants and 1119 IDs

### Test effect of missingness filter - try 0.2, 0.4, 0.6, 0.8, 1 ----
## Calculate Fgrm using different missingness thresholds - ensuring that maf is always kept at 1% 

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.05 --maf 0.01 --ibc --out scratch/maf_0.01_miss_0.05") 
# 3513 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.1 --maf 0.01 --ibc --out scratch/maf_0.01_miss_0.1")   
# 11902 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.4 --maf 0.01 --ibc --out scratch/maf_0.01_miss_0.4")   
# 30925 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.6 --maf 0.01 --ibc --out scratch/maf_0.01_miss_0.6")   
# 36243 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.8 --maf 0.01 --ibc --out scratch/maf_0.01_miss_0.8")   
# 37032 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 1 --maf 0.01 --ibc --out scratch/maf_0.01_miss_1")       
# 37036 variants and 1119 IDs

### Read in all data frames to compare Fgrm estimates with different missingness levels

 # list of files
miss_files <- list.files(path = "scratch", pattern = "miss_.*\\.ibc$", full.names = TRUE)
 
 # sth to start with:
Fgrm_miss <- NULL
 
for (file_path in miss_files) {
  # extract missingness value from filename
  miss_val <- str_extract(file_path, "(?<=miss_)[0-9.]+") |> 
    str_replace("\\.$", "") # get rid of the dot
  
  # read files and select FID and Fhat3 columns
  miss_tbl <- read_table(file_path, show_col_types = FALSE) |>
    mutate(ID = get_ID_plink(.data)) |> 
    select(ID, !!paste0("miss_", miss_val) := Fhat3)  # rename Fhat3 (!! to unquote, := to assign name dynamically)
  
  # join the tables by FID
  if (is.null(Fgrm_miss)) {
    Fgrm_miss <- miss_tbl
  } else {
    Fgrm_miss <- left_join(Fgrm_miss, miss_tbl, by = "ID")
  }
}
 
Fgrm_miss
 
### Correlation table MISS ----
cor_miss <- data.frame(cor(Fgrm_miss[, -1L]))
cor_miss <- round(cor_miss,3)
cor_miss

# Save the table for the manuscript
write.table(cor_miss, file = "figures_and_tables/TableS1.txt", sep = "\t", quote = FALSE)


## 4. Minor allele frequencies (MAF) ---------------------------------------

### Test effect of maf filter -  0.05, 0.1, 0.2, 0.3, 0.4 ----
## Calculate Fgrm using different maf thresholds - ensuring that missingness is always kept at 20% 

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.2 --maf 0.05 --ibc --out scratch/maf_0.05_miss_0.2")  
# 17234 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.2 --maf 0.1 --ibc --out scratch/maf_0.1_miss_0.2")    
# 13622 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.2 --maf 0.2 --ibc --out scratch/maf_0.2_miss_0.2")    
# 8717 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.2 --maf 0.3 --ibc --out scratch/maf_0.3_miss_0.2")    
# 5285 variants and 1119 IDs

system("plink --bfile scratch/ME_filtered_0.01 --allow-extra-chr --geno 0.2 --maf 0.4 --ibc --out scratch/maf_0.4_miss_0.2")    
# 2592 variants and 1119 IDs


### Read in all data frames to compare Fgrm estimates with different maf levels

# list of files
maf_files <- list.files(path = "scratch", pattern = "maf_.*\\_0.2.ibc$", full.names = TRUE)

# sth to start with:
Fgrm_maf <- NULL

for (file_path in maf_files) {
  # extract MAF value from filename
  maf_val <- str_extract(file_path, "(?<=maf_)[0-9.]+") |> 
    str_replace("\\.$", "") # get rid of the dot
  
  # read files and select FID and Fhat3 columns
  maf_tbl <- read_table(file_path, show_col_types = FALSE) |>
    mutate(ID = get_ID_plink(.data)) |> 
    select(ID, !!paste0("maf_", maf_val) := Fhat3)  # rename Fhat3 (!! to unquote, := to assign name dynamically)
  
  # join the tables by FID
  if (is.null(Fgrm_maf)) {
    Fgrm_maf <- maf_tbl
  } else {
    Fgrm_maf <- left_join(Fgrm_maf, maf_tbl, by = "ID")
  }
}

Fgrm_maf


### Correlation table MAF ----
cor_maf <- data.frame(cor(Fgrm_maf[, -1L]))
cor_maf <- round(cor_maf, 3)
cor_maf

# Save the table for the manuscript
write.table(cor_maf, file = "figures_and_tables/TableS2.txt", sep = "\t", quote = FALSE)

### Decision for Fgrms estimates generated using 21929 variants and 1119 IDs (miss 0.2, maf 0.01) ----
# Read in the Fgrms we chose:
# * MAF:  trying to keep rare alleles, so low maf selected since correlated with higher threshold anyhow
# * MISS: trying to keep as many SNPs as possible while maximising coverage
read_table("scratch/maf_0.01_miss_0.2.ibc", col_names = TRUE) |> 
  mutate(ID = get_ID_plink(.data)) |> 
  select(ID, Fhat3) -> ID_Fhat3

### Store this table
write_csv(ID_Fhat3, "data/ID_Fhat3.csv")


## 5. Join the Fhat3 values to our 'complete_table' ------------------------

complete_table |> 
  left_join(ID_Fhat3, by = "ID") |> 
  mutate(keep_genotyped = ID %in% all_genotyped_IDs,
         keep_not_too_negative_Fgrm = Fhat3 > -0.1) -> complete_table_Fhat3

### Checking basic characteristics of 2 individuals with extreme negative Fhat3
### (likely reflecting sequencing errors in our case)
complete_table_Fhat3 |> 
  filter(!keep_not_too_negative_Fgrm) |> 
  select(ID, Fhat3, sex, birth_clan, birth_year, Fped, mother_rank_std_0_1) |> 
  print(n = Inf)

### Store final table
write_csv(complete_table_Fhat3, "data/complete_table_Fhat3.csv") 

## Eliminate scratch folder and remove all temporary files
unlink("./scratch", recursive = TRUE)

message("Done!")
