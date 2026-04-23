## Note: this script cannot be run without the full database (not provided)

message("Building complete_table...")
# Extraction/filtering of population data. 
# Please call '01_run_always.R' first, in case you want to inspect this script separately:
# source("01_run_always.R")

# PLEASE NOTE: Using another 'hyenaR' version than the one provided (data/hyenaR-master_20250308.zip) may not produce the same results as this package evolves frequently. 
# Similarly, a fixed version of the database (data/database_2025_09_03.sqlite) must be used.

# Loading hyenaR --------------------------------

if (requireNamespace("hyenaR", quietly = TRUE) && 
    packageVersion("hyenaR") == "0.10.0.9000") {
  suppressMessages(library("hyenaR")) } else {
    stop(paste("Please rerun the script after installing the R package hyenaR:\n",
               "`remotes::install_local('data/hyenaR-master_20250308.zip')`"))
  }

message("Setting up CPUs for parallel computation...")
Ncpu <- max(c(parallel::detectCores() - 1, 1))
options(hyenaR_CPUcores = Ncpu)

# Loading additional functions not yet integrated in hyenaR ------------------

fetch_id_inbreeding <- function(ID) {
  
  inputID <- hyenaR::check_function_arg.ID(ID, argument.name = "ID", strict = TRUE, .fill = FALSE, arg.max.length = Inf)
  
  input.tbl <- data.frame(ID = inputID)
  
  if (!requireNamespace("ggroups", quietly = TRUE)) {
    stop("Please install R package 'ggroups' to use this function")
  }
  ## build a pedigree to start from
  pedigree <- hyenaR::create_pop_pedigree(pedigree.format = "both")
  
  ## select columns needed for the package ggroups
  pedigree[, c("ID_num", "sire_num", "dam_num")] -> pedigree_ggroups
  
  ## compute inbreeding
  ggroups::inbreed(ped = as.data.frame(pedigree_ggroups)) -> inbr_coeff
  
  ## format to typical hyenaR tibble::tibble
  tibble::tibble(ID = pedigree$ID, inbr_coeff) -> inbr_tbl
  
  hyenaR::check_function_output(input.tbl = input.tbl, output.tbl = inbr_tbl, join.by = "ID",
                                duplicates = "input", output.IDcolumn = "inbr_coeff")
}

find_id_id.ancestor2 <- function(ID, max.depth = Inf) {
  
  if (max.depth <= 0) {
    character(0)
  } else {
    mum <- hyenaR::fetch_id_id.mother.genetic(ID)
    dad <- hyenaR::fetch_id_id.father(ID)
    
    ancestors <- c(mum = mum, dad = dad)
    ancestors_todo <- ancestors[!is.na(ancestors)]
    
    if (length(ancestors_todo) == 0) {
      character(0)
    } else {
      ancestor_list <- lapply(ancestors_todo, function(p) find_id_id.ancestor2(ID = p, max.depth = max.depth - 1))
      ancestors <- c(ancestors, unlist(ancestor_list))
    }
    ancestors
  }
}


# Extraction of life history details for all IDs --------------------------

## Life history data for all IDs ----
starting_table <- create_id_starting.table() # all IDs ever observed

starting_table |>
  mutate(mother = fetch_id_id.mother.genetic(ID),
         mother_social = fetch_id_id.mother.social(ID),
         father = fetch_id_id.father(ID),
         sex = fetch_id_sex(ID),
         birth_clan = fetch_id_clan.birth(ID),
         birth_date = fetch_id_date.birth(ID),
         birth_year = year(birth_date),
         mother_rank_ord = fetch_id_rank.native(ID = mother, at = birth_date), # nb of native in clan is bottom, 1 is top
         mother_rank_std = fetch_id_rank.native.std(ID = mother, at = birth_date), # -1 is bottom, 1 is top
         mother_rank_std_0_1 = (mother_rank_std + 1)/2, # 0 is bottom, 1 is top
         lifespan = fetch_id_duration.lifespan(ID,unit = "year"),
         left_censored = fetch_id_is.censored.left(ID),
         right_censored = fetch_id_is.censored.right(ID),
         Fped = fetch_id_inbreeding(ID),
         total_number_cubs = fetch_id_number.offspring(ID)) |> 
 rowwise() |> 
 mutate(ancestors = list(find_id_id.ancestor2(ID = ID, max.depth = 2))) |> 
 ungroup() |>
 mutate(n_ancestors = sapply(ancestors, function(x) sum(!is.na(x)))) -> complete_table_1

## Remove individuals that have nothing to do with the Crater clans ----
ID_to_keep <- c(create_id_starting.table(main.clans.birth = TRUE)$ID,
                create_id_starting.table(lifestage = "foreigner")$ID,
                create_id_starting.table(lifestage = "founder_male")$ID)

ID_to_keep <- unique(ID_to_keep)
length(ID_to_keep) # 3137

complete_table_1 |> 
  filter(ID %in% ID_to_keep) -> complete_table_2

## Add clan size ----
# (slow, 7min/15 cores --- parallel processing possible)
# note: count same IDs as included in native ranks
message("Adding clan size... (be patient, it may take several minutes)")
fetch_clan_number(clan = "A", at = "2000-01-01") ## Blank run for caching purposes
system.time(
  complete_table_2 |>
    mutate(clansize_subadult    = fetch_clan_number(clan = birth_clan, at = birth_date, lifestage = "subadult"),
           clansize_natal       = fetch_clan_number(clan = birth_clan, at = birth_date, lifestage = "natal"),
           clansize_philopatric = fetch_clan_number(clan = birth_clan, at = birth_date, lifestage = "philopatric"),
           clansize_native      = clansize_subadult + clansize_natal + clansize_philopatric,
           clansize_all         = fetch_clan_number(clan = birth_clan, at = birth_date)) |>
    select(-clansize_natal, -clansize_philopatric, -clansize_subadult) -> complete_table_3
  )


# Selection criteria ------------------------------------------------------

## Here we add 'keep columns' to make the filtering processes clear.
## We will filter the data right before modelling/plotting, etc.
## note: 'TRUE' in these columns means: to keep!


## Keep column explanation ----

# We have 'keep columns' for the following filtering: 
#  1. Born inside the crater
#  2. Not adopted
#  3. All grandparents/parents known
#  4. One unknown grandparent, all others known
#  5. Could have been followed for at least 4 years
#  6. Could have been followed for at least 8 years
#  7. Could have been followed for at least 10 years (our choice for main text, see extract_pop_data_companion_Rscripts/check_censorship.R)
#  8. Could have been followed for at least 12 years
#  9. Did survive to 24 months
# 10. Is not left censored
# 11. Is not right censored
# 12. !Additional Fhat3 column and keep column for all genotyped IDs are created and added in '03_extract_SNP_data.R'!

# Note about Crater born definition: by design in the data, left-censored individuals
# are considered to be born in their clan _only_ if they were females
# are considered to be born in their clan _only_ if they were females
# or males under 4 years at the date of first observation.

complete_table_3 |> 
  mutate(## flag individuals born inside the main clans with TRUE 
         keep_mainclan_born = birth_clan %in% find_clan_name.all(main.clans = TRUE),
         ## flag non-adopted IDs with TRUE 
         keep_non_adopted = mother == mother_social, 
         ## flag IDs which could have been followed for at least 4 years with TRUE (born on or before "2023-09-01") 
         keep_followed_atleast4yrs = birth_date <= (find_pop_date.observation.last() - years(4)),
         ## flag IDs which could have been followed for at least 8 years with TRUE (born on or before "2017-09-01") 
         keep_followed_atleast8yrs = birth_date <= (find_pop_date.observation.last() - years(8)), # 84.31% of pop. has full lifespan
         ## flag IDs which could have been followed for at least 10 years with TRUE (born on or before "2015-09-01") 
         keep_followed_atleast10yrs = birth_date <= (find_pop_date.observation.last() - years(10)), # 89.50% of pop. has full lifespan
         ## flag IDs which could have been followed for at least 12 years with TRUE (born on or before "2013-09-01") 
         keep_followed_atleast12yrs = birth_date <= (find_pop_date.observation.last() - years(12)), # 93.88% of pop. has full lifespan
         ## flag IDs which did survive to 24 months with TRUE
         keep_survived_to_24m = fetch_id_is.alive(ID, at = birth_date + months(24)),
         ## flag IDs that have been born after the first day of observation (born before "1996-04-12") with TRUE
         keep_left_uncensored = !left_censored,
         ## flag IDs that have been alive on the last day of observation with TRUE
         keep_right_uncensored = !right_censored,
         ## flag IDs with max: one unknown grandparent with TRUE 
         keep_sextets_and_septets = n_ancestors > 4) -> complete_table

## Quick checks of rank / clan size consistency ----
# note: there should be no point above the diagonal
complete_table |> 
  filter(keep_mainclan_born) |> 
  plot(mother_rank_ord ~ clansize_native,
       col = ifelse(mother_rank_ord > clansize_native, 2, 1),
       cex = 0.5,
       data  = _)
  abline(0, 1, col = "blue")
# -> few very small error at small clan sizes caused by inconsistent treatment of selector_X within hyenaR (#751) -> acceptable

# Store file 
write_csv(complete_table, "data/complete_table.csv") ## Note: this destroys the ancestor column which is a list column


message("Done!")
  

# Additional demographic statistics provided in Methods ------------------------
  
first_obs <- "1996-04-12" # find_pop_date.observation.first()
last_obs <- "2025-09-01" # find_pop_date.observation.last()

# Paragraph Study population and pedigree ----

## Numbers of hyenas in pop/clans ####
start_tbl_date <- reshape_row_date.seq(tbl = tibble::tibble(date = NA),
                                       from = first_obs,
                                       to = last_obs,
                                       by = "1 day")

complete_table |> 
  select(ID) -> start_tbl_ID

start_tbl_combined <- tidyr::crossing(start_tbl_ID, start_tbl_date)

message("Computing variation in pop size over the study period...")
start_tbl_combined |> 
  mutate(ID_alive = fetch_id_is.alive(ID, at = date)) |> 
  filter(ID_alive == TRUE) |> 
  mutate(clan_current = fetch_id_clan.current(ID, at = date)) |> 
  filter(clan_current %in% find_clan_name.all(main.clans = TRUE)) |> 
  summarize(n_alive = n(), .by = "date") |> 
  summarize(min_alive = min(n_alive),
            max_alive = max(n_alive)) -> Ids_in_pop

message("\nRange of pop size:")
print(Ids_in_pop) # 165 -- 609

message("Mean clan sizes over the study period...")
start_tbl_combined |> 
  filter(fetch_id_is.alive(ID, at = date)) |> 
  mutate(clan_current = fetch_id_clan.current(ID, at = date)) |> 
  filter(clan_current %in% find_clan_name.all(main.clans = TRUE)) |> 
  summarize(n_alive = n(), .by = c("date", "clan_current")) |> 
  summarize(mean_alive = mean(n_alive),
            .by = clan_current) |> 
  arrange(mean_alive) -> IDs_per_clan # runs 1 min

message("\nMean clan size:")
print(IDs_per_clan) # 26 -- 73 (detailed by clan)


## Number of hyenas sampled for DNA/microsats ####
N_DNA <- sum(fetch_id_is.sampled.dna(start_tbl_ID$ID)) # 2181
message("\nNumber of hyenas sampled for DNA:")
print(N_DNA)


## Litter size ####
create_litter_offspring.count() |> 
  select(female, male, unknown) |> 
  rowwise() |> 
  mutate(Noffspring = as.character(sum(c_across(everything())))) |> 
  ungroup() |> 
  filter_out(Noffspring == "0") |> 
  count(Noffspring) |> 
  mutate(prop = n/sum(n)) -> litter_size

litter_size |> 
  filter(Noffspring < 3) |>
  summarise(n = sum(n), prop = sum(prop)) |> 
  mutate(Noffspring = "1 or 2 (pooled)")  -> litter_size_12

litter_size |> 
  summarise(n = sum(n), prop = sum(prop)) |> 
  mutate(Noffspring = "all (pooled)")  -> litter_size_all

litter_size |> 
  bind_rows(litter_size_12) |> 
  bind_rows(litter_size_all) -> litter_size_final

message("\nNumber of litters per genetic offspring number:")
print(litter_size_final)
# A tibble: 5 × 3
#   Noffspring          n    prop
#   <chr>           <int>   <dbl>
# 1 1                1003 0.528  
# 2 2                 888 0.467  
# 3 3                  10 0.00526
# 4 1 or 2 (pooled)  1891 0.995  
# 5 all (pooled)     1901 1 


## Litter paternity ####
create_litter_offspring.table() |> 
  filter(filiation != "mother_social") |> 
  mutate(dadID = fetch_id_id.father(offspringID)) |> 
  filter_out(is.na(dadID)) |> 
  summarise(Ndad = length(unique(dadID)),
            Noffspring = length(unique(offspringID)),
            .by = litterID) |> 
  filter_out(Noffspring == 1) |> 
  summarise(Nsingledad = sum(Ndad == 1),
            prop_singledad = mean(Ndad == 1),
            Nlitter = n()) -> single_sire

message("\nInfo on single-sired genetic litters")
print(single_sire)
# A tibble: 1 × 3
#     Nsingledad prop_singledad Nlitter
#          <int>          <dbl>   <int>
#   1        415          0.849     489


# Paragraph Study population and pedigree ----

## Total number of hyenas in the database ####
N_total <- nrow(complete_table)
message("\nTotal number of hyenas in the database:")
print(N_total) # 3137

message("\nTotal number of hyenas born in main clan:")
print(sum(complete_table$keep_mainclan_born)) # 3012

complete_table |> 
  mutate(lifestage = fetch_id_lifestage(ID, at = first_obs)) |> 
  filter(lifestage == "founder_male") |> 
  nrow() -> N_founder_males
message("\nNumber of founder males")
print(N_founder_males) # 32

complete_table |> 
  mutate(founder_male = fetch_id_lifestage(ID, at = first_obs) == "founder_male") |> 
  filter(!founder_male | is.na(founder_male), !keep_mainclan_born) |> 
  nrow() -> N_outsider_integrated
message("\nNumber of outsider integrating pedigree")
print(N_outsider_integrated) # 93


###### Maximum generations of hyenas #### 
pedigree <- create_pop_pedigree()
pedigree |> 
  filter(ID %in% complete_table$ID) -> pedigree

### 1) Pedigree depth over the entire population ####
generations_all <- kinship2::kindepth(id = pedigree$ID,
                                      dad.id = pedigree$sire,
                                      mom.id = pedigree$dam) 

message("\nMax number of generations in pedigree: ")
print(max(generations_all)) # 9 generations not counting the founding cohort

### 2) Pedigree depth for alive, main clan, individuals (sensu right-censored) ####
generations_alive <- data.frame(ID = pedigree$ID,
                                generation = kinship2::kindepth(id = pedigree$ID,
                                                                dad.id = pedigree$sire,
                                                                mom.id = pedigree$dam)) |> 
  filter(fetch_id_is.alive(ID, at = last_obs) & fetch_id_clan.birth(ID) %in% find_clan_name.all(main.clans = TRUE))

message("\nMean pedigree depth for last generation: ")
print(signif(mean(generations_alive$generation), digits = 3)) # 5.36


## Completeness of the pedigree ####
pedigree |> 
  count(non_missing_mum = !is.na(dam), non_missing_dad = !is.na(sire)) |> 
  mutate(prop = n/sum(n)) -> ped_completeness

message("\nPedigree completeness:")
print(ped_completeness)
# # A tibble: 4 × 4
#   non_missing_mum non_missing_dad     n    prop
#   <lgl>           <lgl>           <int>   <dbl>
# 1 FALSE           FALSE             326 0.104  
# 2 FALSE           TRUE                4 0.00128
# 3 TRUE            FALSE             964 0.307  
# 4 TRUE            TRUE             1843 0.588 


pedigree |> 
  filter(is.na(sire), is.na(dam)) |> 
  mutate(founder = fetch_id_is.founder(ID)) |> 
  count(founder) -> founders_no_parents

message("\nNumber of founders with no known parents:")
print(founders_no_parents)
# # A tibble: 2 × 2
#   founder     n
#   <lgl>   <int>
# 1 FALSE     192
# 2 TRUE      134

## Average age at first sighting ####
create_id_starting.table(clan.birth = find_clan_name.all(main.clans = TRUE), 
                         lifestage = "cub", 
                         lifestage.overlap = "start", 
                         from = first_obs, 
                         to = last_obs) |> 
  mutate(first_sighting = fetch_id_date.observation.first(ID),
         birth_date = fetch_id_date.birth(ID)) |> 
  mutate(age_at_first_sighting = fetch_id_age(ID, at = first_sighting, unit = "months")) |> 
  summarize(avg_age = mean(age_at_first_sighting),
            sd_age = sd(age_at_first_sighting),
            median_age = median(age_at_first_sighting)) -> avg_age_first_obs

message("Average age at first sighting (in months!):")
print(avg_age_first_obs)
# # A tibble: 1 × 3
#     avg_age sd_age median_age
#       <dbl>  <dbl>      <dbl>
#   1    3.11   3.50       1.97

message("Done!")
  
