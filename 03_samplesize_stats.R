# This is a short script to verify the demographic numbers mentioned in the manuscript.
# Please call '01_run_always.R' first, in case you want to inspect this script separately:
# source("01_run_always.R")

message("Computing demographics...")

complete_table <- read_csv("data/complete_table.csv")
first_obs <- "1996-04-12" # find_pop_date.observation.first()
last_obs <- "2025-09-01" # find_pop_date.observation.last()


# Paragraph Data selection ----

## Sample size for individuals born and dead in main clan between first and 10 years before last observation ####
complete_table |> 
  filter(keep_mainclan_born, keep_followed_atleast10yrs, keep_left_uncensored, keep_right_uncensored) |> 
  nrow() -> N_mainclan_10yrs

message("\nNumber of individuals born and dead in main clans between first and 10 years before last observation:")
print(N_mainclan_10yrs) # 1748


## Sample size for individuals born and dead in main clan between first and 12 years before last observation ####
complete_table |> 
  filter(keep_mainclan_born, keep_followed_atleast12yrs, keep_left_uncensored, keep_right_uncensored) |> 
  nrow() -> N_mainclan_12yrs

message("\nNumber of individuals born and dead in main clans between first and 12 years before last observation:")
print(N_mainclan_12yrs) # 1544


## Sample size for individuals born and dead in main clan between first and 8 years before last observation ####
complete_table |> 
  filter(keep_mainclan_born, keep_followed_atleast8yrs, keep_left_uncensored, keep_right_uncensored) |> 
  nrow() -> N_mainclan_8yrs

message("\nNumber of individuals born and dead in main clans between first and 8 years before last observation:")
print(N_mainclan_8yrs)  # 2013


## Sample size for individuals born in main clan between first and 4 years before last observation ####
complete_table |> 
  filter(keep_mainclan_born, keep_followed_atleast4yrs, keep_left_uncensored) |> 
  nrow() -> N_mainclan_4yrs_nodeathrestriction

message("\nNumber of individuals born in main clans between first and 4 years before last observation (no death restriction):")
print(N_mainclan_4yrs_nodeathrestriction) # 2526


## Effect of dropping adoptees in LRS/lifespan dataset ####
complete_table |> 
  filter(keep_mainclan_born, keep_followed_atleast10yrs, keep_left_uncensored, keep_right_uncensored, keep_non_adopted) |> 
  nrow() -> N_mainclan_10yrs_nonadoptees

message("\nNumber of individuals dropped due to adoption in LRS/lifespan dataset:")
print(N_mainclan_10yrs - N_mainclan_10yrs_nonadoptees) # 83


## Effect of dropping adoptees in juv survival dataset ####
complete_table |> 
  filter(keep_mainclan_born, keep_followed_atleast4yrs, keep_left_uncensored, keep_non_adopted) |> 
  nrow() -> N_mainclan_4yrs_nodeathrestriction_nonadoptees

message("\nNumber of individuals dropped due to adoption in juv survival dataset:")
print(N_mainclan_4yrs_nodeathrestriction - N_mainclan_4yrs_nodeathrestriction_nonadoptees) # 137


message("\nFor specific sample sizes used for modelling, see 05_modelling.R")
