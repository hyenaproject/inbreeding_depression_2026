# Fit the models, extract information from fits, and create tables deriving from those
# Note that this script also contains the supplementary models used for sensitivity analyses

# Please call '01_run_always.R' first, in case you want to run this script separately:
# source("01_run_always.R")
# nb_cores <- 50     # adjust depending on system (too much parallelism can be slower here)
# nb_boot  <- 1000   # use 1000 for final use!


if (!dir.exists("data/fitted_models")) {
  dir.create("data/fitted_models")
}


# Data preparation/filtering ---------------------------------------------------

## load full data set ----
complete_table_Fhat3 <- read_csv("data/complete_table_Fhat3.csv")

## Data for LRS and lifespan models ----

### Common datasets for adults ----
complete_table_Fhat3 |>
  filter(keep_mainclan_born & 
         keep_non_adopted & 
         keep_left_uncensored &
         keep_right_uncensored) |> 
  mutate(LRS = total_number_cubs,
         across(c(sex, birth_clan, mother, birth_year), as.factor)) |> 
  select(ID, LRS, lifespan, Fhat3, Fped,
         keep_followed_atleast10yrs, keep_followed_atleast8yrs, keep_followed_atleast12yrs,
         keep_sextets_and_septets, keep_not_too_negative_Fgrm,
         sex, mother_rank_ord, clansize_all, birth_clan, mother, birth_year) |> 
  drop_na(ID, sex, mother_rank_ord, clansize_all, birth_clan, mother, birth_year) -> y_adult

apply(y_adult, 2, \(x) sum(is.na(x))) #check which columns contain NAs (Fhat3 and keep_not_too_negative_Fgrm)

y_adult |> 
  filter(keep_followed_atleast10yrs) |> 
  select(-contains("followed")) |> 
  droplevels() -> y10_adult

y_adult |> 
  filter(keep_followed_atleast8yrs) |> 
  select(-contains("followed")) |> 
  droplevels() -> y8_adult

y_adult |> 
  filter(keep_followed_atleast12yrs) |> 
  select(-contains("followed")) |> 
  droplevels() -> y12_adult

nrow(y10_adult) # 1458
nrow(y8_adult) # 1668
nrow(y12_adult) # 1285

### Transform lifespan for better modelling ----
#apply a rank transformation and store the transformed value in the original data
#note: we're doing this on full data prior to subsetting
#so that individuals phenotype values between Fped and Fgrm are consistent

y10_adult$rt_lifespan <- qnorm(rank(y10_adult$lifespan)/max(rank(y10_adult$lifespan) + 1))
y8_adult$rt_lifespan  <- qnorm(rank(y8_adult$lifespan)/max(rank(y8_adult$lifespan) + 1))
y12_adult$rt_lifespan <- qnorm(rank(y12_adult$lifespan)/max(rank(y12_adult$lifespan) + 1))

# We use a gam to obtain un-transformed lifespan values based on any transformed lifespan
# NB: we use a high k value for a close fit

fit_gam_inverse10y <- mgcv::gam(lifespan ~ s(rt_lifespan, k = 30), data = y10_adult)
saveRDS(fit_gam_inverse10y, file = "data/fitted_models/fit_gam_inverse10y.RDS")
fit_gam_inverse8y  <- mgcv::gam(lifespan ~ s(rt_lifespan, k = 30), data = y8_adult)
saveRDS(fit_gam_inverse8y, file = "data/fitted_models/fit_gam_inverse8y.RDS")
fit_gam_inverse12y <- mgcv::gam(lifespan ~ s(rt_lifespan, k = 30), data = y12_adult)
saveRDS(fit_gam_inverse12y, file = "data/fitted_models/fit_gam_inverse12y.RDS")

## Fgrm data for adults ----
y10_adult |>
  select(-keep_sextets_and_septets) |> 
  filter(keep_not_too_negative_Fgrm) |> 
  drop_na() -> y10_adult_Fgrm

y8_adult |>
  select(-keep_sextets_and_septets) |> 
  filter(keep_not_too_negative_Fgrm) |> 
  drop_na() -> y8_adult_Fgrm

y12_adult |>
  select(-keep_sextets_and_septets) |> 
  filter(keep_not_too_negative_Fgrm) |> 
  drop_na() -> y12_adult_Fgrm

nrow(y10_adult_Fgrm) # 768
nrow(y8_adult_Fgrm)  # 864
nrow(y12_adult_Fgrm) # 715

table(y10_adult_Fgrm$sex) # 389 females + 379 males
table(y8_adult_Fgrm$sex)  # 435 females + 429 males
table(y12_adult_Fgrm$sex) # 364 females + 351 males

signif(mean(y10_adult_Fgrm$LRS), digits = 3) # 2.96
signif(mean(y8_adult_Fgrm$LRS), digits = 3)  # 2.66
signif(mean(y12_adult_Fgrm$LRS), digits = 3) # 3.13

signif(mean(y10_adult_Fgrm$lifespan), digits = 3) # 6.39
signif(mean(y8_adult_Fgrm$lifespan), digits = 3)  # 5.93
signif(mean(y12_adult_Fgrm$lifespan), digits = 3) # 6.55

## Fped data for adults ----
y10_adult |>
  select(-Fhat3, -keep_not_too_negative_Fgrm) |> 
  filter(keep_sextets_and_septets) |> 
  drop_na() -> y10_adult_Fped

nrow(y10_adult_Fped) # 793

table(y10_adult_Fped$sex) # 416 females + 377 males
signif(mean(y10_adult_Fped$LRS), digits = 3) # 2.07
signif(mean(y10_adult_Fped$lifespan), digits = 3) # 5.23

## Create sparse inverse of relatedness matrix for animal model ----
complete_table_Fhat3 |>
  arrange(birth_date) |>
  select(ID, mother, father) |> 
  as.data.frame() |> 
  nadiv::makeA() |> 
  as_precision() -> invAmat

## Common data for juveniles ----
complete_table_Fhat3 |>
  filter(keep_mainclan_born & 
         keep_followed_atleast4yrs &
         keep_non_adopted & 
         keep_left_uncensored &
         !is.na(sex)) |> 
  mutate(js = as.numeric(keep_survived_to_24m), 
         across(c(sex, birth_clan, mother, birth_year), as.factor)) |> 
  select(ID, Fhat3, Fped, js,
         keep_sextets_and_septets, keep_not_too_negative_Fgrm,
         sex, mother_rank_ord, clansize_all, birth_clan, mother, birth_year) |> 
  drop_na(ID, sex, mother_rank_ord, clansize_all, birth_clan, mother, birth_year) -> y4_juv

apply(y4_juv, 2, \(x) sum(is.na(x))) #check which columns contain NAs (Fhat3 and keep_not_too_negative_Fgrm)

nrow(y4_juv) #2093

## Fgrm data for juveniles ----
y4_juv |>
  select(-keep_sextets_and_septets) |> 
  filter(keep_not_too_negative_Fgrm) |> 
  drop_na() -> y4_juv_Fgrm

nrow(y4_juv_Fgrm) # 993
table(y4_juv_Fgrm$sex) # 497 females + 496 males
signif(mean(y4_juv_Fgrm$js), digits = 2) # 0.69

## Fped data for juveniles ----
y4_juv |>
  select(-Fhat3, -keep_not_too_negative_Fgrm) |> 
  filter(keep_sextets_and_septets) |> 
  drop_na() -> y4_juv_Fped

nrow(y4_juv_Fped) # 1254
table(y4_juv_Fped$sex) # 622 females + 632 males
signif(mean(y4_juv_Fped$js), digits = 2) # 0.63

## Consensus Fgrm data used for predictions ----
#note: the predictors need to be complete but not the response vars nor Fhat3
complete_table_Fhat3 |> 
  select(sex, mother_rank_ord, clansize_all, birth_clan, ID, mother, birth_year) |>
  drop_na() -> data_for_predictions
nrow(data_for_predictions) # 2393


# Fitting the models -----------------------------------------------------------

## LRS ----

### LRS Fgrm main models ----
y10_model_lrs_fgrm <- fitme(LRS ~ Fhat3 + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year),
                            family = negbin2,
                            data = y10_adult_Fgrm)
gof(y10_model_lrs_fgrm) # Note: a perfect goodness of fit would lead to W = 1 (and p ~ 1)
summary(y10_model_lrs_fgrm, details = TRUE)
saveRDS(y10_model_lrs_fgrm, file = "data/fitted_models/y10_model_lrs_fgrm.RDS")


### LRS Fgrm IBD x E model ----
y10_model_lrs_fgrm_IBxE <- fitme(LRS ~ Fhat3 * (sex + mother_rank_ord + clansize_all) + 
                                 (1 + Fhat3|birth_clan) + (1|mother) + (1|birth_year),
                                 family = negbin2,
                                 data = y10_adult_Fgrm)
gof(y10_model_lrs_fgrm_IBxE)
summary(y10_model_lrs_fgrm_IBxE, details = TRUE)
saveRDS(y10_model_lrs_fgrm_IBxE, file = "data/fitted_models/y10_model_lrs_fgrm_IBxE.RDS")


### LRS Fgrm IBD x E model no random slope ----
y10_model_lrs_fgrm_IBxE_norandslope <- update(y10_model_lrs_fgrm_IBxE, . ~ . - (1 + Fhat3|birth_clan) + (1|birth_clan),
                                              method = "REML")
gof(y10_model_lrs_fgrm_IBxE_norandslope)
summary(y10_model_lrs_fgrm_IBxE_norandslope, details = TRUE)

### LRS Fgrm SI models ----

#### Year threshold alternatives
y8_model_lrs_fgrm <- update(y10_model_lrs_fgrm, data = y8_adult_Fgrm)
gof(y8_model_lrs_fgrm)
summary(y8_model_lrs_fgrm, details = TRUE)
saveRDS(y8_model_lrs_fgrm, file = "data/fitted_models/y8_model_lrs_fgrm.RDS")

y12_model_lrs_fgrm <- update(y10_model_lrs_fgrm, data = y12_adult_Fgrm)
gof(y12_model_lrs_fgrm)
summary(y12_model_lrs_fgrm, details = TRUE)
saveRDS(y12_model_lrs_fgrm, file = "data/fitted_models/y12_model_lrs_fgrm.RDS")


#### Animal model
#### note: fit took 35 min
if (file.exists("data/fitted_models/y10_model_lrs_fgrm_am.RDS")) {
  y10_model_lrs_fgrm_am <- readRDS("data/fitted_models/y10_model_lrs_fgrm_am.RDS")
  } else {
  y10_model_lrs_fgrm_am <- update(y10_model_lrs_fgrm, . ~ . + corrMatrix(1|ID),
                                  corrMatrix = invAmat,
                                  verbose = c(TRACE = 1)) # to display fit progression
  saveRDS(y10_model_lrs_fgrm_am, file = "data/fitted_models/y10_model_lrs_fgrm_am.RDS")
  }


#### 2 components model
y10_adult_Fgrm$LRS_bin <- y10_adult_Fgrm$LRS > 0
y10_model_lrs_fgrm_0 <- fitme(LRS_bin ~ Fhat3 + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year),
                              family = binomial(link = "logit"),
                              data = y10_adult_Fgrm)
gof(y10_model_lrs_fgrm_0)
summary(y10_model_lrs_fgrm_0, details = TRUE)

y10_model_lrs_fgrm_non0 <- fitme(LRS ~ Fhat3 + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year),
                                 family = negbin2(trunc = 0),
                                 data = y10_adult_Fgrm[y10_adult_Fgrm$LRS > 0, ])
gof(y10_model_lrs_fgrm_non0)
summary(y10_model_lrs_fgrm_non0, details = TRUE)

### LRS Fped main models ----
y10_model_lrs_fped <- fitme(LRS ~ Fped + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year), 
                            family = negbin2,
                            data = y10_adult_Fped)
gof(y10_model_lrs_fped)
summary(y10_model_lrs_fped, details = TRUE)
saveRDS(y10_model_lrs_fped, file = "data/fitted_models/y10_model_lrs_fped.RDS")

### LRS Fped for same subset as Fgrm ----
y10_model_lrs_fped_geno <- fitme(LRS ~ Fped + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year), 
                                 family = negbin2,
                                 data = y10_adult_Fgrm)
gof(y10_model_lrs_fped_geno)
summary(y10_model_lrs_fped_geno, details = TRUE)


## Lifespan ----

### Lifespan Fgrm main models ----
y10_model_lifespan_fgrm <- fitme(rt_lifespan ~ Fhat3 + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year), 
                                 family = gaussian,
                                 data = y10_adult_Fgrm)
gof(y10_model_lifespan_fgrm)
summary(y10_model_lifespan_fgrm, details = TRUE)
saveRDS(y10_model_lifespan_fgrm, file = "data/fitted_models/y10_model_lifespan_fgrm.RDS")

### Lifespan Fgrm IBD x E model ----
y10_model_lifespan_fgrm_IBxE <- fitme(rt_lifespan ~ Fhat3 * (sex + mother_rank_ord + clansize_all) + 
                                      (1 + Fhat3|birth_clan) + (1|mother) + (1|birth_year),
                                      family = gaussian,
                                      data = y10_adult_Fgrm)
gof(y10_model_lifespan_fgrm_IBxE)
summary(y10_model_lifespan_fgrm_IBxE, details = TRUE)
saveRDS(y10_model_lifespan_fgrm_IBxE, file = "data/fitted_models/y10_model_lifespan_fgrm_IBxE.RDS")


### Lifespan Fgrm IBD x E model no random slope ----
y10_model_lifespan_fgrm_IBxE_norandslope <- update(y10_model_lifespan_fgrm_IBxE, . ~ . - (1 + Fhat3|birth_clan) + (1|birth_clan),
                                                   method = "REML")
gof(y10_model_lifespan_fgrm_IBxE_norandslope)
summary(y10_model_lifespan_fgrm_IBxE_norandslope, details = TRUE)

### Lifespan Fgrm SI models ----

#### year threshold alternatives
y8_model_lifespan_fgrm <- update(y10_model_lifespan_fgrm, data = y8_adult_Fgrm)
gof(y8_model_lifespan_fgrm) 
summary(y8_model_lifespan_fgrm, details = TRUE)
saveRDS(y8_model_lifespan_fgrm, file = "data/fitted_models/y8_model_lifespan_fgrm.RDS")

y12_model_lifespan_fgrm <- update(y10_model_lifespan_fgrm, data = y12_adult_Fgrm)
gof(y12_model_lifespan_fgrm) 
summary(y12_model_lifespan_fgrm, details = TRUE)
saveRDS(y12_model_lifespan_fgrm, file = "data/fitted_models/y12_model_lifespan_fgrm.RDS")


#### animal model
#note: fit took 6 min
if (file.exists("data/fitted_models/y10_model_lifespan_fgrm_am.RDS")) {
  y10_model_lifespan_fgrm_am <- readRDS("data/fitted_models/y10_model_lifespan_fgrm_am.RDS")
} else {
  y10_model_lifespan_fgrm_am <- update(y10_model_lifespan_fgrm, . ~ . + corrMatrix(1|ID),
                                       corrMatrix = invAmat,
                                       verbose = c(TRACE = 1)) # to display fit progression
  saveRDS(y10_model_lifespan_fgrm_am, file = "data/fitted_models/y10_model_lifespan_fgrm_am.RDS")
}


### Lifespan Fped main models----
y10_model_lifespan_fped <- fitme(rt_lifespan ~ Fped + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year), 
                                 family = gaussian,
                                 data = y10_adult_Fped)
gof(y10_model_lifespan_fped)
summary(y10_model_lifespan_fped, details = TRUE)
saveRDS(y10_model_lifespan_fped, file = "data/fitted_models/y10_model_lifespan_fped.RDS")

### Lifespan Fped for same subset as Fgrm ----
y10_model_lifespan_fped_geno <- fitme(rt_lifespan ~ Fped + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year), 
                                      family = gaussian,
                                      data = y10_adult_Fgrm)
gof(y10_model_lifespan_fped_geno)
summary(y10_model_lifespan_fped_geno, details = TRUE)


## Juvenile survival ----

### Juvenile survival Fgrm main model -----
y4_model_js_fgrm <- fitme(js ~ Fhat3 + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year),
                          family = binomial(link = "logit"),
                          data = y4_juv_Fgrm)
gof(y4_model_js_fgrm)
summary(y4_model_js_fgrm, details = TRUE)
saveRDS(y4_model_js_fgrm, file = "data/fitted_models/y4_model_js_fgrm.RDS")

### Juvenile survival Fgrm IBD x E model ----
y4_model_js_fgrm_IBxE <- fitme(js ~ Fhat3 * (sex + mother_rank_ord + clansize_all) +
                               (1 + Fhat3|birth_clan) + (1|mother) + (1|birth_year),
                               family = binomial(link = "logit"),
                               data = y4_juv_Fgrm)
gof(y4_model_js_fgrm_IBxE)
summary(y4_model_js_fgrm_IBxE, details = TRUE)
saveRDS(y4_model_js_fgrm_IBxE, file = "data/fitted_models/y4_model_js_fgrm_IBxE.RDS")

### Juvenile survival Fgrm IBD x E model no random slope ----
y4_model_js_fgrm_IBxE_norandslope <- update(y4_model_js_fgrm_IBxE, . ~ . - (1 + Fhat3|birth_clan) + (1|birth_clan),
                                            method = "REML")
gof(y4_model_js_fgrm_IBxE_norandslope)
summary(y4_model_js_fgrm_IBxE_norandslope, details = TRUE)

### Juvenile survival Fgrm SI model ----
#note: fit took 19 min
if (file.exists("data/fitted_models/y4_model_js_fgrm_am.RDS")) {
  y4_model_js_fgrm_am <- readRDS("data/fitted_models/y4_model_js_fgrm_am.RDS")
} else {
  y4_model_js_fgrm_am <- update(y4_model_js_fgrm, . ~ . + corrMatrix(1|ID),
                                corrMatrix = invAmat,
                                verbose = c(TRACE = 1)) # to display fit progression
  saveRDS(y4_model_js_fgrm_am, file = "data/fitted_models/y4_model_js_fgrm_am.RDS")
}


### Juvenile survival Fped main model ----
y4_model_js_fped <- fitme(js ~ Fped + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year),
                          family = binomial(link = "logit"),
                          data = y4_juv_Fped)
gof(y4_model_js_fped)
summary(y4_model_js_fped, details = TRUE)
saveRDS(y4_model_js_fped, file = "data/fitted_models/y4_model_js_fped.RDS")


### Juvenile survival Fped for same subset as Fgrm ----
y4_model_js_fped_geno <- fitme(js ~ Fped + sex + mother_rank_ord + clansize_all + (1|birth_clan) + (1|mother) + (1|birth_year),
                               family = binomial(link = "logit"),
                               data = y4_juv_Fgrm)
gof(y4_model_js_fped_geno)
summary(y4_model_js_fped_geno, details = TRUE)

## Saving all objects created above ----
#save.image(file = "data/fitted_models/datapred_and_model_fits.RData")


# LRT --------------------------------------------------------------------------

## Computing LRT of fixed effects for all models ----
#note: took 8 hours using 50 CPUs
(timing_LRT_fixef <- system.time({
  fixedef_lrs                <- compute_LRT(y10_model_lrs_fgrm,           boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_IBxE           <- compute_LRT(y10_model_lrs_fgrm_IBxE,      boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_fped           <- compute_LRT(y10_model_lrs_fped,           boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_y8             <- compute_LRT(y8_model_lrs_fgrm,            boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_y12            <- compute_LRT(y12_model_lrs_fgrm,           boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_fped_geno      <- compute_LRT(y10_model_lrs_fped_geno,      boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_0              <- compute_LRT(y10_model_lrs_fgrm_0,         boot.repl = nb_boot, cores = nb_cores)
  fixedef_lrs_non0           <- compute_LRT(y10_model_lrs_fgrm_non0,      boot.repl = nb_boot, cores = nb_cores)
  fixedef_lifespan           <- compute_LRT(y10_model_lifespan_fgrm,      boot.repl = nb_boot, cores = nb_cores)
  fixedef_lifespan_IBxE      <- compute_LRT(y10_model_lifespan_fgrm_IBxE, boot.repl = nb_boot, cores = nb_cores)
  fixedef_lifespan_fped      <- compute_LRT(y10_model_lifespan_fped,      boot.repl = nb_boot, cores = nb_cores)
  fixedef_lifespan_y8        <- compute_LRT(y8_model_lifespan_fgrm,       boot.repl = nb_boot, cores = nb_cores)
  fixedef_lifespan_y12       <- compute_LRT(y12_model_lifespan_fgrm,      boot.repl = nb_boot, cores = nb_cores)
  fixedef_lifespan_fped_geno <- compute_LRT(y10_model_lifespan_fped_geno, boot.repl = nb_boot, cores = nb_cores)
  fixedef_js                 <- compute_LRT(y4_model_js_fgrm,             boot.repl = nb_boot, cores = nb_cores)
  fixedef_js_IBxE            <- compute_LRT(y4_model_js_fgrm_IBxE,        boot.repl = nb_boot, cores = nb_cores)
  fixedef_js_fped            <- compute_LRT(y4_model_js_fped,             boot.repl = nb_boot, cores = nb_cores)
  fixedef_js_fped_geno       <- compute_LRT(y4_model_js_fped_geno,        boot.repl = nb_boot, cores = nb_cores)
}))

## Computing LRT of random slopes for all models ----
#note: took around 25 min using 50 CPUs
spaMM.options(nb_cores = nb_cores)
timing_LRT_randslope <- system.time({
  set.seed(123)
  LRS_LRTtests_IBxE_rr <- anova(update(y10_model_lrs_fgrm_IBxE, method = "REML"),
                                y10_model_lrs_fgrm_IBxE_norandslope,
                                boot.repl = nb_boot)
  set.seed(123)
  lifespan_LRTtests_IBxE_rr <- anova(update(y10_model_lifespan_fgrm_IBxE, method = "REML"),
                                     y10_model_lifespan_fgrm_IBxE_norandslope,
                                     boot.repl = nb_boot)
  set.seed(123)
  js_LRTtests_IBxE_rr <- anova(update(y4_model_js_fgrm_IBxE, method = "REML"),
                               y4_model_js_fgrm_IBxE_norandslope,
                               boot.repl = nb_boot)
})


# Inbreeding depression results ------------------------------------------------

## Fgrm effect as slopes ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y10_model_lrs_fgrm","y4_model_js_fgrm"),
                  LRT_names = c("fixedef_lifespan", "fixedef_lrs", "fixedef_js")) |> 
  filter(Parameter == "FGRM")
#               Trait Effect Parameter Estimate    SE   Chi2 p-value CensusYears                   Model   N
# 1     Log(Lifespan)  Fixed      FGRM   -2.658 0.914  8.325   0.006          10 y10_model_lifespan_fgrm 768
# 2               LRS  Fixed      FGRM   -9.036 2.173 15.777   0.001          10      y10_model_lrs_fgrm 768
# 3 Juvenile survival  Fixed      FGRM   -2.987 2.481  1.323   0.261           4        y4_model_js_fgrm 993

### Fgrm effect on lifespan ----
Fgrm_lifespan <- pdep_lifespan(y10_model_lifespan_fgrm,
                               focal_var = "Fhat3",
                               focal_values = c(0.25, 0),
                               newdata = data_for_predictions,
                               gam_inverse = fit_gam_inverse10y) 

Fgrm_lifespan |> 
  mutate(across(everything(), \(x) signif(x, digits = 3)))
#   focal_var pointp   low    up
# 1      0.25   1.74 0.771  4.03
# 2      0      5.26 4.17   6.46

signif(Fgrm_lifespan[2, "pointp"] - Fgrm_lifespan[1, "pointp"], digits = 3) # 3.52

### Fped effect on lifespan ----

Fped_lifespan <- pdep_lifespan(y10_model_lifespan_fped,
                               focal_var = "Fped",
                               focal_values = c(0.25, 0),
                               newdata = data_for_predictions,
                               gam_inverse = fit_gam_inverse10y) 

Fped_lifespan |> 
  mutate(across(everything(), \(x) signif(x, digits = 3)))
#   focal_var pointp   low    up
# 1      0.25   0.73 0.273  2.15
# 2      0      4.73 3.650  5.98

signif(Fped_lifespan[2, "pointp"] - Fped_lifespan[1, "pointp"], digits = 3) # 4

### Fgrm effect on LRS ----
Fgrm_LRS <- pdep_simple(y10_model_lrs_fgrm, 
                        focal_var = "Fhat3",
                        focal_values = c(0, 0.25),
                        newdata = data_for_predictions)
Fgrm_LRS |> 
  mutate(across(everything(), \(x) signif(x, digits = 3))) |> 
  select(focal_var, pointp, low, up)
#   focal_var pointp    low    up
# 1      0.25  0.241 0.0764 0.759
# 2      0.00  2.300 1.5600 3.410

signif(Fgrm_LRS[2, "pointp"] - Fgrm_LRS[1, "pointp"], digits = 3) # 2.06

### Fped effect on LRS ----
Fped_LRS <- pdep_simple(y10_model_lrs_fped, 
                        focal_var = "Fped",
                        focal_values = c(0, 0.25),
                        newdata = data_for_predictions)

Fped_LRS |> 
  mutate(across(everything(), \(x) signif(x, digits = 3))) |> 
  select(focal_var, pointp, low, up)
#   focal_var pointp    low    up
# 1      0.25  0.116 0.0173 0.781
# 2      0.00  2.170 1.4900 3.180

signif(Fped_LRS[2, "pointp"] - Fped_LRS[1, "pointp"], digits = 3) # 2.05

### Fgrm effect on js ----
Fgrm_js <- pdep_simple(y4_model_js_fgrm,
                       focal_var = "Fhat3",
                       focal_values = c(0, 0.25),
                       newdata = data_for_predictions)
Fgrm_js |> 
  mutate(across(everything(), \(x) signif(x, digits = 3))) |> 
  select(focal_var, pointp, low, up)
#   focal_var pointp   low    up
# 1      0.25  0.516 0.258 0.762
# 2      0.00  0.666 0.575 0.747

signif(Fgrm_js[2, "pointp"] - Fgrm_js[1, "pointp"], digits = 3) # 0.149

### Fped effect on js ----
Fped_js <- pdep_simple(y4_model_js_fped,
                       focal_var = "Fped",
                       focal_values = c(0, 0.25),
                       newdata = data_for_predictions)
Fped_js |> 
  mutate(across(everything(), \(x) signif(x, digits = 3))) |> 
  select(focal_var, pointp, low, up)
#   focal_var pointp   low    up
# 1      0.25  0.544 0.252 0.808
# 2      0.00  0.642 0.535 0.738

signif(Fped_js[2, "pointp"] - Fped_js[1, "pointp"], digits = 3) # 0.0982

## Other effects ----

### Social rank ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y10_model_lrs_fgrm","y4_model_js_fgrm"),
                  LRT_names = c("fixedef_lifespan", "fixedef_lrs", "fixedef_js")) |> 
  filter(Parameter == "Social rank")
#               Trait Effect   Parameter Estimate    SE   Chi2 p-value CensusYears                   Model   N
# 1     Log(Lifespan)  Fixed Social rank   -0.011 0.003 12.068   0.002          10 y10_model_lifespan_fgrm 768
# 2               LRS  Fixed Social rank   -0.030 0.006 23.734   0.001          10      y10_model_lrs_fgrm 768
# 3 Juvenile survival  Fixed Social rank   -0.047 0.007 43.149   0.001           4        y4_model_js_fgrm 993

### Sex ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y10_model_lrs_fgrm","y4_model_js_fgrm"),
                  LRT_names = c("fixedef_lifespan", "fixedef_lrs", "fixedef_js")) |> 
  filter(Parameter == "SexM")
#               Trait Effect Parameter Estimate    SE  Chi2 p-value CensusYears                   Model   N
# 1     Log(Lifespan)  Fixed      SexM   -0.126 0.058 4.600   0.028          10 y10_model_lifespan_fgrm 768
# 2               LRS  Fixed      SexM   -0.368 0.127 7.691   0.007          10      y10_model_lrs_fgrm 768
# 3 Juvenile survival  Fixed      SexM   -0.009 0.154 0.003   0.964           4        y4_model_js_fgrm 993

### Clan size ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y10_model_lrs_fgrm","y4_model_js_fgrm"),
                  LRT_names = c("fixedef_lifespan", "fixedef_lrs", "fixedef_js")) |> 
  filter(Parameter == "Clan size")
#               Trait Effect Parameter Estimate    SE   Chi2 p-value CensusYears                   Model   N
# 1     Log(Lifespan)  Fixed Clan size   -0.008 0.002 10.349   0.003          10 y10_model_lifespan_fgrm 768
# 2               LRS  Fixed Clan size   -0.017 0.004  8.107   0.004          10      y10_model_lrs_fgrm 768
# 3 Juvenile survival  Fixed Clan size   -0.002 0.004  0.178   0.686           4        y4_model_js_fgrm 993

## Inbreeding depression across social environments ----

### Fgrm:Sex ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm_IBxE", "y10_model_lrs_fgrm_IBxE","y4_model_js_fgrm_IBxE"),
                  LRT_names = c("fixedef_lifespan_IBxE", "fixedef_lrs_IBxE", "fixedef_js_IBxE")) |> 
  filter(Parameter %in% c("FGRM", "FGRM:SexM"))
#               Trait Effect Parameter Estimate    SE  Chi2 p-value CensusYears                        Model   N
# 1     Log(Lifespan)  Fixed      FGRM   -6.165 3.180    NA      NA          10 y10_model_lifespan_fgrm_IBxE 768
# 2     Log(Lifespan)  Fixed FGRM:SexM    0.944 1.841 0.241   0.615          10 y10_model_lifespan_fgrm_IBxE 768
# 3               LRS  Fixed      FGRM   -7.455 6.369    NA      NA          10      y10_model_lrs_fgrm_IBxE 768
# 4               LRS  Fixed FGRM:SexM   -6.317 4.440 2.019   0.142          10      y10_model_lrs_fgrm_IBxE 768
# 5 Juvenile survival  Fixed      FGRM  -11.854 7.024    NA      NA           4        y4_model_js_fgrm_IBxE 993
# 6 Juvenile survival  Fixed FGRM:SexM   -6.707 4.952 1.728   0.204           4        y4_model_js_fgrm_IBxE 993

### Fgrm:Clan size ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm_IBxE", "y10_model_lrs_fgrm_IBxE","y4_model_js_fgrm_IBxE"),
                  LRT_names = c("fixedef_lifespan_IBxE", "fixedef_lrs_IBxE", "fixedef_js_IBxE")) |> 
  filter(Parameter %in% c("FGRM:Clan size"))
#               Trait Effect      Parameter Estimate    SE  Chi2 p-value CensusYears                        Model   N
# 1     Log(Lifespan)  Fixed FGRM:Clan size    0.047 0.054 0.643   0.437          10 y10_model_lifespan_fgrm_IBxE 768
# 2               LRS  Fixed FGRM:Clan size    0.014 0.111 0.011   0.916          10      y10_model_lrs_fgrm_IBxE 768
# 3 Juvenile survival  Fixed FGRM:Clan size    0.172 0.114 2.127   0.159           4        y4_model_js_fgrm_IBxE 993

### Fgrm:Social rank ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm_IBxE", "y10_model_lrs_fgrm_IBxE","y4_model_js_fgrm_IBxE"),
                  LRT_names = c("fixedef_lifespan_IBxE", "fixedef_lrs_IBxE", "fixedef_js_IBxE")) |> 
  filter(Parameter %in% c("FGRM:Social rank"))
#               Trait Effect        Parameter Estimate    SE  Chi2 p-value CensusYears                        Model   N
# 1     Log(Lifespan)  Fixed FGRM:Social rank    0.032 0.092 0.119   0.740          10 y10_model_lifespan_fgrm_IBxE 768
# 2               LRS  Fixed FGRM:Social rank    0.122 0.221 0.295   0.566          10      y10_model_lrs_fgrm_IBxE 768
# 3 Juvenile survival  Fixed FGRM:Social rank    0.158 0.225 0.455   0.498           4        y4_model_js_fgrm_IBxE 993

### Random slope ----
VarCorr(y10_model_lrs_fgrm_IBxE)
LRS_LRTtests_IBxE_rr$rawBootLRT
ranef(y10_model_lrs_fgrm_IBxE)

VarCorr(y10_model_lifespan_fgrm_IBxE)
lifespan_LRTtests_IBxE_rr$rawBootLRT

VarCorr(y4_model_js_fgrm_IBxE)
js_LRTtests_IBxE_rr$rawBootLRT



# Tables -----------------------------------------------------------------------

## Table 1 ----
Table1 <- summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y10_model_lrs_fgrm","y4_model_js_fgrm"),
                            LRT_names = c("fixedef_lifespan", "fixedef_lrs", "fixedef_js"))

write.table(Table1, "figures_and_tables/Table1.txt", sep = "\t", quote = FALSE, row.names = FALSE) 


## Table 2 ----
Table2 <- summary_fit_table(model_names = c("y10_model_lifespan_fgrm_IBxE", "y10_model_lrs_fgrm_IBxE","y4_model_js_fgrm_IBxE"),
                            LRT_names = c("fixedef_lifespan_IBxE", "fixedef_lrs_IBxE", "fixedef_js_IBxE"))
#add 3 random slopes LRT not handled by summary_fit_table()
Table2[Table2$Trait == "LRS" & Table2$Parameter == "FGRM|Clan", c("Chi2", "p-value")] <- 
  round(LRS_LRTtests_IBxE_rr$rawBootLRT[, c("chi2_LR" , "p_value")], digit = 3)
Table2[Table2$Trait == "Juvenile survival" & Table2$Parameter == "FGRM|Clan", c("Chi2", "p-value")] <- 
  round(js_LRTtests_IBxE_rr$rawBootLRT[, c("chi2_LR" , "p_value")], digit = 3)
Table2[Table2$Trait == "Log(Lifespan)" & Table2$Parameter == "FGRM|Clan", c("Chi2", "p-value")] <- 
  round(lifespan_LRTtests_IBxE_rr$rawBootLRT[, c("chi2_LR" , "p_value")], digit = 3)

write.table(Table2, "figures_and_tables/Table2.txt", sep = "\t", quote = FALSE, row.names = FALSE) 


## Table 3 ----
Table3 <- summary_fit_table(model_names = c("y10_model_lifespan_fped", "y10_model_lrs_fped","y4_model_js_fped"),
                            LRT_names = c("fixedef_lifespan_fped", "fixedef_lrs_fped", "fixedef_js_fped"))

write.table(Table3, "figures_and_tables/Table3.txt", sep = "\t", quote = FALSE, row.names = FALSE) 


## Table S3 (output from double modelling approach for LRS: proba to reproduce, and zero truncated LRS) ----
TableS3 <- summary_fit_table(model_names = c("y10_model_lrs_fgrm", "y10_model_lrs_fgrm_0","y10_model_lrs_fgrm_non0"),
                                        LRT_names = c("fixedef_lrs", "fixedef_lrs_0", "fixedef_lrs_non0"))
TableS3 |> 
  filter(Parameter == "FGRM")

write.table(TableS3, "figures_and_tables/TableS3.txt", sep = "\t", quote = FALSE, row.names = FALSE) 




## Table S4 (Fgrm with VA) ----
TableS4 <- summary_fit_table(model_names = c("y10_model_lifespan_fgrm_am", "y10_model_lrs_fgrm_am","y4_model_js_fgrm_am"))

write.table(TableS4, "figures_and_tables/TableS4.txt", sep = "\t", quote = FALSE, row.names = FALSE) 



## Table S5 (Fped for same dataset as Fgrm)  ----
TableS5 <- summary_fit_table(model_names = c("y10_model_lifespan_fped_geno", "y10_model_lrs_fped_geno","y4_model_js_fped_geno"),
                             LRT_names = c("fixedef_lifespan_fped_geno", "fixedef_lrs_fped_geno", "fixedef_js_fped_geno"))
if (!file.exists("figures_and_tables/TableS5.txt")) {
  write.table(TableS5, "figures_and_tables/TableS5.txt", sep = "\t", quote = FALSE, row.names = FALSE) 
}

## Table comparing Fgrm estimates in the present or absence of VA ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y10_model_lifespan_fgrm_am",
                                  "y10_model_lrs_fgrm", "y10_model_lrs_fgrm_am",
                                  "y4_model_js_fgrm", "y4_model_js_fgrm_am")) |> 
  filter(Parameter == "FGRM") |> 
  mutate(VA = if_else(grepl("_am", Model), "With VA", "WithoutVA"), .after = "Parameter") -> TableVA

write.table(TableVA, "figures_and_tables/TableVA.txt", sep = "\t", quote = FALSE, row.names = FALSE) 



## Table comparing Fgrm estimates for different cohort cut-offs ----
summary_fit_table(model_names = c("y10_model_lifespan_fgrm", "y8_model_lifespan_fgrm", "y12_model_lifespan_fgrm",
                                  "y10_model_lrs_fgrm", "y8_model_lrs_fgrm", "y12_model_lrs_fgrm")) |> 
  filter(Parameter == "FGRM") -> TableCohorts

write.table(TableCohorts, "figures_and_tables/TableCohorts.txt", sep = "\t", quote = FALSE, row.names = FALSE) 


# Saving all -------------------------------------------------------------------
#save.image(file = "data/fitted_models/full_modelrun.RData")
