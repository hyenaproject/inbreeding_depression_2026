# Create Figure 2 and related
# Please call '01_run_always.R' first, in case you want to run this script separately:
# source("01_run_always.R")

# Figure 2 ---------------------------------------------------------------------

## Data preparation/filtering ----

### load full data set ----
complete_table_Fhat3 <- read_csv("data/complete_table_Fhat3.csv")

### Load fitted models ----
fit_gam_inverse10y      <- readRDS("data/fitted_models/fit_gam_inverse10y.RDS")

y10_model_lifespan_fgrm <- readRDS("data/fitted_models/y10_model_lifespan_fgrm.RDS")
y10_model_lrs_fgrm      <- readRDS("data/fitted_models/y10_model_lrs_fgrm.RDS")
y4_model_js_fgrm        <- readRDS("data/fitted_models/y4_model_js_fgrm.RDS")

### Consensus Fgrm data used for predictions ----
complete_table_Fhat3 |> 
  select(sex, mother_rank_ord, clansize_all, birth_clan, ID, mother, birth_year) |>
  drop_na() -> data_for_predictions
nrow(data_for_predictions) # 2393

## Loading data used to fit the models ----
y10_adult_Fgrm <- y10_model_lifespan_fgrm$data
y4_juv_Fgrm <- y4_model_js_fgrm$data

## Predictions ----

Fhat3_min <- min(y10_adult_Fgrm$Fhat3, na.rm = TRUE)
Fhat3_max <- max(y10_adult_Fgrm$Fhat3, na.rm = TRUE)

### Predict partial dependence effects for lifespan ----
Fgrm_lifespan_for_plot <- pdep_lifespan(y10_model_lifespan_fgrm,
                                        focal_var = "Fhat3",
                                        focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                        newdata = data_for_predictions,
                                        gam_inverse = fit_gam_inverse10y)

### Predict partial dependence effects for LRS ----
Fgrm_LRS_for_plot <- pdep_simple(y10_model_lrs_fgrm,
                                 focal_var = "Fhat3",
                                 focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                 newdata = data_for_predictions)

### Predict partial dependence effects for juvenile survival ----
Fgrm_js_for_plot <- pdep_simple(y4_model_js_fgrm,
                                focal_var = "Fhat3",
                                focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                newdata = data_for_predictions)


## Plot ----

### Draw sub-plots ----

set.seed(123) # for consistent jittering
Figure_2_A <- create_pdep_plot(plot_data = Fgrm_js_for_plot,
                               raw_data  = y4_juv_Fgrm,
                               x_raw     = "Fhat3",
                               y_raw     = "js",
                               n_model   = y4_model_js_fgrm,
                               jitter_h  = 0.05,
                               n_label_y = 0.1,
                               point_alpha = 0.2)

Figure_2_B <- create_pdep_plot(plot_data = Fgrm_lifespan_for_plot,
                               raw_data  = y10_adult_Fgrm,
                               x_raw     = "Fhat3",
                               y_raw     = "lifespan",
                               n_model   = y10_model_lifespan_fgrm,
                               point_alpha = 0.2)

Figure_2_C <- create_pdep_plot(plot_data = Fgrm_LRS_for_plot,
                               raw_data  = y10_adult_Fgrm,
                               x_raw     = "Fhat3",
                               y_raw     = "LRS",
                               n_model   = y10_model_lrs_fgrm,
                               point_alpha = 0.2)

Figure2 <- Figure_2_A / Figure_2_B / Figure_2_C +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

## Save plot ----
ggsave("figures_and_tables/Figure2.pdf", plot = Figure2, width = 5, height = 12, dpi = 300)
ggsave("figures_and_tables/Figure2.png", plot = Figure2, width = 5, height = 12, dpi = 300)


# Figure S3 (same as Fig. 2BC using threshold 12yrs rather than Fgrm) --------------
y12_model_lifespan_fgrm <- readRDS("data/fitted_models/y12_model_lifespan_fgrm.RDS")
y12_model_lrs_fgrm      <- readRDS("data/fitted_models/y12_model_lrs_fgrm.RDS")
fit_gam_inverse8y      <- readRDS("data/fitted_models/fit_gam_inverse8y.RDS")
y12_adult_Fgrm <- y12_model_lifespan_fgrm$data

Fgrm_lifespan_for_plot8y <- pdep_lifespan(y12_model_lifespan_fgrm,
                                          focal_var = "Fhat3",
                                          focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                          newdata = data_for_predictions,
                                          gam_inverse = fit_gam_inverse8y)

Fgrm_LRS_for_plot8y <- pdep_simple(y12_model_lrs_fgrm,
                                   focal_var = "Fhat3",
                                   focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                   newdata = data_for_predictions)

FigureS3A <- create_pdep_plot(plot_data = Fgrm_lifespan_for_plot8y,
                              raw_data  = y12_adult_Fgrm,
                              x_raw     = "Fhat3",
                              y_raw     = "lifespan",
                              n_model   = y12_model_lifespan_fgrm,
                              point_alpha = 0.2)

FigureS3B <- create_pdep_plot(plot_data = Fgrm_LRS_for_plot8y,
                              raw_data  = y12_adult_Fgrm,
                              x_raw     = "Fhat3",
                              y_raw     = "LRS",
                              n_model   = y12_model_lrs_fgrm,
                              point_alpha = 0.2)

FigureS3 <- FigureS3A + FigureS3B +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

## Save plot ----
ggsave("figures_and_tables/FigureS3.pdf", plot = FigureS3, width = 10, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS3.png", plot = FigureS3, width = 10, height = 5, dpi = 300)


# Figure S4 (same as Fig. 2BC using threshold 8yrs rather than Fgrm) --------------
y8_model_lifespan_fgrm <- readRDS("data/fitted_models/y8_model_lifespan_fgrm.RDS")
y8_model_lrs_fgrm      <- readRDS("data/fitted_models/y8_model_lrs_fgrm.RDS")
fit_gam_inverse8y      <- readRDS("data/fitted_models/fit_gam_inverse8y.RDS")
y8_adult_Fgrm <- y8_model_lifespan_fgrm$data

Fgrm_lifespan_for_plot8y <- pdep_lifespan(y8_model_lifespan_fgrm,
                                          focal_var = "Fhat3",
                                          focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                          newdata = data_for_predictions,
                                          gam_inverse = fit_gam_inverse8y)

Fgrm_LRS_for_plot8y <- pdep_simple(y8_model_lrs_fgrm,
                                   focal_var = "Fhat3",
                                   focal_values = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                   newdata = data_for_predictions)

FigureS4A <- create_pdep_plot(plot_data = Fgrm_lifespan_for_plot8y,
                              raw_data  = y8_adult_Fgrm,
                              x_raw     = "Fhat3",
                              y_raw     = "lifespan",
                              n_model   = y8_model_lifespan_fgrm,
                              point_alpha = 0.2)

FigureS4B <- create_pdep_plot(plot_data = Fgrm_LRS_for_plot8y,
                              raw_data  = y8_adult_Fgrm,
                              x_raw     = "Fhat3",
                              y_raw     = "LRS",
                              n_model   = y8_model_lrs_fgrm,
                              point_alpha = 0.2)

FigureS4 <- FigureS4A + FigureS4B +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

## Save plot ----
ggsave("figures_and_tables/FigureS4.pdf", plot = FigureS4, width = 10, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS4.png", plot = FigureS4, width = 10, height = 5, dpi = 300)


# Figure S7 (same approach as Fig. 2 showing Fgrm effect per sex) --------------
y10_model_lifespan_fgrm_IBxE <- readRDS("data/fitted_models/y10_model_lifespan_fgrm_IBxE.RDS")
y10_model_lrs_fgrm_IBxE      <- readRDS("data/fitted_models/y10_model_lrs_fgrm_IBxE.RDS")
y4_model_js_fgrm_IBxE        <- readRDS("data/fitted_models/y4_model_js_fgrm_IBxE.RDS")

sex_to_do <- c("female", "male")

Fgrm_sex_for_plot <- pdep_by_group(group_var        = "sex",
                                   group_values     = sex_to_do,
                                   newdata          = data_for_predictions,
                                   focal_var        = "Fhat3",
                                   focal_values     = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                   model_simple_js  = y4_model_js_fgrm_IBxE,
                                   model_simple_lrs = y10_model_lrs_fgrm_IBxE,
                                   model_lifespan   = y10_model_lifespan_fgrm_IBxE,
                                   gam_inverse      = fit_gam_inverse10y)

set.seed(123) # for consistent jittering
FigureS7_A <- create_pdep_plot(plot_data   = Fgrm_sex_for_plot$js,
                                raw_data    = y4_juv_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "js",
                                n_model     = y4_model_js_fgrm_IBxE,
                                jitter_h    = 0.05,
                                n_label_y   = 0.1,
                                group_var   = "sex",
                                line_width  = c(0.6, 1),
                                colours     = c("#8624F5", "#1FC3AA"),
                                shapes      = c(1, 16),
                                alphas      = c(0.7, 1))

FigureS7_B <- create_pdep_plot(plot_data   = Fgrm_sex_for_plot$lifespan,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "lifespan",
                                n_model     = y10_model_lifespan_fgrm_IBxE,
                                group_var   = "sex",
                                line_width  = c(0.6, 1),
                                colours     = c("#8624F5", "#1FC3AA"),
                                shapes      = c(1, 16),
                                alphas      = c(0.7, 1))

FigureS7_C <- create_pdep_plot(plot_data   = Fgrm_sex_for_plot$LRS,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "LRS",
                                n_model     = y10_model_lrs_fgrm_IBxE,
                                group_var   = "sex",
                                line_width  = c(0.6, 1),
                                colours     = c("#8624F5", "#1FC3AA"),
                                shapes      = c(1, 16),
                                alphas      = c(0.7, 1),
                                ylim        = c(NA, 45))

FigureS7 <- FigureS7_A + FigureS7_B + FigureS7_C +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures_and_tables/FigureS7.pdf", plot = FigureS7, width = 14, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS7.png", plot = FigureS7, width = 14, height = 5, dpi = 300)


# Figure S8 (same approach as Fig. 2 showing Fgrm effect per rank) -------------

data_for_predictions |> 
  summarise(max_rank = max(mother_rank_ord), .by = "birth_clan") |> 
  filter(birth_clan != "U") |> # remove non-main clan
  summarise(max_rank = mean(max_rank)) |> # avg of max per main clan
  mutate(min_rank = 1) |> 
  as.numeric() |> 
  setNames(c("bottom ranking", "top ranking")) |> 
  rev() -> ranks_to_do
ranks_to_do

Fgrm_rank_for_plot <- pdep_by_group(group_var        = "mother_rank_ord",
                                    group_values     = ranks_to_do,
                                    newdata          = data_for_predictions,
                                    focal_var        = "Fhat3",
                                    focal_values     = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                    model_simple_js  = y4_model_js_fgrm_IBxE,
                                    model_simple_lrs = y10_model_lrs_fgrm_IBxE,
                                    model_lifespan   = y10_model_lifespan_fgrm_IBxE,
                                    gam_inverse      = fit_gam_inverse10y)

set.seed(123) # for consistent jittering
FigureS8_A <- create_pdep_plot(plot_data   = Fgrm_rank_for_plot$js,
                                raw_data    = y4_juv_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "js",
                                n_model     = y4_model_js_fgrm_IBxE,
                                jitter_h    = 0.05,
                                n_label_y   = 0.1,
                                group_var   = "mother_rank_ord",
                                group_points = FALSE,
                                line_width  = c(1, 0.6),
                                colours     = c("darkgoldenrod2", "coral"),
                                alphas      = c(1, 0.7),
                                breaks = ranks_to_do)

FigureS8_B <- create_pdep_plot(plot_data   = Fgrm_rank_for_plot$lifespan,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "lifespan",
                                n_model     = y10_model_lifespan_fgrm_IBxE,
                                group_var   = "mother_rank_ord",
                                group_points = FALSE,
                                line_width  = c(1, 0.6),
                                colours     = c("darkgoldenrod2", "coral"),
                                alphas      = c(1, 0.7),
                                breaks = ranks_to_do)

FigureS8_C <- create_pdep_plot(plot_data   = Fgrm_rank_for_plot$LRS,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "LRS",
                                n_model     = y10_model_lrs_fgrm_IBxE,
                                group_var   = "mother_rank_ord",
                                group_points = FALSE,
                                line_width  = c(1, 0.6),
                                colours     = c("darkgoldenrod2", "coral"),
                                alphas      = c(1, 0.7),
                                ylim        = c(NA, 45),
                                breaks = ranks_to_do)

FigureS8 <- FigureS8_A + FigureS8_B + FigureS8_C +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures_and_tables/FigureS8.pdf", plot = FigureS8, width = 14, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS8.png", plot = FigureS8, width = 14, height = 5, dpi = 300)


# Figure S9 (same approach as Fig. 2 showing Fgrm effect per clan size) --------

data_for_predictions |> 
  summarise(max_clansize = max(clansize_all),
            min_clansize = min(clansize_all),
            .by = "birth_clan") |> 
  filter(birth_clan != "U") |> # remove non-main clan
  summarise(max_clansize = mean(max_clansize),# avg per main clan
            min_clansize = mean(min_clansize)) |>
  as.numeric() |> 
  setNames(c("large clan", "small clan")) -> clansizes_to_do
clansizes_to_do
# large clan small clan 
#     96.750     19.375 

Fgrm_clansize_for_plot <- pdep_by_group(group_var        = "clansize_all",
                                        group_values     = clansizes_to_do,
                                        newdata          = data_for_predictions,
                                        focal_var        = "Fhat3",
                                        focal_values     = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                        model_simple_js  = y4_model_js_fgrm_IBxE,
                                        model_simple_lrs = y10_model_lrs_fgrm_IBxE,
                                        model_lifespan   = y10_model_lifespan_fgrm_IBxE,
                                        gam_inverse      = fit_gam_inverse10y)

set.seed(123) # for consistent jittering
FigureS9_A <- create_pdep_plot(plot_data   = Fgrm_clansize_for_plot$js,
                                raw_data    = y4_juv_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "js",
                                n_model     = y4_model_js_fgrm_IBxE,
                                jitter_h    = 0.05,
                                n_label_y   = 0.1,
                                group_var   = "clansize_all",
                                group_points = FALSE,
                                line_width  = c(1, 0.6),
                                colours     = c("#c46666", "#298c8c"),
                                alphas      = c(1, 0.7),
                                breaks = clansizes_to_do)

FigureS9_B <- create_pdep_plot(plot_data   = Fgrm_clansize_for_plot$lifespan,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "lifespan",
                                n_model     = y10_model_lifespan_fgrm_IBxE,
                                group_var   = "clansize_all",
                                group_points = FALSE,
                                line_width  = c(1, 0.6),
                                colours     = c("#c46666", "#298c8c"),
                                alphas      = c(1, 0.7),
                                breaks = clansizes_to_do)

FigureS9_C <- create_pdep_plot(plot_data   = Fgrm_clansize_for_plot$LRS,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "LRS",
                                n_model     = y10_model_lrs_fgrm_IBxE,
                                group_var   = "clansize_all",
                                group_points = FALSE,
                                line_width  = c(1, 0.6),
                                colours     = c("#c46666", "#298c8c"),
                                alphas      = c(1, 0.7),
                                ylim        = c(NA, 45),
                                breaks = clansizes_to_do)

FigureS9 <- FigureS9_A + FigureS9_B + FigureS9_C +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures_and_tables/FigureS9.pdf", plot = FigureS9, width = 14, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS9.png", plot = FigureS9, width = 14, height = 5, dpi = 300)


# Figure S10 (same approach as Fig. 2 showing Fgrm effect per clan size) --------

main_clans <- c("A", "E", "F", "L", "M", "N", "S", "T")

Fgrm_birthclan_for_plot <- pdep_by_group(group_var        = "birth_clan",
                                         group_values     = main_clans,
                                         newdata          = data_for_predictions,
                                         focal_var        = "Fhat3",
                                         focal_values     = seq(Fhat3_min, Fhat3_max, length.out = 25),
                                         model_simple_js  = y4_model_js_fgrm_IBxE,
                                         model_simple_lrs = y10_model_lrs_fgrm_IBxE,
                                         model_lifespan   = y10_model_lifespan_fgrm_IBxE,
                                         gam_inverse      = fit_gam_inverse10y)

set.seed(123) # for consistent jittering
FigureS10_A <- create_pdep_plot(plot_data   = Fgrm_birthclan_for_plot$js,
                                raw_data    = y4_juv_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "js",
                                n_model     = y4_model_js_fgrm_IBxE,
                                jitter_h    = 0.05,
                                n_label_y   = 0.1,
                                group_var   = "birth_clan",
                                group_points = TRUE,
                                colours     = clan_colors_solid,
                                line_width = rev(seq(0.5, 1, length = 8)),
                                alphas = rev(seq(0.5, 1, length = 8)),
                                breaks = setNames(main_clans, convert_clan_name(main_clans)))

FigureS10_B <- create_pdep_plot(plot_data   = Fgrm_birthclan_for_plot$lifespan,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "lifespan",
                                n_model     = y10_model_lifespan_fgrm_IBxE,
                                group_var   = "birth_clan",
                                group_points = TRUE,
                                colours     = clan_colors_solid,
                                line_width = rev(seq(0.5, 1, length = 8)),
                                alphas = rev(seq(0.5, 1, length = 8)),
                                breaks = setNames(main_clans, convert_clan_name(main_clans)))

FigureS10_C <- create_pdep_plot(plot_data   = Fgrm_birthclan_for_plot$LRS,
                                raw_data    = y10_adult_Fgrm,
                                x_raw       = "Fhat3",
                                y_raw       = "LRS",
                                n_model     = y10_model_lrs_fgrm_IBxE,
                                group_var   = "birth_clan",
                                group_points = TRUE,
                                colours     = clan_colors_solid,
                                line_width = rev(seq(0.5, 1, length = 8)),
                                alphas = rev(seq(0.5, 1, length = 8)),
                                breaks = setNames(main_clans, convert_clan_name(main_clans)))

FigureS10 <- FigureS10_A + FigureS10_B + FigureS10_C +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures_and_tables/FigureS10.pdf", plot = FigureS10, width = 14, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS10.png", plot = FigureS10, width = 14, height = 5, dpi = 300)


# Figure Extra (same as Fig. S10 but focus on LRS) ---------------------------------

Fgrm_birthclan_for_plot_coarser <- pdep_by_group(group_var        = "birth_clan",
                                                 group_values     = main_clans,
                                                 newdata          = data_for_predictions,
                                                 focal_var        = "Fhat3",
                                                 focal_values     = seq(Fhat3_min, Fhat3_max, length.out = 10),
                                                 model_simple_js  = y4_model_js_fgrm_IBxE,
                                                 model_simple_lrs = y10_model_lrs_fgrm_IBxE,
                                                 model_lifespan   = y10_model_lifespan_fgrm_IBxE,
                                                 gam_inverse      = fit_gam_inverse10y)

FigureS11 <- ggplot(Fgrm_birthclan_for_plot_coarser$LRS) +
  aes(x = focal_var) +
  geom_point(aes(y = log(pointp), colour = birth_clan, shape = birth_clan),
             size = 3) +
  geom_line(aes(y = log(pointp), colour = birth_clan), linetype = "dotted") +
  scale_shape_manual(values = 1:8,
                     breaks = setNames(main_clans, convert_clan_name(main_clans))) +
  scale_colour_manual(values = clan_colors_solid,
                      breaks = setNames(main_clans, convert_clan_name(main_clans))) +
  labs(colour = "", shape = "", y = "log LRS", x = expression(F[GRM])) +
  theme_homemade()


# Figure S11 (same approach as Fig. 2 using Fped rather than Fgrm) --------------
y10_model_lifespan_fped <- readRDS("data/fitted_models/y10_model_lifespan_fped.RDS")
y10_model_lrs_fped      <- readRDS("data/fitted_models/y10_model_lrs_fped.RDS")
y4_model_js_fped        <- readRDS("data/fitted_models/y4_model_js_fped.RDS")

y10_adult_Fped <- y10_model_lifespan_fped$data
y4_juv_Fped <- y4_model_js_fped$data

Fped_min <- min(y10_adult_Fped$Fped, na.rm = TRUE)
Fped_max <- max(y10_adult_Fped$Fped, na.rm = TRUE)

Fped_lifespan_for_plot <- pdep_lifespan(y10_model_lifespan_fped,
                                        focal_var = "Fped",
                                        focal_values = seq(Fped_min, Fped_max, length.out = 25),
                                        newdata = data_for_predictions,
                                        gam_inverse = fit_gam_inverse10y)

Fped_LRS_for_plot <- pdep_simple(y10_model_lrs_fped,
                                 focal_var = "Fped",
                                 focal_values = seq(Fped_min, Fped_max, length.out = 25),
                                 newdata = data_for_predictions)

Fped_js_for_plot <- pdep_simple(y4_model_js_fped,
                                focal_var = "Fped",
                                focal_values = seq(Fped_min, Fped_max, length.out = 25),
                                newdata = data_for_predictions)

set.seed(123) # for consistent jittering
FigureS11_A <- create_pdep_plot(plot_data = Fped_js_for_plot,
                                raw_data  = y4_juv_Fped,
                                x_raw     = "Fped",
                                y_raw     = "js",
                                n_model   = y4_model_js_fped,
                                x_lab       = expression(F[PED]),
                                jitter_h  = 0.05,
                                n_label_y = 0.1,
                                point_alpha = 0.2)

FigureS11_B <- create_pdep_plot(plot_data = Fped_lifespan_for_plot,
                                raw_data  = y10_adult_Fped,
                                x_raw     = "Fped",
                                y_raw     = "lifespan",
                                n_model   = y10_model_lifespan_fped,
                                x_lab       = expression(F[PED]),
                                point_alpha = 0.2)

FigureS11_C <- create_pdep_plot(plot_data = Fped_LRS_for_plot,
                                raw_data  = y10_adult_Fped,
                                x_raw     = "Fped",
                                y_raw     = "LRS",
                                n_model   = y10_model_lrs_fped,
                                x_lab       = expression(F[PED]),
                                point_alpha = 0.2)

FigureS11 <- FigureS11_A + FigureS11_B + FigureS11_C +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures_and_tables/FigureS11.pdf", plot = FigureS11, width = 11, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS11.png", plot = FigureS11, width = 11, height = 5, dpi = 300)


