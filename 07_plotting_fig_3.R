# Create Figure 3
# Please call '01_run_always.R' first, in case you want to run this script separately:
# source("01_run_always.R")

# Data preparation/filtering ---------------------------------------------------

## load full data set ----
complete_table_Fhat3 <- read_csv("data/complete_table_Fhat3.csv")

## Load fitted models ----
y10_model_lifespan_fgrm <- readRDS("data/fitted_models/y10_model_lifespan_fgrm.RDS")
y10_model_lrs_fgrm      <- readRDS("data/fitted_models/y10_model_lrs_fgrm.RDS")
y4_model_js_fgrm        <- readRDS("data/fitted_models/y4_model_js_fgrm.RDS")

y10_model_lifespan_fped <- readRDS("data/fitted_models/y10_model_lifespan_fped.RDS")
y10_model_lrs_fped      <- readRDS("data/fitted_models/y10_model_lrs_fped.RDS")
y4_model_js_fped        <- readRDS("data/fitted_models/y4_model_js_fped.RDS")

fit_gam_inverse10y      <- readRDS("data/fitted_models/fit_gam_inverse10y.RDS")

## Consensus Fgrm data used for predictions ----
complete_table_Fhat3 |> 
  select(sex, mother_rank_ord, clansize_all, birth_clan, ID, mother, birth_year) |>
  drop_na() -> data_for_predictions
nrow(data_for_predictions) # 2393


# Compute predictions at F = 0 and F = 0.25 ------------------------------------

Fgrm_lifespan <- pdep_lifespan(y10_model_lifespan_fgrm,
                               focal_var = "Fhat3", 
                               focal_values = c(0, 0.25),
                               newdata = data_for_predictions,
                               gam_inverse = fit_gam_inverse10y)

Fped_lifespan <- pdep_lifespan(y10_model_lifespan_fped,
                               focal_var = "Fped",
                               focal_values = c(0, 0.25),
                               newdata = data_for_predictions,
                               gam_inverse = fit_gam_inverse10y)

Fgrm_LRS <- pdep_simple(y10_model_lrs_fgrm,
                        focal_var = "Fhat3",
                        focal_values = c(0, 0.25),
                        newdata = data_for_predictions)

Fped_LRS <- pdep_simple(y10_model_lrs_fped,
                        focal_var = "Fped", 
                        focal_values = c(0, 0.25),
                        newdata = data_for_predictions)

Fgrm_js  <- pdep_simple(y4_model_js_fgrm, 
                        focal_var = "Fhat3",
                        focal_values = c(0, 0.25),
                        newdata = data_for_predictions)

Fped_js  <- pdep_simple(y4_model_js_fped,
                        focal_var = "Fped", 
                        focal_values = c(0, 0.25),
                        newdata = data_for_predictions)

# Reshape to long format with % change standardisation -------------------------
# F = 0 is the reference (0% by construction); CIs are propagated using
# the F = 0 point estimate as the fixed denominator throughout.

make_row_long <- function(df, Fmetric, trait) {
  ref <- df[df$focal_var == 0, ]
  val <- df[df$focal_var == 0.25, ]
  rbind(
    data.frame(Fmetric, trait, Fvalue = "0",
               pointp     = 0,
               low        = (ref$low    - ref$pointp) / ref$pointp * 100,
               up         = (ref$up     - ref$pointp) / ref$pointp * 100,
               raw_pointp = ref$pointp),
    data.frame(Fmetric, trait, Fvalue = "0.25",
               pointp     = (val$pointp - ref$pointp) / ref$pointp * 100,
               low        = (val$low    - ref$pointp) / ref$pointp * 100,
               up         = (val$up     - ref$pointp) / ref$pointp * 100,
               raw_pointp = val$pointp)
  )
}

plot_df_long <- rbind(
  make_row_long(Fgrm_lifespan, "FGRM", "Lifespan (years)"),
  make_row_long(Fped_lifespan, "FPED", "Lifespan (years)"),
  make_row_long(Fgrm_LRS,      "FGRM", "LRS"),
  make_row_long(Fped_LRS,      "FPED", "LRS"),
  make_row_long(Fgrm_js,       "FGRM", "Juvenile survival"),
  make_row_long(Fped_js,       "FPED", "Juvenile survival")
)

plot_df_long |>
  mutate(trait  = factor(trait,  levels = c("Juvenile survival", "Lifespan (years)", "LRS")),
         Fmetric = factor(Fmetric, levels = c("FGRM", "FPED")),
         Fvalue  = factor(Fvalue,  levels = c("0", "0.25")),
         FF      = factor(interaction(Fmetric, Fvalue),
                          levels = c("FGRM.0", "FGRM.0.25", "FPED.0", "FPED.0.25")),
         x_pos   = as.numeric(trait) + 
           case_when(Fmetric == "FGRM" & Fvalue == "0"    ~ -0.25,
                     Fmetric == "FGRM" & Fvalue == "0.25" ~ -0.1,
                     Fmetric == "FPED" & Fvalue == "0"    ~  0.1,
                     Fmetric == "FPED" & Fvalue == "0.25" ~  0.25,
                     .default = NA),
         label_xpos = x_pos + c(-0.15, 0.15),
         label_ypos = pointp + c(3, 0)
  ) -> plot_df_long


# Plot -------------------------------------------------------------------------

labels_legend <- c("FGRM.0.25" = expression(F[GRM]),
                   "FPED.0.25" = expression(F[PED]))

Figure3 <- ggplot(plot_df_long |> filter(Fvalue == "0.25")) + 
  aes(x = x_pos, y = pointp, fill = FF, shape = FF) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_col(aes(alpha = FF)) +
  geom_errorbar(aes(ymin = low, ymax = up),
                width     = 0.05,
                linewidth = 0.3) +
  scale_fill_manual(values = c("FGRM.0.25" = "#56B4E9", "FPED.0.25" = "#E69F00"),
                      labels = labels_legend) +
  scale_x_continuous(breaks = seq_along(levels(plot_df_long$trait)),
                     labels = levels(plot_df_long$trait)) +
  scale_alpha_manual(values = c("FGRM.0.25" = 0.9, "FPED.0.25" = 0.5),
                     labels = labels_legend) +
  labs(x = NULL, y = "Fitness % change, relative to outbred",
       fill = "", shape = "", alpha = "") +
  theme_homemade()

ggsave("figures_and_tables/Figure3.pdf", plot = Figure3, width = 7, height = 5, dpi = 300)
ggsave("figures_and_tables/Figure3.png", plot = Figure3, width = 7, height = 5, dpi = 300)

