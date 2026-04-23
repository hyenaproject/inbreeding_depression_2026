# Fig 1 and S1 created for the manuscript
# Please call '01_run_always.R' first, in case you want to inspect this script separately:
# source("01_run_always.R")

# Load data set
complete_table_Fhat3 <- read_csv("data/complete_table_Fhat3.csv")

# Figure 1 ----------------------------------------------------------------

## Fgrm vs Fped ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3), !is.na(Fped)) |> 
  create_fig_Fped.vs.Fhat3(loess_span = 0.85) -> Figure1_FgrmFped # increase span if no smoothed line

## Fgrm vs birth year ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3)) |> 
  create_fig_F.vs.pred(var_F = "Fhat3", var_pred = "birth_year",
                       xlab = "Birth year") -> Figure1_FgrmBirthYear

## Fgrm vs social rank ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3), !is.na(mother_rank_ord), keep_non_adopted) |> 
  create_fig_F.vs.pred(var_F = "Fhat3", var_pred = "mother_rank_ord", 
                       xlab = "Social rank at birth") -> Figure1_FgrmSocialRank

## Fgrm vs clan size ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3), !is.na(clansize_all)) |> 
  create_fig_F.vs.pred(var_F = "Fhat3", var_pred = "clansize_all", 
                       xlab = "Clan size") +
  scale_x_continuous(breaks = seq(0, 200, 25), limits = c(0, NA)) -> Figure1_FgrmClanSize

## Fgrm vs birth clan ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3)) |> 
  create_fig_F.vs.pred(var_F = "Fhat3", var_pred = "birth_clan", 
                       xlab = "Birth clan", geom = "box") +
  scale_fill_manual(values = clan_colors, guide = "none") -> Figure1_FgrmBirthClan

## Fgrm vs sex ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3), !is.na(sex)) |> 
  create_fig_F.vs.pred(var_F = "Fhat3", var_pred = "sex", 
                       xlab = "Sex", geom = "box") +
  scale_fill_manual(values = c("#8700f9", "#00c5a9"), guide = "none") -> Figure1_FgrmSex

## Combine to one plot ----
Figure1 <- (Figure1_FgrmFped + Figure1_FgrmSex + Figure1_FgrmSocialRank) /
           (Figure1_FgrmClanSize + Figure1_FgrmBirthClan + Figure1_FgrmBirthYear) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

# Save plot
ggsave("figures_and_tables/Figure1.pdf", plot = Figure1, width = 12, height = 8, dpi = 300)
ggsave("figures_and_tables/Figure1.png", plot = Figure1, width = 12, height = 8, dpi = 300)


# Figure S1 ---------------------------------------------------------------

## Fgrm vs Fped ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3), !is.na(Fped)) |> 
  create_fig_Fped.vs.Fhat3() -> FigureS1_FgrmFped

## Fped vs birth year ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped)) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "birth_year",
                       xlab = "Birth year") -> FigureS1_FpedBirthYear

## Fped vs social rank ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), !is.na(mother_rank_ord), keep_non_adopted) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "mother_rank_ord", 
                       xlab = "Social rank at birth") -> FigureS1_FpedSocialRank

## Fped vs clan size ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), !is.na(clansize_all)) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "clansize_all", 
                       xlab = "Clan size") +
  scale_x_continuous(breaks = seq(0, 200, 25), limits = c(0, NA)) -> FigureS1_FpedClanSize

## Fped vs birth clan ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped)) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "birth_clan", 
                       xlab = "Birth clan", geom = "box") +
  scale_fill_manual(values = clan_colors, guide = "none") -> FigureS1_FpedBirthClan

## Fped vs sex ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), !is.na(sex)) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "sex", 
                       xlab = "Sex", geom = "box") +
  scale_fill_manual(values = c("#8700f9", "#00c5a9"), guide = "none") -> FigureS1_FpedSex

## Combine to one plot ----
FigureS1 <- (FigureS1_FgrmFped + FigureS1_FpedSex + FigureS1_FpedSocialRank) /
  (FigureS1_FpedClanSize + FigureS1_FpedBirthClan + FigureS1_FpedBirthYear) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))

# Save plot
ggsave("figures_and_tables/FigureS1.pdf", plot = FigureS1, width = 12, height = 8, dpi = 300)
ggsave("figures_and_tables/FigureS1.png", plot = FigureS1, width = 12, height = 8, dpi = 300)


# Figure S2 ---------------------------------------------------------------

## Fgrm vs Fped ----
complete_table_Fhat3 |> 
  filter(!is.na(Fhat3), !is.na(Fped), keep_sextets_and_septets) |> 
  create_fig_Fped.vs.Fhat3() -> FigureS2_FgrmFped

## Fped vs birth year ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), keep_sextets_and_septets) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "birth_year",
                       xlab = "Birth year") -> FigureS2_FpedBirthYear

## Fped vs social rank ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), keep_sextets_and_septets, !is.na(mother_rank_ord), keep_non_adopted) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "mother_rank_ord", 
                       xlab = "Social rank at birth") -> FigureS2_FpedSocialRank

## Fped vs clan size ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), keep_sextets_and_septets, !is.na(clansize_all)) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "clansize_all", 
                       xlab = "Clan size") +
  scale_x_continuous(breaks = seq(0, 200, 25), limits = c(0, NA)) -> FigureS2_FpedClanSize

## Fped vs birth clan ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), keep_sextets_and_septets) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "birth_clan", 
                       xlab = "Birth clan", geom = "box") +
  scale_fill_manual(values = clan_colors, guide = "none") -> FigureS2_FpedBirthClan

## Fped vs sex ----
complete_table_Fhat3 |> 
  filter(!is.na(Fped), !is.na(sex), keep_sextets_and_septets) |> 
  create_fig_F.vs.pred(var_F = "Fped", var_pred = "sex", 
                       xlab = "Sex", geom = "box") +
  scale_fill_manual(values = c("#8700f9", "#00c5a9"), guide = "none") -> FigureS2_FpedSex

## Combine to one plot ----
FigureS2 <- (FigureS2_FgrmFped + FigureS2_FpedSex + FigureS2_FpedSocialRank) /
  (FigureS2_FpedClanSize + FigureS2_FpedBirthClan + FigureS2_FpedBirthYear) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(face = "bold"))


# Save plot
ggsave("figures_and_tables/FigureS2.pdf", plot = FigureS2, width = 12, height = 8, dpi = 300)
ggsave("figures_and_tables/FigureS2.png", plot = FigureS2, width = 12, height = 8, dpi = 300)
