# All plots created for the manuscript
# Please call '01_run_always.R' first, in case you want to inspect this script separately:
# source("01_run_always.R")

# Figure S5 -------------------------------------------------------------------

## Load data sets 
TableCohorts <- read.table("figures_and_tables/TableCohorts.txt", sep = "\t", header = TRUE) 

FigureS5 <- ggplot(TableCohorts, 
                   aes(x = Trait, y = Estimate,
                       color = factor(CensusYears),
                       linewidth = factor(CensusYears))) +
  geom_point(position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = Estimate - SE, 
                    ymax = Estimate + SE, 
                    group = CensusYears), 
                width = 0.2, position = position_dodge(width = 0.3)) +
  scale_color_manual(breaks = c(8, 10, 12),
                     values = c("coral3","deepskyblue3","gold2"),
                     labels = c("2017/09/01","2015/09/01","2013/09/01")) +
  scale_linewidth_manual(values = c(0.5, 0.7, 1),
                         labels = c("2017/09/01","2015/09/01","2013/09/01")) +
  scale_x_discrete(labels = c("Lifespan", "LRS")) +
  labs(colour = "Most recent\n birthdate included:",
       linewidth = "Most recent\n birthdate included:",
       y = expression(paste("Fixed-effect estimate", ~F[GRM]~(+-SE))),
       x = NULL) +
  theme_homemade()

ggsave("figures_and_tables/FigureS5.pdf",
       plot = FigureS5, width = 7, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS5.png",
       plot = FigureS5, width = 7, height = 5, dpi = 300)


# Figure S6  ------------------------------------------------------------------

TableVA <- read.table("figures_and_tables/TableVA.txt", sep = "\t", header = TRUE) 

FigureS6 <- ggplot(TableVA) +
              aes(x = Trait, y = Estimate, colour = VA, linewidth = VA) +
              geom_point(position = position_dodge(width = 0.3)) +
              geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE, group = VA),
                            width = 0.2,
                            position = position_dodge(width = 0.3)) +
              geom_hline(yintercept = 0, linetype = "dashed", colour = "darkgrey") +
              scale_x_discrete(labels = c("Juvenile survival", "Lifespan", "LRS")) +
              scale_color_manual(values = c("mediumorchid3","goldenrod2"),
                                 labels = c(bquote("without V"[A]), bquote("with V"[A]))) +
              scale_linewidth_manual(values = c(0.5, 1),
                                     labels = c(bquote("without V"[A]), bquote("with V"[A]))) +
              labs(y = expression(paste("Fixed-effect estimate", ~F[GRM]~(+-SE))),
                   x = NULL, linewidth = NULL, colour = NULL) +
              theme_homemade()
  
ggsave("figures_and_tables/FigureS6.pdf", plot = FigureS6, width = 7, height = 5, dpi = 300)
ggsave("figures_and_tables/FigureS6.png", plot = FigureS6, width = 7, height = 5, dpi = 300)
