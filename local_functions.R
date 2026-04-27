# no need to run this script outside of the 00_main_script.R
# Helper functions

# Functions for handling IDs in plink ----
get_ID_plink <- function(x) {
  ifelse(grepl("CroCro", x$FID),
         paste0(x$FID, "_", x$IID),
         x$FID)
}

fetch_FID_plink <- function(ID) {
  ifelse(grepl("CroCro", ID), "CroCro", ID)
}

fetch_IID_plink <- function(ID) {
  ifelse(grepl("CroCro", ID), sapply(strsplit(ID, split = "_"), \(x) x[2]), ID)
}

complete_IID_CroCro <- function(ID) {
  ifelse(grepl("-", ID), ID, paste0("CroCro_", ID))
}

# Functions for modelling ----
compute_LRT <- function(fit, boot.repl = 1000, seed = 123, cores = 2) {
  spaMM.options(nb_cores = cores)
  set.seed(seed) # initialize the random generator for parametric bootstrap
  fixedef <- names(fixef(fit))[-1] # create vector of fixed effects to test
  if (any(fixedef == "Fhat3:sexmale") | any(fixedef == "Fped:sexmale")) {
    # remove main term test since meaningless when interaction is kept
    fixedef <- fixedef[fixedef %in% c("Fhat3:sexmale", "Fhat3:mother_rank_ord", "Fhat3:clansize_all",
                                      "Fped:sexmale", "Fped:mother_rank_ord", "Fped:clansize_all")]
  }
  fixedef[fixedef == "sexmale"]       <- "sex"       # for easier model updating later
  fixedef[fixedef == "Fhat3:sexmale"] <- "Fhat3:sex" # for easier model updating later
  fixedef[fixedef == "Fped:sexmale"]  <- "Fped:sex"  # for easier model updating later
  
  LRTtests <- list()
  for (i in seq_along(fixedef)) {
    LRTtests[[i]] <- anova(fit,
                           update(fit, as.formula(paste(". ~ . -", noquote(fixedef[[i]])))),
                           boot.repl = boot.repl)
  }
  names(LRTtests) <- fixedef
  message("done")
  LRTtests
}

extract_spamm_fixef <- function(model, model_name) {
  stopifnot(inherits(model, "HLfit"))
  summary_model <- summary(model, verbose = FALSE)
  coefs <- as.data.frame(summary_model$beta_table)
  coefs$Parameter <- rownames(coefs)
  rownames(coefs) <- NULL
  coefs$Effect <- "Fixed"
  coefs$Model <- model_name
  coefs$N <- nrow(model$data)
  coefs
}


extract_spamm_rand <- function(model, model_name) {
  stopifnot(inherits(model, "HLfit"))
  rand <- VarCorr(model)[, c("Group", "Term", "Variance")]
  rand |> 
    mutate(Parameter = case_when(Term == "(Intercept)" ~ Group,
                                 Term != "(Intercept)" ~ paste(Term, Group, sep = "|"),
                                 TRUE ~ NA),
           Effect = "Random",
           Model = !!model_name,
           N = nrow(!!model$data)) |> 
    select(Model, Effect, Parameter, Estimate = Variance, N)
}

extract_spamm_LRT <- function(LRT_list, model_name) {
  stopifnot(inherits(LRT_list[[1]], "fixedLRT"))
  res <- t(sapply(LRT_list, \(LRT) unlist(LRT$rawBootLRT)))
  res <- as.data.frame(res)
  res$Parameter <- rownames(res)
  res$Parameter[res$Parameter == "sex"] <- "sexmale"
  res$Parameter[res$Parameter == "Fhat3:sex"] <- "Fhat3:sexmale"
  rownames(res) <- NULL
  res$Model <- model_name
  res
}

summary_fit_table <- function(model_names, LRT_names = character(0)) {
  models <- setNames(lapply(model_names, get), model_names)
  if (length(LRT_names) > 0) {
    LRTobj <- setNames(lapply(LRT_names, get), model_names) # note: use model names not LRT names for second arg 
  }
  
  res <- bind_rows(lapply(model_names, function(name) {
    if (length(LRT_names) > 0) {
      estimates <- full_join(extract_spamm_fixef(models[[name]], name),
                             extract_spamm_LRT(LRTobj[[name]], name),
                             by = c("Parameter", "Model"))
    } else {
      estimates <- extract_spamm_fixef(models[[name]], name)
      estimates$chi2_LR <- NA
      estimates$p_value <- NA
    }
    random <- extract_spamm_rand(models[[name]], name)
    results <- full_join(estimates, random,
                         by = c("Estimate", "Parameter", "Effect", "Model", "N"))
    results |> 
      mutate(CensusYears = case_when(grepl("y4", name) ~ 4,
                                     grepl("y8", name) ~ 8,
                                     grepl("y10", name) ~ 10,
                                     grepl("y12", name) ~ 12,
                                     TRUE ~ NA),
             CensusYears = factor(CensusYears))  |> 
      select(Model, Effect, Parameter, Estimate, "Cond. SE", chi2_LR, p_value, N, CensusYears) -> results
  }))
  res |> 
    mutate(Trait = case_when(grepl("lrs|LRS", Model) ~ "LRS",
                             grepl("js", Model) ~ "Juvenile survival",
                             
                             grepl("lifespan", Model) ~ "Log(Lifespan)",
                             TRUE ~ NA),
           Trait = factor(Trait, levels = c("LRS", "Juvenile survival", "Log(Lifespan)")),
           Parameter = case_when(grepl("(Intercept)", Parameter) ~ "Intercept",
                                 grepl("^Fhat3:sexmale$", Parameter) ~ "FGRM:SexM",
                                 grepl("^Fhat3:mother_rank_ord$", Parameter) ~ "FGRM:Social rank",
                                 grepl("^Fhat3:clansize_all$", Parameter) ~ "FGRM:Clan size",
                                 grepl("^Fhat3$", Parameter) ~ "FGRM",
                                 grepl("^Fped:sexmale$", Parameter) ~ "FPED:SexM",
                                 grepl("^Fped:mother_rank_ord$", Parameter) ~ "FPED:Social rank",
                                 grepl("^Fped:clansize_all$", Parameter) ~ "FPED:Clan size",
                                 grepl("^Fped$", Parameter) ~ "FPED",
                                 grepl("^sexmale", Parameter) ~ "SexM",
                                 grepl("^mother_rank_ord:clansize_all$", Parameter) ~ "Social rank:Clan size",
                                 grepl("^mother_rank_ord$", Parameter) ~ "Social rank",
                                 grepl("^clansize_all$", Parameter) ~ "Clan size",
                                 grepl("Fhat3|birth_clan", Parameter, fixed = TRUE) ~ "FGRM|Clan",
                                 grepl("Fped|birth_clan", Parameter, fixed = TRUE) ~ "FPED|Clan",
                                 grepl("^birth_clan$", Parameter) ~ "Clan",
                                 grepl("^mother$", Parameter) ~ "Mother",
                                 grepl("^birth_year$", Parameter) ~ "Birth year",
                                 grepl("^ID$", Parameter) ~ "Individual",
                                 grepl("Residual", Parameter) ~ "Residual variance",
                                 TRUE ~ NA)) |> 
    mutate(across(where(is.numeric), \(x) round(x, digits = 3))) |> 
    select(Trait, Effect, Parameter, Estimate, SE = 'Cond. SE', Chi2 = chi2_LR, `p-value` = p_value, CensusYears, Model, N)
}

pdep_simple <- function(model, focal_var, focal_values, newdata) {
  newdata[[focal_var]] <- NA
  result <- pdep_effects(model,
                         focal_var    = focal_var,
                         focal_values = focal_values,
                         intervals    = "fixefVar",
                         newdata      = newdata)
  result <- result[order(result$focal_var, decreasing = TRUE), ]
  rownames(result) <- NULL
  result
}

pdep_lifespan <- function(model,
                          focal_var,
                          focal_values,
                          newdata,
                          gam_inverse) {
  #note: more complex than for LRS or juvenile survival since we cannot rely on
  #pdep_effects() due to the transformation of the response var which must occur
  #before the averaging
  rows <- lapply(focal_values, \(fval) {
    newdata[[focal_var]] <- fval
    pred_obj <- predict(model, newdata = newdata, intervals = "fixefVar")
    pred <- cbind(pred = pred_obj[, 1], attr(pred_obj, "intervals"))
    pred_untransformed <- apply(pred, 2, \(x) 
      mgcv::predict.gam(gam_inverse, newdata = data.frame(rt_lifespan = x)))
    c(focal_var = fval, colMeans(pred_untransformed))
  })
  result <- do.call(rbind, rows)
  result <- as.data.frame(result)
  result[, c("focal_var", "pred", "fixefVar_0.025", "fixefVar_0.975")] |>
    setNames(c("focal_var", "pointp", "low", "up"))
}

pdep_by_group <- function(group_var,
                          group_values,
                          newdata,
                          focal_var,
                          focal_values,
                          model_simple_js,
                          model_simple_lrs,
                          model_lifespan,
                          gam_inverse) {
  
  result_list <- lapply(group_values, \(grp) {
    newdata[[group_var]] <- grp
    
    list(js       = cbind(pdep_simple(model_simple_js,
                                      focal_var    = focal_var,
                                      focal_values = focal_values,
                                      newdata      = newdata),
                          setNames(data.frame(grp), group_var)),
         lifespan = cbind(pdep_lifespan(model_lifespan,
                                        focal_var    = focal_var,
                                        focal_values = focal_values,
                                        newdata      = newdata,
                                        gam_inverse  = gam_inverse),
                          setNames(data.frame(grp), group_var)),
         LRS      = cbind(pdep_simple(model_simple_lrs,
                                      focal_var    = focal_var,
                                      focal_values = focal_values,
                                      newdata      = newdata),
                          setNames(data.frame(grp), group_var)))})
  
  lapply(do.call(Map, c(list(list), result_list)), \(x) do.call(rbind, x))
}

# Plotting functions and objects -----------------------------------------------

## common objects ----
theme_homemade <- function(...) {
  theme_classic(...) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 18, margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        aspect.ratio = 1,
        plot.margin = unit(c(1, 1, 1, 1), "lines"))
}

clan_colors       <- c("#85C2FF", "#ffd937", "#AAE4D6", "#ff5757", "#6ED161", "#c0c0c0", "#FFFEE5", "#E38BFF", "#333333")
clan_colors_solid <- c("#0000FF", "#ff9900", "#009999", "#ff0000", "#339966", "#303030", "#c0c0c0", "#9933ff", "#333333")

## helpers to draw fig 1 ----
create_fig_Fped.vs.Fhat3 <- function(data, ylim = c(-0.16, 0.32), loess_span = 0.9) {
  
  data_outliers <- filter(data, !keep_not_too_negative_Fgrm)
  data_no_outliers <- filter(data, keep_not_too_negative_Fgrm)
  
  ggplot(data = data_no_outliers, 
         aes(x = Fped, y = Fhat3)) + ## to flag the outliers < -0.1
    geom_point(shape = 8, colour = "red", alpha = 0.8,
               data = data_outliers) +
    geom_point(alpha = 0.3, size = 2) +
    geom_segment(x = 0, xend = 0.25, y = 0, yend = 0.25, color = 'orange', linewidth = 1) +
    geom_smooth(method = "loess", se = FALSE, 
                linewidth = 0.8, color = "blue", alpha = 1/2,
                span = loess_span) +
    geom_line(aes(y = stage(start = Fhat3, after_stat = ymin)),
              stat = "smooth", method = "loess", lty = "dashed", alpha = 1/2, 
              color = "blue", span = loess_span) +
    geom_line(aes(y = stage(start = Fhat3, ymax)),
              stat = "smooth", method = "loess", lty = "dashed", alpha = 1/2, 
              color = "blue", span = loess_span) +
    theme_homemade() +
    labs(x = expression(F[PED]),
         y = expression(F[GRM])) +
    xlim(0, 0.25) +
    ylim(ylim[1], ylim[2]) +
    scale_color_manual(values = c(`FALSE` = 2, `TRUE` = 1), guide = "none") +
    annotation_custom(grid::grobTree(grid::textGrob(label = paste0("N = ", nrow(data)),
                                                    just = "right",
                                                    x = 0.9, y = 0.1, gp = grid::gpar(fontsize = 12))))
}

create_fig_F.vs.pred <- function(data, var_F, var_pred, xlab = "", ylim = c(-0.16, 0.32), geom = "point") {
  
  data$y <- data[[var_F]]
  data$x <- data[[var_pred]]
  data_outliers <- filter(data, !keep_not_too_negative_Fgrm)
  data_no_outliers <- filter(data, keep_not_too_negative_Fgrm) 
  
  ggplot(data_no_outliers, 
         aes(x = x, y = y)) +
    geom_point(shape = 8, colour = "red", alpha = 0.8,
               data = data_outliers) +
    annotation_custom(grid::grobTree(grid::textGrob(label = paste0("N = ", nrow(data)),
                                                    just = "right",
                                                    x = 0.9, y = 0.1, gp = grid::gpar(fontsize = 12)))) -> gg
  
  if (geom == "point") {
    gg +
      geom_point(alpha = 0.25, size = 2) +
      geom_hline(yintercept = 0, color = 'orange', linewidth = 1) +
      geom_smooth(method = "loess", se = FALSE, 
                  linewidth = 0.8, color = "blue", alpha = 1/2) +
      geom_line(aes(y = stage(start = y, after_stat = ymin)),
                stat = "smooth", method = "loess", lty = "dashed", alpha = 1/2, 
                color = "blue") +
      geom_line(aes(y = stage(start = y, ymax)),
                stat = "smooth", method = "loess", lty = "dashed", alpha = 1/2, 
                color = "blue") -> gg
  } else if (geom == "box") {
    gg +
      geom_boxplot(aes(fill = x),
                   data = data |> filter(keep_not_too_negative_Fgrm),
                   outlier.alpha = 0.25,
                   width = 0.5) -> gg
  } else {
    stop("geom not recognised")
  }
  
  gg +
    ylim(ylim[1], ylim[2]) +
    scale_shape_manual(values = c(19, 17), guide = "none") +
    scale_color_manual(values = c(1, 2), guide = "none") +
    theme_homemade() +
    labs(x = xlab,
         y = ifelse(var_F == "Fped", expression(F[PED]), expression(F[GRM])))
}

## helper to draw fig 2 and related figs ----
create_pdep_plot <- function(plot_data,
                             raw_data,
                             x_raw,
                             y_raw,
                             n_model,
                             x_lab       = expression(F[GRM]),
                             n_label_x   = 0.9,
                             n_label_y   = 0.9,
                             jitter_h    = 0,
                             jitter_w    = 0,
                             point_alpha = 0.3,
                             point_size  = 2,
                             line_width = NULL,
                             group_var   = NULL,
                             group_points = TRUE,
                             colours     = NULL,
                             shapes      = NULL,
                             alphas      = NULL,
                             breaks = waiver(),
                             ylim = c(NA, NA)) {
  
  y_lab_lookup <- c(LRS      = "LRS",
                    js       = "Juvenile survival",
                    lifespan = "Lifespan (years)")
  
  y_lab <- y_lab_lookup[[y_raw]]
  
  n_label <- annotation_custom(
    grid::grobTree(grid::textGrob(label = paste0("N = ", nrow(n_model$data)),
                                  just  = "right",
                                  x     = n_label_x,
                                  y     = n_label_y,
                                  gp    = grid::gpar(fontsize = 12)))
  )
  
  jitter_aes <- aes(x = .data[[x_raw]], y = .data[[y_raw]])
  
  if (!is.null(group_var)) {

    if (!is.null(colours) & group_points) {
      jitter_aes |> 
        modifyList(aes(colour = factor(.data[[group_var]]))) -> jitter_aes
    }
    if (!is.null(shapes) & group_points) {
      jitter_aes |> 
        modifyList(aes(shape = factor(.data[[group_var]]))) -> jitter_aes
    }
  }
  
  line_aes <- if (!is.null(group_var)) {
    aes(colour    = factor(.data[[group_var]]),
        linewidth = factor(.data[[group_var]]),
        alpha     = factor(.data[[group_var]]))
  } else {
    aes()
  }
  
  colour_scale <- if (!is.null(group_var) && !is.null(colours))
    scale_colour_manual(values = colours, breaks = breaks)
  
  shape_scale <- if (!is.null(group_var) && !is.null(shapes))
    scale_shape_manual(values = shapes, breaks = breaks)
  
  alpha_scale <- if (!is.null(group_var) && !is.null(alphas))
    scale_alpha_manual(values = alphas, breaks = breaks)
  
  linewidth_scale <- if (!is.null(group_var) && !is.null(line_width))
    scale_linewidth_manual(values = line_width, breaks = breaks)
  
  labs <- labs(x = x_lab, y = y_lab)
  if (!is.null(colours)) {
    labs |> 
      modifyList(labs(colour = "")) -> labs
  }
  if (!is.null(shapes)) {
    labs |> 
      modifyList(labs(shape = "")) -> labs
  }
  if (!is.null(alphas)) {
    labs |> 
      modifyList(labs(alpha = "")) -> labs
  }
  if (!is.null(line_width)) {
    labs |> 
      modifyList(labs(linewidth = "")) -> labs
  }

  ggplot(plot_data) +
    aes(x = focal_var) +
    geom_jitter(jitter_aes,
                data   = raw_data,
                height = jitter_h,
                width  = jitter_w,
                alpha  = point_alpha,
                size   = point_size) +
    geom_line(line_aes |> modifyList(aes(y = pointp))) +
    geom_line(line_aes |> modifyList(aes(y = up)),  linetype = "dashed") +
    geom_line(line_aes |> modifyList(aes(y = low)), linetype = "dashed") +
    n_label +
    colour_scale +
    shape_scale +
    alpha_scale +
    linewidth_scale +
    labs +
    coord_cartesian(ylim = ylim) +
    theme_homemade()
}

convert_clan_name <- function(clans) {
  sapply(clans, \(clan) {
    switch(clan,
           A = "Airstrip",
           E = "Engitati",
           F = "Forest",
           L = "Lemala",
           M = "Munge",
           N = "Ngoitokitok",
           S = "Shamba",
           T = "Triangle")
  }) |> unname()
}
