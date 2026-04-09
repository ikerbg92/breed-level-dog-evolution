# ------------------------------------------------------------
# Posterior Predictive Checks (PPC)
# ------------------------------------------------------------
library(MCMCglmm)
library(ggplot2)
library(bayesplot)

ppc_mcmcglmm <- function(model, data, outcome_col, n_sims = 100, label = NULL) {
  # model: an MCMCglmm object
  # data: the dataframe used in the model
  # outcome_col: the name of the outcome variable (e.g. "OScv_N")
  # n_sims: number of posterior predictive datasets to simulate
  # label: optional text label for the model
  
  if (!(outcome_col %in% names(data))) {
    stop(paste("Outcome column", outcome_col, "not found in data."))
  }
  
  n_iter <- length(model$Sol[,1])
  sim_idx <- sample(1:n_iter, n_sims)
  y_obs <- data[[outcome_col]]
  y_rep <- matrix(NA, nrow = n_sims, ncol = length(y_obs))
  is_multinomial <- any(grepl("multinomial", model$family))
  X <- model$X
  vcv_names <- colnames(model$VCV)
  unit_col_pattern <- paste0(outcome_col, ".*units")
  unit_col_idx <- grep(unit_col_pattern, vcv_names)
  
  if (length(unit_col_idx) == 0 && !is_multinomial) {
    unit_col_idx <- grep("^units$", vcv_names)
  }
  
  for (i in 1:n_sims) {
    idx <- sim_idx[i]
    beta <- model$Sol[idx,]
    
    # Calculate fixed effects contribution
    # This version follows Achaz's original logic
    sol_names <- colnames(model$Sol)
    outcome_sol_idx <- grep(paste0("trait", outcome_col, "$"), sol_names)
    
    if (length(outcome_sol_idx) > 0) {
      mu <- model$Sol[idx, outcome_sol_idx[1]]
    } else {
      mu <- (X %*% beta)[,1]
    }
    
    if (is_multinomial) {
      n_vec <- if ("n" %in% names(data)) data$n else rep(round(median(y_obs, na.rm = TRUE)), length(y_obs))
      prob <- plogis(mu)
      y_rep[i,] <- rbinom(length(y_obs), size = round(n_vec), prob = prob)
    } else {
      if (length(unit_col_idx) == 0) {
        stop("Could not find residual variance column in model$VCV. Please check model structure or if outcome_col matches model traits.")
      }
      res_var <- model$VCV[idx, unit_col_idx[1]]
      y_rep[i,] <- rnorm(length(y_obs), mean = mu, sd = sqrt(res_var))
    }
  }
  
  # Visualize PPC
  ppc_data <- data.frame(y = y_obs, type = "Observed")
  for (j in 1:n_sims) {
    ppc_data <- rbind(ppc_data, data.frame(y = y_rep[j,], type = paste0("Simulated_", j)))
  }
  
  title_text <- if (is.null(label)) "Posterior Predictive Check" else paste("PPC:", label)
  
  ggplot(ppc_data, aes(x = y, group = type, color = type == "Observed")) +
    geom_density(alpha = 0.1, linewidth = 0.5, na.rm = TRUE) +
    scale_color_manual(values = c("gray", "red")) +
    theme_minimal() +
    labs(title = title_text,
         subtitle = paste("Outcome:", outcome_col, "| Observed (red) vs. Simulations (gray)"),
         x = outcome_col, y = "Density") +
    theme(legend.position = "none")
}

# ------------------------------------------------------------
# Assemble all PPC plots for one model
# ------------------------------------------------------------

assemble_all_ppcs <- function(model, data, categories, n_sims = 50) {
  plot_list <- list()
  
  for (cat in categories) {
    message("Plotting: ", cat)
    
    if (cat %in% names(data)) {
      p <- ppc_mcmcglmm(model, data, cat, n_sims = n_sims)
      p <- p + labs(title = cat, subtitle = NULL) +
        theme(plot.title = element_text(size = 10))
      plot_list[[cat]] <- p
    } else {
      message("Skipping ", cat, " - column not found in data.")
    }
  }
  
  if (length(plot_list) == 0) {
    stop("No valid category columns found in data.")
  }
  
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    install.packages("gridExtra")
  }
  
  do.call(gridExtra::grid.arrange, c(plot_list, ncol = 3))
}

# ------------------------------------------------------------
# Categories
# OS models use the 8 focal categories shown in Achaz's figures
# PP models must match the fitted response and therefore use PPothers_N
# ------------------------------------------------------------

categories_os <- c(
  "OScv_N", "OSendo_N", "OSgi_N", "OShem_N",
  "OSms_N", "OSneuro_N", "OSresp_N", "OSuro_N", "OSothers_N"
)

categories_pp <- c(
  "PPcongen_N", "PPdegen_N", "PPinfect_N", "PPinflam_N",
  "PPmetab_N", "PPneopl_N", "PPtraum_N", "PPothers_N"
)

# ------------------------------------------------------------
# PPCs already done by Achaz, but corrected to use the proper subsets
# ------------------------------------------------------------

# 1) OS no predictor
full_ppc_OS <- assemble_all_ppcs(OS_multinomial, fleming, categories_os)
ggsave("PPC_OS_All.pdf", full_ppc_OS, width = 12, height = 12)

# 4) PP no predictor
full_ppc_PP <- assemble_all_ppcs(PP_multinomial, fleming, categories_pp)
ggsave("PPC_PP_all.pdf", full_ppc_PP, width = 12, height = 12)

# 2) OS weight
full_ppc_OS_weight <- assemble_all_ppcs(OS_multinomial_weight, dfweight, categories_os)
ggsave("PPC_OS_weight.pdf", full_ppc_OS_weight, width = 12, height = 12)

# 3) OS growth
full_ppc_OS_growth <- assemble_all_ppcs(OS_multinomial_growth, dfgrowth, categories_os)
ggsave("PPC_OS_growth.pdf", full_ppc_OS_growth, width = 12, height = 12)

# ------------------------------------------------------------

# 5) PP weight
full_ppc_PP_weight <- assemble_all_ppcs(PP_multinomial_weight, dfweight, categories_pp)
ggsave("PPC_PP_weight.pdf", full_ppc_PP_weight, width = 12, height = 12)

# 6) PP growth
full_ppc_PP_growth <- assemble_all_ppcs(PP_multinomial_growth, dfgrowth, categories_pp)
ggsave("PPC_PP_growth.pdf", full_ppc_PP_growth, width = 12, height = 12)

# 7) OS reproductive investment
full_ppc_OS_bwkgxls <- assemble_all_ppcs(OS_multinomial_bwkgxls, dfbwkgxls, categories_os)
ggsave("PPC_OS_bwkgxls.pdf", full_ppc_OS_bwkgxls, width = 12, height = 12)

# 8) PP reproductive investment
full_ppc_PP_bwkgxls <- assemble_all_ppcs(PP_multinomial_bwkgxls, dfbwkgxls, categories_pp)
ggsave("PPC_PP_bwkgxls.pdf", full_ppc_PP_bwkgxls, width = 12, height = 12)

# 9) OS activity
full_ppc_OS_activity <- assemble_all_ppcs(OS_multinomial_activity, dfactivity, categories_os)
ggsave("PPC_OS_activity.pdf", full_ppc_OS_activity, width = 12, height = 12)

# 10) PP activity
full_ppc_PP_activity <- assemble_all_ppcs(PP_multinomial_activity, dfactivity, categories_pp)
ggsave("PPC_PP_activity.pdf", full_ppc_PP_activity, width = 12, height = 12)

# 11) OS aggressiveness
full_ppc_OS_aggressiveness <- assemble_all_ppcs(OS_multinomial_aggressiveness, dfaggressivenessCareau2010, categories_os)
ggsave("PPC_OS_aggressiveness.pdf", full_ppc_OS_aggressiveness, width = 12, height = 12)

# 12) PP aggressiveness
full_ppc_PP_aggressiveness <- assemble_all_ppcs(PP_multinomial_aggressiveness, dfaggressivenessCareau2010, categories_pp)
ggsave("PPC_PP_aggressiveness.pdf", full_ppc_PP_aggressiveness, width = 12, height = 12)

# 13) OS trainability
full_ppc_OS_trainability <- assemble_all_ppcs(OS_multinomial_trainability, dftrainability, categories_os)
ggsave("PPC_OS_trainability.pdf", full_ppc_OS_trainability, width = 12, height = 12)

# 14) PP trainability
full_ppc_PP_trainability <- assemble_all_ppcs(PP_multinomial_trainability, dftrainability, categories_pp)
ggsave("PPC_PP_trainability.pdf", full_ppc_PP_trainability, width = 12, height = 12)
