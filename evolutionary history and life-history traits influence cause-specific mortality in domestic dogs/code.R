
setwd()

# ------------------------------------------------------------
# Load libraries used across the workflow:
# - Data import/export (readxl, writexl, read.csv)
# - Data wrangling and reshaping (dplyr, tidyr, reshape2, stringr, tibble)
# - Plotting and diagnostic utilities (ggplot2, PerformanceAnalytics, coda)
# - Statistical helper packages (orthopolynom, jtools)
# ------------------------------------------------------------

library(writexl)
library(readxl)
library(orthopolynom)
library(PerformanceAnalytics)
library(dplyr) 
library(jtools)
library(tibble)
library(ggplot2)
library(reshape2)
library(coda)  # Needed for HPDinterval() used to summarize posterior distributions
library(tidyr)
library(reshape2)
library(stringr)

# ------------------------------------------------------------
# (Optional) This commented-out block shows how the final fitted models
# could be saved into a single .RData object for later reuse.
# In this workflow, the models are instead loaded from finalresults.RData.
# ------------------------------------------------------------

# save(OS_multinomial,
#      PP_multinomial,
#      OS_multinomial_weight,
#      PP_multinomial_weight,
#      OS_multinomial_growth,
#      PP_multinomial_growth,
#      OS_multinomial_bwkgxls,
#      PP_multinomial_bwkgxls,
#      OS_multinomial_activity,
#      PP_multinomial_activity,
#      OS_multinomial_aggressiveness,
#      PP_multinomial_aggressiveness,
#      OS_multinomial_trainability,
#      PP_multinomial_trainability,
#      
#      file = "finalresults.RData")

# Load previously fitted models/results (the full set of MCMCglmm objects)
load("finalresults.RData")

# ------------------------------------------------------------
# MCMCglmm is the Bayesian mixed-model framework used for the analyses.
# NOTE: Package version 2.29 is required to include both random structures
# (common ancestry + haplotype sharing) simultaneously (as in Garamszegi et al. 2020).
# ------------------------------------------------------------
library(MCMCglmm)

# ------------------------------------------------------------
# Input data:
# 1) Breed-level life-history / predictor database
# 2) Genetic similarity matrices from Garamszegi et al. 2020:
#    - SNP-based similarity (interpreted as common ancestry)
#    - Haplotype-based similarity (interpreted as haplotype sharing / gene flow)
# ------------------------------------------------------------

# Data: life history database
data<-read.csv("data.csv", header=TRUE)
# Data: genetic similarity matrices based on Garamszegi et al. 2020
snps<-read.csv("snps.csv", header=TRUE)
haplotypes<-read.csv("haplotypes.csv", header=TRUE)

# ------------------------------------------------------------
# Priors for Bayesian models (MCMCglmm)
# Priors are based on Hadfield (2010), with variants used depending on:
# - whether the response is univariate vs multinomial
# - number of response categories (OS vs PP)
# ------------------------------------------------------------

# Based on Hadfield (2010): basic structure for R (residual) and two G matrices
prior1<-list(
  R=list(V=1,nu=0.002),
  G=list(
    G1=list(V=1,nu=1),
    G2=list(V=1,nu=1)))

# Inverse-Wishart prior for both random effects (more diffuse settings for G components)
prior2<-list(
  R=list(V=1,nu=0.002),
  G=list(
    G1=list(V=1,nu=0.002),
    G2=list(V=1,nu=0.002)))

# For multinomial models with 8 variables (OS)
# (Here the dimension corresponds to the number of non-reference categories used by MCMCglmm)
prior3 <- list(
  R = list(V = diag(7), nu = 0.002),
  G = list(
    G1 = list(V = diag(7), nu = 0.002),
    G2 = list(V = diag(7), nu = 0.002)))

# For multinomial models with 7 variables (PP)
prior4 <- list(
  R = list(V = diag(6), nu = 0.002),
  G = list(
    G1 = list(V = diag(6), nu = 0.002),
    G2 = list(V = diag(6), nu = 0.002)))

# For multinomial models with 11 variables (OS)
prior5 <- list(
  R = list(V = diag(8), nu = 0.002),
  G = list(
    G1 = list(V = diag(8), nu = 0.002),
    G2 = list(V = diag(8), nu = 0.002)))

# ------------------------------------------------------------
# Diagnostic helper: autocorrelation plots for MCMC chains
# Input x is typically an mcmc object/matrix with posterior samples (columns = parameters).
# This function loops across parameters and plots ACF up to lag 100.
# ------------------------------------------------------------
plot.acfs <- function(x) {
  n <- dim(x)[2]
  par(mfrow=c(1,1))
  for (i in 1:n) {
    acf(x[,i], 
        lag.max=100, 
        main=colnames(x)[i])
    grid()}}

# ------------------------------------------------------------
# Exploratory data visualization:
# Histograms are used to inspect raw, scaled, and square-root-transformed distributions
# for key predictors (behaviour, life-history traits, physiological markers).
# This supports decisions about transformations and standardization.
# ------------------------------------------------------------

hist(data$activityCareau2010)
hist(scale(data$activityCareau2010))
hist(sqrt(data$activityCareau2010))

hist(data$aggressivenessCareau2010)
hist(scale(data$aggressivenessCareau2010))
hist(sqrt(data$aggressivenessCareau2010))

hist(data$trainabilityCareau2010)
hist(scale(data$trainabilityCareau2010))
hist(sqrt(data$trainabilityCareau2010))

hist(data$weight)
hist(scale(data$weight))
hist(sqrt(data$weight))

hist(data$lifespan)
hist(scale(data$lifespan))
hist(sqrt(data$lifespan))

hist(data$littersize)
hist(scale(data$littersize))
hist(sqrt(data$littersize))

hist(data$bw)
hist(scale(data$bw))
hist(sqrt(data$bw))

hist(data$bwkgxls)
hist(scale(data$bwkgxls))
hist(sqrt(data$bwkgxls))

hist(data$growth)
hist(scale(data$growth))
hist(sqrt(data$growth))

# ------------------------------------------------------------
# Standardization step ("Data transformation just in case"):
# Each key predictor is scaled (mean=0, SD=1) and stored as a new column in `data`.
# This creates standardized versions used in modelling to improve comparability
# and MCMC mixing across predictors with different scales.
# ------------------------------------------------------------

a23<-scale(data$activityCareau2010)
a24<-scale(data$aggressivenessCareau2010)
a25<-scale(data$trainabilityCareau2010)
a26<-scale(data$weight)
a27<-scale(data$lifespan)
a29<-scale(data$littersize)
a30<-scale(data$bwkgxls)
a31<-scale(data$growth)
a38<-scale(data$bw)

# Store scaled predictors in the main dataset with explicit column names
data["scaleactivityCareau2010"]<-a23
data["scaleaggressivenessCareau2010"]<-a24
data["scaletrainabilityCareau2010"]<-a25
data["scaleweight"]<-a26
data["scalelifespan"]<-a27
data["scalelittersize"]<-a29
data["scalebwkgxls"]<-a30
data["scalegrowth"]<-a31
data["scalebw"]<-a38

# ------------------------------------------------------------
# Subset construction for mortality outcomes (Fleming et al. 2011):
# Select only columns required for causes of death by:
# - Organ System (OS)
# - Pathophysiological Process (PP)
# plus breed identifier and sample size (n).
#
# The subset corresponds to the 72 breeds used in the core analyses.
# ------------------------------------------------------------

fleming<-subset(data, select=c(Garamszegi2020,
                               n,
                               OScardiovascular,
                               OSdermatologic,
                               OSendocrine,
                               OSgastrointestinal,
                               OShematopoietic,
                               OShepatic,
                               OSmusculoskeletal,
                               OSneurologic,
                               OSophthalmologic,
                               OSrespiratory,
                               OSurogenital,
                               OSunclear,
                               PPcongenital,
                               PPdegenenerative,
                               PPinfectious,
                               PPinflammatory,
                               PPmetabolic,
                               PPneoplastic,
                               PPtoxic,
                               PPtraumatic,
                               PPvascular,
                               PPunclear
))

# Remove any breeds with missing information in the selected fields
fleming<- na.omit(fleming)

# Ensure breed IDs are treated as character strings and used as row names
fleming$Garamszegi2020<-as.character(fleming$Garamszegi2020)
row.names(fleming)<-fleming$Garamszegi2020

# ------------------------------------------------------------
# Convert Fleming et al. 2011 proportional mortality by OS into counts:
# For each breed, multiply proportion by total n and round to the nearest integer.
# These *_N variables represent breed-level death counts by category.
# ------------------------------------------------------------

# Get the number of dogs per breed for every disease (OS)
fleming$OScv_N<-round(fleming$OScardiovascular*fleming$n)
fleming$OSderm_N<-round(fleming$OSdermatologic*fleming$n)
fleming$OSendo_N<-round(fleming$OSendocrine*fleming$n)
fleming$OSgi_N<-round(fleming$OSgastrointestinal*fleming$n)
fleming$OShem_N<-round(fleming$OShematopoietic*fleming$n)
fleming$OShep_N<-round(fleming$OShepatic*fleming$n)
fleming$OSms_N<-round(fleming$OSmusculoskeletal*fleming$n)
fleming$OSneuro_N<-round(fleming$OSneurologic*fleming$n)
fleming$OSophth_N<-round(fleming$OSophthalmologic*fleming$n)
fleming$OSresp_N<-round(fleming$OSrespiratory*fleming$n)
fleming$OSuro_N<-round(fleming$OSurogenital*fleming$n)
fleming$OSuncl_N<-round(fleming$OSunclear*fleming$n)

# ------------------------------------------------------------
# Convert Fleming et al. 2011 proportional mortality by PP into counts (rounded).
# ------------------------------------------------------------

# Get the number of dogs per breed for every disease (PP)
fleming$PPcongen_N<-round(fleming$PPcongenital*fleming$n)
fleming$PPdegen_N<-round(fleming$PPdegenenerative*fleming$n)
fleming$PPinfect_N<-round(fleming$PPinfectious*fleming$n)
fleming$PPinflam_N<-round(fleming$PPinflammatory*fleming$n)
fleming$PPmetab_N<-round(fleming$PPmetabolic*fleming$n)
fleming$PPneopl_N<-round(fleming$PPneoplastic*fleming$n)
fleming$PPtoxic_N<-round(fleming$PPtoxic*fleming$n)
fleming$PPtraum_N<-round(fleming$PPtraumatic*fleming$n)
fleming$PPvasc_N<-round(fleming$PPvascular*fleming$n)
fleming$PPuncl_N<-round(fleming$PPunclear*fleming$n)

# ------------------------------------------------------------
# Breed matching between mortality data (Fleming subset) and genetic matrices:
# Fleming et al. (2011) and Garamszegi et al. (2020) do not include exactly the same breeds.
# Therefore, the intersection is taken to keep only breeds with both:
#  - mortality outcome data, and
#  - available genetic similarity information.
# ------------------------------------------------------------

# common ancestry data matrix
snps_fleming_intersect <- intersect(fleming$Garamszegi2020, snps$Garamszegi2020)
snps_fleming_common <- snps[snps$Garamszegi2020 %in% snps_fleming_intersect, ]
fleming_filtered <- fleming[fleming$Garamszegi2020 %in% snps_fleming_intersect, ]

# Optional diagnostic checks (commented out): verifies breed order identity after filtering
# same <- identical(snps_fleming_common$Garamszegi2020, fleming_filtered$Garamszegi2020)
# if (same) {
#   print("same")
# } else {
#   print("not the same")
# }

# The SNP similarity matrix is transposed so that breeds become both row and column names.
# This step converts the stored format into a square matrix-like data frame indexed by breed.
snps_fleming_common_transposed <- as.data.frame(t(snps_fleming_common))
colnames(snps_fleming_common_transposed) <- snps_fleming_common_transposed[1,]
snps_fleming_common_transposed <- snps_fleming_common_transposed[-1,]

# Extract the breed ID column and retain only the intersecting set again (for safety)
snps2 <- rownames_to_column(snps_fleming_common_transposed, var = "Garamszegi2020")
snps_fleming_intersect <- intersect(fleming$Garamszegi2020, snps2$Garamszegi2020)
snps_fleming_common <- snps2[snps2$Garamszegi2020 %in% snps_fleming_intersect, ]
row.names(snps_fleming_common)<-snps_fleming_common$Garamszegi2020

# `ancestry` is the numeric SNP-based similarity matrix used as a random-effect structure proxy
ancestry <- subset(snps_fleming_common, select = -Garamszegi2020)

# gene flow data matrix
haplotypes_fleming_intersect <- intersect(fleming$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_fleming_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_fleming_intersect, ]
fleming_filtered <- fleming[fleming$Garamszegi2020 %in% haplotypes_fleming_intersect, ]

# Optional diagnostic checks (commented out): verifies breed order identity after filtering
# same <- identical(haplotypes_fleming_common$Garamszegi2020, fleming_filtered$Garamszegi2020)
# if (same) {
#   print("same")
# } else {
#   print("not the same")
# }

# Transpose haplotype sharing matrix so breeds become both row and column names
haplotypes_fleming_common_transposed <- as.data.frame(t(haplotypes_fleming_common))
colnames(haplotypes_fleming_common_transposed) <- haplotypes_fleming_common_transposed[1,]
haplotypes_fleming_common_transposed <- haplotypes_fleming_common_transposed[-1,]

haplotypes2 <- rownames_to_column(haplotypes_fleming_common_transposed, var = "Garamszegi2020")
haplotypes_fleming_intersect <- intersect(fleming$Garamszegi2020, haplotypes2$Garamszegi2020)
haplotypes_fleming_common <- haplotypes2[haplotypes2$Garamszegi2020 %in% haplotypes_fleming_intersect, ]
row.names(haplotypes_fleming_common)<-haplotypes_fleming_common$Garamszegi2020

# `geneflow` is the numeric haplotype-based similarity matrix (haplotype sharing)
geneflow <- subset(haplotypes_fleming_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# SVD transformation to include genetic matrices in MCMCglmm:
# MCMCglmm can incorporate these structures via a set of latent variables
# obtained through singular value decomposition (SVD), following Garamszegi et al. (2020).
# These are stored in `fleming$phylo` and `fleming$haplo`.
# ------------------------------------------------------------

# To include genetic data as a random factor in the MCMCglmm model:
ancestry<-data.frame(sapply(ancestry, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(ancestry))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
fleming$phylo = Snpsvd1

geneflow<-data.frame(sapply(geneflow, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(geneflow))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
fleming$haplo = Hapsvd2

# Re-assert breed IDs as row names and store the SVD-derived random-effect covariates
row.names(fleming)<-fleming$Garamszegi2020
fleming$phylo = Snpsvd1
fleming$haplo = Hapsvd2

###############################################################################
# Analysis: influence of common ancestry and haplotype sharing (gene flow proxy)
# on among-breed variation in causes of death by Organ System (OS)
#
# MODEL TYPE:
# - Bayesian mixed-effects MULTINOMIAL model (MCMCglmm)
# - NO fixed predictors (intercept-free trait parameterization only)
#
# BIOLOGICAL/STATISTICAL PURPOSE:
# - Partition among-breed variance in OS mortality categories into components
#   attributable to:
#   (i) common ancestry (SNP-based similarity; random effect: phylo)
#   (ii) haplotype sharing / gene flow (haplotype-based similarity; random effect: haplo)
#   (iii) residual ("units") variance at the breed level
#
# NOTE:
# - The response is a vector of breed-level counts across OS categories.
# - "trait - 1" estimates a separate baseline (intercept) for each OS category
#   (no overall reference category intercept in the fixed effects).
###############################################################################

# This line prints a subset of column names in `fleming` (columns 26, 27, 30, 33, 36).
# Practical purpose: to verify which OS count columns are being pooled into "others".
names(fleming)[c(26, 27, 30, 33, 36)]

# Create an aggregated OS category ("OSothers_N") by summing selected OS count columns.
# Purpose: reduce sparse/low-frequency categories by pooling them into a single "others"
# category, so the multinomial response has 9 categories in total (multinomial9).
fleming$OSothers_N<-rowSums(fleming[, c(26, 30, 33, 36)])

# Fit the multinomial MCMCglmm model:
OS_multinomial <- MCMCglmm(cbind(OScv_N,
                                 OSendo_N,
                                 OSgi_N, 
                                 OShem_N, 
                                 OSms_N, 
                                 OSneuro_N,
                                 OSresp_N, 
                                 OSuro_N,
                                 OSothers_N) ~ trait-1, 
                           random = ~idh(trait):phylo + idh(trait):haplo, 
                           rcov = ~idh(trait):units, 
                           data = fleming, 
                           family = "multinomial9",
                           prior = prior5,
                           thin = 400,
                           burnin = 100000,
                           nitt = 1800000,
                           verbose = TRUE)

# Basic model inspection:
# - plot(): trace/density style plots for quick visual diagnostics
# - summary(): posterior means, credible intervals, and variance components summaries
plot(OS_multinomial)
summary(OS_multinomial)

# Autocorrelation diagnostics for:
# - VCV: posterior samples of (co)variance components (random effects + residual)
# - Sol: posterior samples of fixed effects (trait-specific intercepts here)
autocorr.diag(OS_multinomial$VCV)
autocorr.diag(OS_multinomial$Sol)

# Custom ACF plots for each variance component (lag up to 100; defined previously as plot.acfs)
plot.acfs(OS_multinomial$VCV)

# ------------------------------------------------------------
# Posterior summarization: variance partitioning by category/component
#
# Goal:
# For each column in OS_multinomial$VCV, compute the posterior distribution of the
# proportion of total variance explained:
#   VCV[, i] / rowSums(VCV)
# Then summarize that posterior proportion using the 95% HPD interval.
# ------------------------------------------------------------

# Calculate the HPD intervals of variance proportions and store them in a list.
# Each list element corresponds to one variance component column in OS_multinomial$VCV.
OS_multinomial_HPD_intervals <- lapply(1:ncol(OS_multinomial$VCV), function(i) {
  HPDinterval(OS_multinomial$VCV[, i] / rowSums(OS_multinomial$VCV))
})

# Label each list element using the original VCV column names (e.g., trait:phylo, trait:haplo, trait:units by category)
names(OS_multinomial_HPD_intervals) <- colnames(OS_multinomial$VCV)

# Convert the list of HPD intervals into a single data.frame for easier inspection/export.
# Output columns:
# - Variable: original variance component label
# - Lower/Upper: 95% HPD bounds for the variance proportion
OS_multinomial_HPD_intervals_df <- do.call(rbind, lapply(names(OS_multinomial_HPD_intervals), function(name) {
  data.frame(
    Variable = name,
    Lower = OS_multinomial_HPD_intervals[[name]][, "lower"],
    Upper = OS_multinomial_HPD_intervals[[name]][, "upper"]
  )
}))

# View the HPD interval table (variance proportions) in the console
print(OS_multinomial_HPD_intervals_df)

# ------------------------------------------------------------
# Extended summary table: HPD + median for each variance proportion
#
# Goal:
# For each variance component column i:
# - compute posterior proportions (VCV[,i]/total)
# - compute HPD (95%) and median of that posterior distribution
# - store results as one row per component
# ------------------------------------------------------------

# Create a list to store row-wise summary results
OS_multinomial_results_list <- lapply(1:ncol(OS_multinomial$VCV), function(i) {
  # Posterior distribution of the proportion explained by component i
  OS_multinomial_variable_proportion <- OS_multinomial$VCV[, i] / rowSums(OS_multinomial$VCV)
  
  # Summaries of the posterior distribution:
  # - HPDinterval(): 95% highest posterior density interval
  # - median(): robust central tendency of the posterior distribution
  OS_multinomial_hpd <- HPDinterval(OS_multinomial_variable_proportion)
  OS_multinomial_median_val <- median(OS_multinomial_variable_proportion)
  
  # Store one row with component name + HPD bounds + median
  data.frame(
    Variable = colnames(OS_multinomial$VCV)[i],
    Lower = OS_multinomial_hpd[1, "lower"],
    Upper = OS_multinomial_hpd[1, "upper"],
    Median = OS_multinomial_median_val
  )
})

# Combine all rows into a single data.frame (one row per variance component)
OS_multinomial_results <- do.call(rbind, OS_multinomial_results_list)

# View the resulting summary table
print(OS_multinomial_results)

# ------------------------------------------------------------
# Human-readable labeling
#
# The VCV column names produced by MCMCglmm are typically technical and encode:
# - OS category
# - random effect source (phylo vs haplo vs units)
#
# Here, those labels are replaced by a manually curated vector of readable names.
# NOTE: This improves interpretability for tables/figures and manuscript reporting.
# ------------------------------------------------------------

# Correct variable names list
correct_variable_names <- c(
  "Cardiovascular (common ancestry)", "Endocrine (common ancestry)", "Gastrointestinal (common ancestry)", "Hematopoietic (common ancestry)", "Musculoskeletal (common ancestry)", "Neurologic (common ancestry)",
  "Respiratory (common ancestry)", "Urogenital (common ancestry)", "Cardiovascular (hybridization)", "Endocrine (hybridization)", "Gastrointestinal (hybridization)", "Hematopoietic (hybridization)",
  "Musculoskeletal (hybridization)", "Neurologic (hybridization)", "Respiratory (hybridization)", "Urogenital (hybridization)",
  "Cardiovascular (units)", "Endocrine (units)", "Gastrointestinal (units)", "Hematopoietic (units)", "Musculoskeletal (units)",
  "Neurologic (units)", "Respiratory (units)", "Urogenital (units)"
)

# Replace the technical "Variable" column with the human-readable labels above.
OS_multinomial_results$Variable <- correct_variable_names

# View the renamed table (now ready for export or manuscript table assembly)
print(OS_multinomial_results)

# Optional export (commented out): write the final summary table to CSV for downstream use
# write.csv(OS_multinomial_results, "OS_multinomial_results_final.csv", row.names = FALSE)

# ------------------------------------------------------------
# Convenience filter: remove "units" rows
#
# Rationale:
# For some reporting contexts (e.g., a table focused only on genetic components),
# it is useful to exclude the residual ("units") component and keep only:
# - common ancestry (phylo)
# - haplotype sharing / hybridization (haplo)
# ------------------------------------------------------------
OS_multinomial_results_filtered <- OS_multinomial_results[!grepl("units", OS_multinomial_results$Variable), ]
print(OS_multinomial_results_filtered)

###############################################################################
# Analysis: influence of common ancestry and haplotype sharing (gene flow proxy)
# on among-breed variation in causes of death by Organ System (OS)
#
# MODEL TYPE:
# - Bayesian mixed-effects MULTINOMIAL model (MCMCglmm)
# - NO fixed predictors (intercept-free trait parameterization only)
#
# BIOLOGICAL/STATISTICAL PURPOSE:
# - Partition among-breed variance in OS mortality categories into components
#   attributable to:
#   (i) common ancestry (SNP-based similarity; random effect: phylo)
#   (ii) haplotype sharing / gene flow (haplotype-based similarity; random effect: haplo)
#   (iii) residual ("units") variance at the breed level
#
# NOTE:
# - The response is a vector of breed-level counts across OS categories.
# - "trait - 1" estimates a separate baseline (intercept) for each OS category
#   (no overall reference category intercept in the fixed effects).
###############################################################################

# This line prints a subset of column names in `fleming` (columns 26, 27, 30, 33, 36).
# Practical purpose: to verify which OS count columns are being pooled into "others".
names(fleming)[c(26, 27, 30, 33, 36)]

# Create an aggregated OS category ("OSothers_N") by summing selected OS count columns.
# Purpose: reduce sparse/low-frequency categories by pooling them into a single "others"
# category, so the multinomial response has 9 categories in total (multinomial9).
fleming$OSothers_N<-rowSums(fleming[, c(26, 30, 33, 36)])

# Fit the multinomial MCMCglmm model:
# Response: cbind(...) of OS-specific death counts per breed (including pooled "others")
# Fixed effects: ~ trait - 1  (separate category-specific intercepts)
# Random effects:
#   - idh(trait):phylo  : category-specific random effects structured by common ancestry
#   - idh(trait):haplo  : category-specific random effects structured by haplotype sharing
# Residual structure:
#   - idh(trait):units  : category-specific residual variance (breed-level)
#
# family = "multinomial9" indicates 9 response categories.
# prior = prior5 corresponds to an OS multinomial structure with the appropriate dimension.
# thin/burnin/nitt define MCMC sampling and chain length for posterior inference.
OS_multinomial <- MCMCglmm(cbind(OScv_N,
                                 OSendo_N,
                                 OSgi_N, 
                                 OShem_N, 
                                 OSms_N, 
                                 OSneuro_N,
                                 OSresp_N, 
                                 OSuro_N,
                                 OSothers_N) ~ trait-1, 
                           random = ~idh(trait):phylo + idh(trait):haplo, 
                           rcov = ~idh(trait):units, 
                           data = fleming, 
                           family = "multinomial9",
                           prior = prior5,
                           thin = 400,
                           burnin = 100000,
                           nitt = 1800000,
                           verbose = TRUE)

# Basic model inspection:
# - plot(): trace/density style plots for quick visual diagnostics
# - summary(): posterior means, credible intervals, and variance components summaries
plot(OS_multinomial)
summary(OS_multinomial)

# Autocorrelation diagnostics for:
# - VCV: posterior samples of (co)variance components (random effects + residual)
# - Sol: posterior samples of fixed effects (trait-specific intercepts here)
autocorr.diag(OS_multinomial$VCV)
autocorr.diag(OS_multinomial$Sol)

# Custom ACF plots for each variance component (lag up to 100; defined previously as plot.acfs)
plot.acfs(OS_multinomial$VCV)

# ------------------------------------------------------------
# Posterior summarization: variance partitioning by category/component
#
# Goal:
# For each column in OS_multinomial$VCV, compute the posterior distribution of the
# proportion of total variance explained:
#   VCV[, i] / rowSums(VCV)
# Then summarize that posterior proportion using the 95% HPD interval.
# ------------------------------------------------------------

# Calculate the HPD intervals of variance proportions and store them in a list.
# Each list element corresponds to one variance component column in OS_multinomial$VCV.
OS_multinomial_HPD_intervals <- lapply(1:ncol(OS_multinomial$VCV), function(i) {
  HPDinterval(OS_multinomial$VCV[, i] / rowSums(OS_multinomial$VCV))
})

# Label each list element using the original VCV column names (e.g., trait:phylo, trait:haplo, trait:units by category)
names(OS_multinomial_HPD_intervals) <- colnames(OS_multinomial$VCV)

# Convert the list of HPD intervals into a single data.frame for easier inspection/export.
# Output columns:
# - Variable: original variance component label
# - Lower/Upper: 95% HPD bounds for the variance proportion
OS_multinomial_HPD_intervals_df <- do.call(rbind, lapply(names(OS_multinomial_HPD_intervals), function(name) {
  data.frame(
    Variable = name,
    Lower = OS_multinomial_HPD_intervals[[name]][, "lower"],
    Upper = OS_multinomial_HPD_intervals[[name]][, "upper"]
  )
}))

# View the HPD interval table (variance proportions) in the console
print(OS_multinomial_HPD_intervals_df)

# ------------------------------------------------------------
# Extended summary table: HPD + median for each variance proportion
#
# Goal:
# For each variance component column i:
# - compute posterior proportions (VCV[,i]/total)
# - compute HPD (95%) and median of that posterior distribution
# - store results as one row per component
# ------------------------------------------------------------

# Create a list to store row-wise summary results
OS_multinomial_results_list <- lapply(1:ncol(OS_multinomial$VCV), function(i) {
  # Posterior distribution of the proportion explained by component i
  OS_multinomial_variable_proportion <- OS_multinomial$VCV[, i] / rowSums(OS_multinomial$VCV)
  
  # Summaries of the posterior distribution:
  # - HPDinterval(): 95% highest posterior density interval
  # - median(): robust central tendency of the posterior distribution
  OS_multinomial_hpd <- HPDinterval(OS_multinomial_variable_proportion)
  OS_multinomial_median_val <- median(OS_multinomial_variable_proportion)
  
  # Store one row with component name + HPD bounds + median
  data.frame(
    Variable = colnames(OS_multinomial$VCV)[i],
    Lower = OS_multinomial_hpd[1, "lower"],
    Upper = OS_multinomial_hpd[1, "upper"],
    Median = OS_multinomial_median_val
  )
})

# Combine all rows into a single data.frame (one row per variance component)
OS_multinomial_results <- do.call(rbind, OS_multinomial_results_list)

# View the resulting summary table
print(OS_multinomial_results)

# ------------------------------------------------------------
# Human-readable labeling
#
# The VCV column names produced by MCMCglmm are typically technical and encode:
# - OS category
# - random effect source (phylo vs haplo vs units)
#
# Here, those labels are replaced by a manually curated vector of readable names.
# NOTE: This improves interpretability for tables/figures and manuscript reporting.
# ------------------------------------------------------------

# Correct variable names list
correct_variable_names <- c(
  "Cardiovascular (common ancestry)", "Endocrine (common ancestry)", "Gastrointestinal (common ancestry)", "Hematopoietic (common ancestry)", "Musculoskeletal (common ancestry)", "Neurologic (common ancestry)",
  "Respiratory (common ancestry)", "Urogenital (common ancestry)", "Cardiovascular (hybridization)", "Endocrine (hybridization)", "Gastrointestinal (hybridization)", "Hematopoietic (hybridization)",
  "Musculoskeletal (hybridization)", "Neurologic (hybridization)", "Respiratory (hybridization)", "Urogenital (hybridization)",
  "Cardiovascular (units)", "Endocrine (units)", "Gastrointestinal (units)", "Hematopoietic (units)", "Musculoskeletal (units)",
  "Neurologic (units)", "Respiratory (units)", "Urogenital (units)"
)

# Replace the technical "Variable" column with the human-readable labels above.
OS_multinomial_results$Variable <- correct_variable_names

# View the renamed table (now ready for export or manuscript table assembly)
print(OS_multinomial_results)

# Optional export (commented out): write the final summary table to CSV for downstream use
# write.csv(OS_multinomial_results, "OS_multinomial_results_final.csv", row.names = FALSE)

# ------------------------------------------------------------
# Convenience filter: remove "units" rows
#
# Rationale:
# For some reporting contexts (e.g., a table focused only on genetic components),
# it is useful to exclude the residual ("units") component and keep only:
# - common ancestry (phylo)
# - haplotype sharing / hybridization (haplo)
# ------------------------------------------------------------
OS_multinomial_results_filtered <- OS_multinomial_results[!grepl("units", OS_multinomial_results$Variable), ]
print(OS_multinomial_results_filtered)

###############################################################################
###############################################################################
# ------------------------------------------------------------
# Subset construction: breeds with causes of death + body weight data
# (71 breeds after removing missing values)
#
# PURPOSE:
# - Create a dataset restricted to breeds with:
#   (i) mortality outcomes (OS and PP proportions + sample size n), and
#   (ii) standardized body weight (scaleweight)
# - Convert proportional mortality values into breed-level death counts
# - Match breeds with available genetic matrices (SNPs and haplotypes)
# - Create SVD-based random-effect covariates (phylo and haplo) for MCMCglmm
# ------------------------------------------------------------

# Subset: only breeds with causes of death + weight data (71 breeds)  
dfweight<-subset(data, select=c(Garamszegi2020,
                                n,
                                scaleweight,
                                OScardiovascular,
                                OSdermatologic,
                                OSendocrine,
                                OSgastrointestinal,
                                OShematopoietic,
                                OShepatic,
                                OSmusculoskeletal,
                                OSneurologic,
                                OSophthalmologic,
                                OSrespiratory,
                                OSurogenital,
                                OSunclear,
                                PPcongenital,
                                PPdegenenerative,
                                PPinfectious,
                                PPinflammatory,
                                PPmetabolic,
                                PPneoplastic,
                                PPtoxic,
                                PPtraumatic,
                                PPvascular,
                                PPunclear
))

# Remove rows with missing values in any of the selected columns
# (ensures the subset used in models has complete information)
dfweight<- na.omit(dfweight)

# Ensure breed IDs are treated as character strings and used as row names
dfweight$Garamszegi2020<-as.character(dfweight$Garamszegi2020)
row.names(dfweight)<-dfweight$Garamszegi2020 

# ------------------------------------------------------------
# Convert OS proportional mortality into breed-level death counts:
# For each breed, count = round(proportion * n).
# These *_N columns are used as the multinomial response (counts).
# ------------------------------------------------------------

# Get the number of dogs per breed for every disease (OS)
dfweight$OScv_N<-round(dfweight$OScardiovascular*dfweight$n)
dfweight$OSderm_N<-round(dfweight$OSdermatologic*dfweight$n)
dfweight$OSendo_N<-round(dfweight$OSendocrine*dfweight$n)
dfweight$OSgi_N<-round(dfweight$OSgastrointestinal*dfweight$n)
dfweight$OShem_N<-round(dfweight$OShematopoietic*dfweight$n)
dfweight$OShep_N<-round(dfweight$OShepatic*dfweight$n)
dfweight$OSms_N<-round(dfweight$OSmusculoskeletal*dfweight$n)
dfweight$OSneuro_N<-round(dfweight$OSneurologic*dfweight$n)
dfweight$OSophth_N<-round(dfweight$OSophthalmologic*dfweight$n)
dfweight$OSresp_N<-round(dfweight$OSrespiratory*dfweight$n)
dfweight$OSuro_N<-round(dfweight$OSurogenital*dfweight$n)
dfweight$OSuncl_N<-round(dfweight$OSunclear*dfweight$n)

# ------------------------------------------------------------
# Convert PP proportional mortality into breed-level death counts:
# For each breed, count = round(proportion * n).
# These *_N columns are used as the multinomial response (counts).
# ------------------------------------------------------------

# Get the number of dogs per breed for every disease (PP)
dfweight$PPcongen_N<-round(dfweight$PPcongenital*dfweight$n)
dfweight$PPdegen_N<-round(dfweight$PPdegenenerative*dfweight$n)
dfweight$PPinfect_N<-round(dfweight$PPinfectious*dfweight$n)
dfweight$PPinflam_N<-round(dfweight$PPinflammatory*dfweight$n)
dfweight$PPmetab_N<-round(dfweight$PPmetabolic*dfweight$n)
dfweight$PPneopl_N<-round(dfweight$PPneoplastic*dfweight$n)
dfweight$PPtoxic_N<-round(dfweight$PPtoxic*dfweight$n)
dfweight$PPtraum_N<-round(dfweight$PPtraumatic*dfweight$n)
dfweight$PPvasc_N<-round(dfweight$PPvascular*dfweight$n)
dfweight$PPuncl_N<-round(dfweight$PPunclear*dfweight$n)

# ------------------------------------------------------------
# Match breeds between dfweight and SNP similarity matrix (common ancestry):
# - Keep only breeds with both weight+mortality data AND SNP matrix entries
# - Transpose and reformat the matrix so breeds are indexed consistently
# ------------------------------------------------------------

# common ancestry data matrix
snps_dfweight_intersect <- intersect(dfweight$Garamszegi2020, snps$Garamszegi2020)
snps_dfweight_common <- snps[snps$Garamszegi2020 %in% snps_dfweight_intersect, ]
dfweight_filtered <- dfweight[dfweight$Garamszegi2020 %in% snps_dfweight_intersect, ]

# Transpose to obtain a square matrix-like format with breed IDs as column names
snps_dfweight_common_transposed <- as.data.frame(t(snps_dfweight_common))
colnames(snps_dfweight_common_transposed) <- snps_dfweight_common_transposed[1,]
snps_dfweight_common_transposed <- snps_dfweight_common_transposed[-1,]

# Move row names to a "Garamszegi2020" column and re-filter to the intersecting breed set
snps2 <- rownames_to_column(snps_dfweight_common_transposed, var = "Garamszegi2020")
snps_dfweight_intersect <- intersect(dfweight$Garamszegi2020, snps2$Garamszegi2020)
snps_dfweight_common <- snps2[snps2$Garamszegi2020 %in% snps_dfweight_intersect, ]
row.names(snps_dfweight_common)<-snps_dfweight_common$Garamszegi2020

# Remove the ID column to obtain the numeric similarity matrix used downstream
ancestry_weight <- subset(snps_dfweight_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# Match breeds between dfweight and haplotype similarity matrix (haplotype sharing):
# - Keep only breeds with both weight+mortality data AND haplotype matrix entries
# - Transpose and reformat the matrix so breeds are indexed consistently
# ------------------------------------------------------------

# gene flow data matrix
haplotypes_dfweight_intersect <- intersect(dfweight$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_dfweight_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_dfweight_intersect, ]
dfweight_filtered <- dfweight[dfweight$Garamszegi2020 %in% haplotypes_dfweight_intersect, ]

# Transpose to obtain a square matrix-like format with breed IDs as column names
haplotypes_dfweight_common_transposed <- as.data.frame(t(haplotypes_dfweight_common))
colnames(haplotypes_dfweight_common_transposed) <- haplotypes_dfweight_common_transposed[1,]
haplotypes_dfweight_common_transposed <- haplotypes_dfweight_common_transposed[-1,]

# Move row names to a "Garamszegi2020" column and re-filter to the intersecting breed set
haplotypes2 <- rownames_to_column(haplotypes_dfweight_common_transposed, var = "Garamszegi2020")
haplotypes_dfweight_intersect <- intersect(dfweight$Garamszegi2020, haplotypes2$Garamszegi2020)
haplotypes_dfweight_common <- haplotypes2[haplotypes2$Garamszegi2020 %in% haplotypes_dfweight_intersect, ]
row.names(haplotypes_dfweight_common)<-haplotypes_dfweight_common$Garamszegi2020

# Remove the ID column to obtain the numeric similarity matrix used downstream
geneflow_weight <- subset(haplotypes_dfweight_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# SVD transformation to include genetic matrices as random factors in MCMCglmm:
# - Convert all entries to numeric
# - Apply SVD and compute the transformed matrix used as latent covariates
# - Store the resulting matrices in dfweight$phylo and dfweight$haplo
#   (following the same approach used for the no-predictor models)
# ------------------------------------------------------------

# To include genetic data as a random factor in the MCMCglmm model:
ancestry_weight<-data.frame(sapply(ancestry_weight, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(ancestry_weight))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfweight$phylo = Snpsvd1

geneflow_weight<-data.frame(sapply(geneflow_weight, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(geneflow_weight))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfweight$haplo = Hapsvd2

# Re-assert breed IDs as row names and store the SVD-derived random-effect covariates
row.names(dfweight)<-dfweight$Garamszegi2020
dfweight$phylo = Snpsvd1
dfweight$haplo = Hapsvd2

# ------------------------------------------------------------
# OS multinomial model including body weight as a predictor
#
# MODEL:
# - Response: OS category-specific death counts per breed (9 categories incl. "others")
# - Fixed effects:
#   ~ trait - 1 + trait:scaleweight
#   This estimates a separate weight slope for each OS category (no global slope).
# - Random effects:
#   idh(trait):phylo + idh(trait):haplo (genetic influence components by category)
# - Residual:
#   idh(trait):units
# ------------------------------------------------------------

# Analysis of the influence of common ancestry and gene flow on the among-breed variation in diseases + weight data (OS)
# WEIGHT DATA MULTINOMIAL MODEL (OS)

# Pool selected OS count columns into a single "others" category
dfweight$OSothers_N<-rowSums(dfweight[, c(27, 31, 34, 37)])

# Fit the multinomial model with category-specific weight effects
OS_multinomial_weight <- MCMCglmm(cbind(OScv_N,
                                        OSendo_N,
                                        OSgi_N, 
                                        OShem_N, 
                                        OSms_N, 
                                        OSneuro_N,
                                        OSresp_N, 
                                        OSuro_N,
                                        OSothers_N) ~ trait-1+trait:scaleweight,
                                  random = ~idh(trait):phylo + idh(trait):haplo, 
                                  rcov = ~idh(trait):units, 
                                  data = dfweight, 
                                  family = "multinomial9",
                                  prior = prior5,
                                  thin = 400,
                                  burnin = 100000,
                                  nitt = 1800000,
                                  verbose = TRUE)

# Optional quick visual inspection of traces/densities (commented out in this workflow)
# plot(OS_multinomial_weight)

# Summary of posterior estimates for fixed effects (including trait:scaleweight) and variance components
summary(OS_multinomial_weight)

# Autocorrelation diagnostics for posterior chains:
# - VCV: variance components (phylo, haplo, units) by OS category
# - Sol: fixed effects (trait-specific intercepts + trait:scaleweight slopes)
autocorr.diag(OS_multinomial_weight$VCV)
autocorr.diag(OS_multinomial_weight$Sol)

# Custom ACF plots for each variance component column (function defined earlier: plot.acfs)
plot.acfs(OS_multinomial_weight$VCV)

# ------------------------------------------------------------
# Extract category-specific weight effects from the fitted model
#
# Goal:
# - Retrieve the posterior summaries for "trait:scaleweight" coefficients
# - Create a tidy data frame with:
#   Trait (OS category), Estimate (posterior mean), and 95% CI bounds
# - This object is typically used to build coefficient plots or tables
# ------------------------------------------------------------

# Extract interaction coefficients (scaleweight * trait) from OS_multinomial_weight
summary_OS <- summary(OS_multinomial_weight)$solutions
coef_data_OS <- as.data.frame(summary_OS) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scaleweight", trait)) %>%
  mutate(Trait = gsub("trait|:scaleweight", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Set plotting order of cause-of-death categories (preserves current order in the data frame)
coef_data_OS$Trait <- factor(coef_data_OS$Trait, levels = unique(coef_data_OS$Trait))

# ------------------------------------------------------------
# Analysis of the influence of common ancestry and gene flow
# on among-breed variation in diseases + body weight data (PP)
#
# MODEL:
# - Bayesian mixed-effects MULTINOMIAL model (MCMCglmm)
# - Response: PP category-specific death counts per breed
# - Fixed effects:
#   ~ trait - 1 + trait:scaleweight
#   This specification estimates a separate body-weight slope
#   for each PP category (no global intercept or slope).
# - Random effects:
#   idh(trait):phylo  -> variance structured by common ancestry
#   idh(trait):haplo  -> variance structured by haplotype sharing (gene flow proxy)
# - Residual:
#   idh(trait):units
# ------------------------------------------------------------

# WEIGHT DATA MULTINOMIAL MODEL (PP)

# Pool selected PP count columns into a single "others" category.
# This reduces sparsity in low-frequency PP categories and yields
# a total of 8 multinomial response categories.
dfweight$PPothers_N<-rowSums(dfweight[, c(44, 46, 47)])

# Fit the multinomial MCMCglmm model with category-specific weight effects:
# - Response: PP death counts per breed (including pooled "others")
# - Fixed effects: trait-specific intercepts and trait-specific slopes for scaleweight
# - Random effects: genetic influence via common ancestry (phylo) and haplotype sharing (haplo)
# - family = "multinomial8" specifies an 8-category multinomial response
PP_multinomial_weight=MCMCglmm(cbind(PPcongen_N, 
                                     PPdegen_N,
                                     PPinflam_N,
                                     PPinfect_N,
                                     PPmetab_N, 
                                     PPneopl_N,
                                     PPtraum_N,
                                     PPothers_N) ~ trait-1+trait:scaleweight, 
                               random = ~idh(trait):phylo + idh(trait):haplo, 
                               rcov = ~idh(trait):units, 
                               data = dfweight, 
                               family = "multinomial8",
                               prior = prior3,
                               thin = 400,
                               burnin = 100000,
                               nitt = 1800000,
                               verbose = TRUE)

# Optional visual inspection of MCMC traces and posterior densities
# plot(PP_multinomial_weight)

# Summary of posterior estimates:
# - Fixed effects (trait-specific intercepts and weight slopes)
# - Variance components for phylo, haplo, and units
summary(PP_multinomial_weight)

# Autocorrelation diagnostics for posterior chains:
# - VCV: variance components by PP category
# - Sol: fixed effects (trait and trait:scaleweight terms)
autocorr.diag(PP_multinomial_weight$VCV)
autocorr.diag(PP_multinomial_weight$Sol)

# Custom autocorrelation plots for each variance component
# (function plot.acfs defined earlier)
plot.acfs(PP_multinomial_weight$VCV)

# ------------------------------------------------------------
# Extraction of category-specific weight effects
#
# Goal:
# - Isolate posterior summaries for the interaction terms
#   (trait:scaleweight)
# - Build a tidy data frame with posterior means and 95% CIs
#   for each PP category
# ------------------------------------------------------------

# Extract coefficients for the interaction between scaleweight and trait
summary_PP <- summary(PP_multinomial_weight)$solutions
coef_data_PP <- as.data.frame(summary_PP) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scaleweight", trait)) %>%
  mutate(Trait = gsub("trait|:scaleweight", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Define plotting order of PP categories
# (preserves the order in which categories appear in the model output)
coef_data_PP$Trait <- factor(coef_data_PP$Trait, levels = unique(coef_data_PP$Trait))

# ------------------------------------------------------------
# Subset: breeds with causes of death + (scaled) growth data
# (56 breeds after removing missing values)
#
# PURPOSE:
# - Create a dataset restricted to breeds with:
#   (i) mortality outcomes (OS and PP proportions + sample size n), and
#   (ii) standardized growth predictor (scalegrowth)
# - Convert proportional mortality values into breed-level death counts
# - Match breeds with available genetic matrices (SNPs and haplotypes)
# - Create SVD-based random-effect covariates (phylo and haplo) for MCMCglmm
# ------------------------------------------------------------

# Subset: only breeds with causes of death + (sqrt) growth data (56 breeds) 
dfgrowth<-subset(data, select=c(Garamszegi2020,
                                n,
                                scalegrowth,
                                OScardiovascular,
                                OSdermatologic,
                                OSendocrine,
                                OSgastrointestinal,
                                OShematopoietic,
                                OShepatic,
                                OSmusculoskeletal,
                                OSneurologic,
                                OSophthalmologic,
                                OSrespiratory,
                                OSurogenital,
                                OSunclear,
                                PPcongenital,
                                PPdegenenerative,
                                PPinfectious,
                                PPinflammatory,
                                PPmetabolic,
                                PPneoplastic,
                                PPtoxic,
                                PPtraumatic,
                                PPvascular,
                                PPunclear
))

# Remove rows with missing values in any selected column
dfgrowth<- na.omit(dfgrowth)

# Ensure breed IDs are treated as character strings and used as row names
dfgrowth$Garamszegi2020<-as.character(dfgrowth$Garamszegi2020)
row.names(dfgrowth)<-dfgrowth$Garamszegi2020 

# ------------------------------------------------------------
# Convert OS proportional mortality into breed-level death counts:
# For each breed, count = round(proportion * n).
# These *_N columns are used as the multinomial response (counts).
# ------------------------------------------------------------

# Get the number of dogs per breed for every disease (OS)
dfgrowth$OScv_N<-round(dfgrowth$OScardiovascular*dfgrowth$n)
dfgrowth$OSderm_N<-round(dfgrowth$OSdermatologic*dfgrowth$n)
dfgrowth$OSendo_N<-round(dfgrowth$OSendocrine*dfgrowth$n)
dfgrowth$OSgi_N<-round(dfgrowth$OSgastrointestinal*dfgrowth$n)
dfgrowth$OShem_N<-round(dfgrowth$OShematopoietic*dfgrowth$n)
dfgrowth$OShep_N<-round(dfgrowth$OShepatic*dfgrowth$n)
dfgrowth$OSms_N<-round(dfgrowth$OSmusculoskeletal*dfgrowth$n)
dfgrowth$OSneuro_N<-round(dfgrowth$OSneurologic*dfgrowth$n)
dfgrowth$OSophth_N<-round(dfgrowth$OSophthalmologic*dfgrowth$n)
dfgrowth$OSresp_N<-round(dfgrowth$OSrespiratory*dfgrowth$n)
dfgrowth$OSuro_N<-round(dfgrowth$OSurogenital*dfgrowth$n)
dfgrowth$OSuncl_N<-round(dfgrowth$OSunclear*dfgrowth$n)

# ------------------------------------------------------------
# Convert PP proportional mortality into breed-level death counts:
# For each breed, count = round(proportion * n).
# These *_N columns are used as the multinomial response (counts).
# ------------------------------------------------------------

# Get the number of dogs per breed for every disease (PP)
dfgrowth$PPcongen_N<-round(dfgrowth$PPcongenital*dfgrowth$n)
dfgrowth$PPdegen_N<-round(dfgrowth$PPdegenenerative*dfgrowth$n)
dfgrowth$PPinfect_N<-round(dfgrowth$PPinfectious*dfgrowth$n)
dfgrowth$PPinflam_N<-round(dfgrowth$PPinflammatory*dfgrowth$n)
dfgrowth$PPmetab_N<-round(dfgrowth$PPmetabolic*dfgrowth$n)
dfgrowth$PPneopl_N<-round(dfgrowth$PPneoplastic*dfgrowth$n)
dfgrowth$PPtoxic_N<-round(dfgrowth$PPtoxic*dfgrowth$n)
dfgrowth$PPtraum_N<-round(dfgrowth$PPtraumatic*dfgrowth$n)
dfgrowth$PPvasc_N<-round(dfgrowth$PPvascular*dfgrowth$n)
dfgrowth$PPuncl_N<-round(dfgrowth$PPunclear*dfgrowth$n)

# ------------------------------------------------------------
# Match breeds between dfgrowth and SNP similarity matrix (common ancestry):
# - Keep only breeds with both growth+mortality data AND SNP matrix entries
# - Transpose and reformat the matrix so breeds are indexed consistently
# ------------------------------------------------------------

# common ancestry data matrix
snps_dfgrowth_intersect <- intersect(dfgrowth$Garamszegi2020, snps$Garamszegi2020)
snps_dfgrowth_common <- snps[snps$Garamszegi2020 %in% snps_dfgrowth_intersect, ]
dfgrowth_filtered <- dfgrowth[dfgrowth$Garamszegi2020 %in% snps_dfgrowth_intersect, ]

# Transpose to obtain a square matrix-like format with breed IDs as column names
snps_dfgrowth_common_transposed <- as.data.frame(t(snps_dfgrowth_common))
colnames(snps_dfgrowth_common_transposed) <- snps_dfgrowth_common_transposed[1,]
snps_dfgrowth_common_transposed <- snps_dfgrowth_common_transposed[-1,]

# Move row names to a "Garamszegi2020" column and re-filter to the intersecting breed set
snps2 <- rownames_to_column(snps_dfgrowth_common_transposed, var = "Garamszegi2020")
snps_dfgrowth_intersect <- intersect(dfgrowth$Garamszegi2020, snps2$Garamszegi2020)
snps_dfgrowth_common <- snps2[snps2$Garamszegi2020 %in% snps_dfgrowth_intersect, ]
row.names(snps_dfgrowth_common)<-snps_dfgrowth_common$Garamszegi2020

# Remove the ID column to obtain the numeric similarity matrix used downstream
ancestry_growth <- subset(snps_dfgrowth_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# Match breeds between dfgrowth and haplotype similarity matrix (haplotype sharing):
# - Keep only breeds with both growth+mortality data AND haplotype matrix entries
# - Transpose and reformat the matrix so breeds are indexed consistently
# ------------------------------------------------------------

# gene flow data matrix
haplotypes_dfgrowth_intersect <- intersect(dfgrowth$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_dfgrowth_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_dfgrowth_intersect, ]
dfgrowth_filtered <- dfgrowth[dfgrowth$Garamszegi2020 %in% haplotypes_dfgrowth_intersect, ]

# Transpose to obtain a square matrix-like format with breed IDs as column names
haplotypes_dfgrowth_common_transposed <- as.data.frame(t(haplotypes_dfgrowth_common))
colnames(haplotypes_dfgrowth_common_transposed) <- haplotypes_dfgrowth_common_transposed[1,]
haplotypes_dfgrowth_common_transposed <- haplotypes_dfgrowth_common_transposed[-1,]

# Move row names to a "Garamszegi2020" column and re-filter to the intersecting breed set
haplotypes2 <- rownames_to_column(haplotypes_dfgrowth_common_transposed, var = "Garamszegi2020")
haplotypes_dfgrowth_intersect <- intersect(dfgrowth$Garamszegi2020, haplotypes2$Garamszegi2020)
haplotypes_dfgrowth_common <- haplotypes2[haplotypes2$Garamszegi2020 %in% haplotypes_dfgrowth_intersect, ]
row.names(haplotypes_dfgrowth_common)<-haplotypes_dfgrowth_common$Garamszegi2020

# Remove the ID column to obtain the numeric similarity matrix used downstream
geneflow_growth <- subset(haplotypes_dfgrowth_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# SVD transformation to include genetic matrices as random factors in MCMCglmm:
# - Convert all entries to numeric
# - Apply SVD and compute the transformed matrix used as latent covariates
# - Store the resulting matrices in dfgrowth$phylo and dfgrowth$haplo
# ------------------------------------------------------------

# To include genetic data as a random factor in the MCMCglmm model:
ancestry_growth<-data.frame(sapply(ancestry_growth, function(x) as.numeric(as.character(x))))
Snpsvd1 = svd(as.matrix(ancestry_growth))
Snpsvd1 = Snpsvd1$v%*%(t(Snpsvd1$u)*sqrt(Snpsvd1$d))
dfgrowth$phylo = Snpsvd1

geneflow_growth<-data.frame(sapply(geneflow_growth, function(x) as.numeric(as.character(x))))
Hapsvd2 = svd(as.matrix(geneflow_growth))
Hapsvd2 = Hapsvd2$v%*%(t(Hapsvd2$u)*sqrt(Hapsvd2$d))
dfgrowth$haplo = Hapsvd2

# Re-assert breed IDs as row names and store the SVD-derived random-effect covariates
row.names(dfgrowth)<-dfgrowth$Garamszegi2020
dfgrowth$phylo = Snpsvd1
dfgrowth$haplo = Hapsvd2

# ------------------------------------------------------------
# OS multinomial model including growth as a predictor
#
# MODEL:
# - Response: OS category-specific death counts per breed (9 categories incl. "others")
# - Fixed effects:
#   ~ trait - 1 + trait:scalegrowth
#   This estimates a separate growth slope for each OS category (no global slope).
# - Random effects:
#   idh(trait):phylo + idh(trait):haplo (genetic influence components by category)
# - Residual:
#   idh(trait):units
# ------------------------------------------------------------

# Analysis of the influence of common ancestry and gene flow on the among-breed variation in diseases + growth data (OS)
# GROWTH DATA MULTINOMIAL MODEL (OS)

# Pool selected OS count columns into a single "others" category
dfgrowth$OSothers_N <- rowSums(dfgrowth[, c(27, 31, 34, 37)])

# Fit the multinomial model with category-specific growth effects
OS_multinomial_growth <- MCMCglmm(cbind(OScv_N,
                                        OSendo_N,
                                        OSgi_N, 
                                        OShem_N, 
                                        OSms_N, 
                                        OSneuro_N,
                                        OSresp_N, 
                                        OSuro_N,
                                        OSothers_N) ~ trait-1+trait:scalegrowth,
                                  random = ~idh(trait):phylo + idh(trait):haplo, 
                                  rcov = ~idh(trait):units, 
                                  data = dfgrowth, 
                                  family = "multinomial9",
                                  prior = prior5,
                                  thin = 400,
                                  burnin = 100000,
                                  nitt = 1800000,
                                  verbose = TRUE)

# Optional visual inspection of MCMC traces and posterior densities
# plot(OS_multinomial_growth)

# Summary of posterior estimates for fixed effects (including trait:scalegrowth) and variance components
summary(OS_multinomial_growth)

# Autocorrelation diagnostics for posterior chains:
# - VCV: variance components (phylo, haplo, units) by OS category
# - Sol: fixed effects (trait-specific intercepts + trait:scalegrowth slopes)
autocorr.diag(OS_multinomial_growth$VCV)
autocorr.diag(OS_multinomial_growth$Sol)

# Custom autocorrelation plots for each variance component column
plot.acfs(OS_multinomial_growth$VCV)

# ------------------------------------------------------------
# Extract category-specific growth effects from the fitted model
#
# Goal:
# - Retrieve posterior summaries for "trait:scalegrowth" coefficients
# - Create a tidy data frame with:
#   Trait (OS category), Estimate (posterior mean), and 95% CI bounds
# - This object is typically used to build coefficient plots or tables
# ------------------------------------------------------------

# Extract interaction coefficients (scalegrowth * trait) from OS_multinomial_growth
summary_OS_growth <- summary(OS_multinomial_growth)$solutions
coef_data_OS_growth <- as.data.frame(summary_OS_growth) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scalegrowth", trait)) %>%
  mutate(Trait = gsub("trait|:scalegrowth", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Define plotting order of cause-of-death categories
# (preserves the order in which categories appear in the model output)
coef_data_OS_growth$Trait <- factor(coef_data_OS_growth$Trait, levels = unique(coef_data_OS_growth$Trait))

# ------------------------------------------------------------
# Analysis of the influence of common ancestry and gene flow
# on among-breed variation in diseases + growth data (PP)
#
# MODEL:
# - Bayesian mixed-effects MULTINOMIAL model (MCMCglmm)
# - Response: PP category-specific death counts per breed
# - Fixed effects:
#   ~ trait - 1 + trait:scalegrowth
#   This specification estimates a separate growth slope
#   for each PP category (no global intercept or slope).
# - Random effects:
#   idh(trait):phylo  -> variance structured by common ancestry
#   idh(trait):haplo  -> variance structured by haplotype sharing (gene flow proxy)
# - Residual:
#   idh(trait):units
# ------------------------------------------------------------

# GROWTH DATA MULTINOMIAL MODEL (PP)

# Pool selected PP count columns into a single "others" category.
# This reduces sparsity in low-frequency PP categories and yields
# a total of 8 multinomial response categories.
dfgrowth$PPothers_N <- rowSums(dfgrowth[, c(44, 46, 47)])

# Fit the multinomial MCMCglmm model with category-specific growth effects:
# - Response: PP death counts per breed (including pooled "others")
# - Fixed effects: trait-specific intercepts and trait-specific slopes for scalegrowth
# - Random effects: genetic influence via common ancestry (phylo) and haplotype sharing (haplo)
# - family = "multinomial8" specifies an 8-category multinomial response
PP_multinomial_growth = MCMCglmm(cbind(PPinfect_N, 
                                       PPcongen_N, 
                                       PPdegen_N, 
                                       PPinflam_N, 
                                       PPmetab_N, 
                                       PPneopl_N,
                                       PPtraum_N,
                                       PPothers_N) ~ trait-1+trait:scalegrowth, 
                                 random = ~idh(trait):phylo + idh(trait):haplo, 
                                 rcov = ~idh(trait):units, 
                                 data = dfgrowth, 
                                 family = "multinomial8",
                                 prior = prior3,
                                 thin = 400,
                                 burnin = 100000,
                                 nitt = 1800000,
                                 verbose = TRUE)

# Optional visual inspection of MCMC traces and posterior densities
# plot(PP_multinomial)

# Summary of posterior estimates:
# - Fixed effects (trait-specific intercepts and growth slopes)
# - Variance components for phylo, haplo, and units
summary(PP_multinomial_growth)

# Autocorrelation diagnostics for posterior chains:
# - VCV: variance components by PP category
# - Sol: fixed effects (trait and trait:scalegrowth terms)
autocorr.diag(PP_multinomial_growth$VCV)
autocorr.diag(PP_multinomial_growth$Sol)

# Custom autocorrelation plots for each variance component
# (function plot.acfs defined earlier)
plot.acfs(PP_multinomial_growth$VCV)

# ------------------------------------------------------------
# Extraction of category-specific growth effects
#
# Goal:
# - Isolate posterior summaries for the interaction terms
#   (trait:scalegrowth)
# - Build a tidy data frame with posterior means and 95% CIs
#   for each PP category
# ------------------------------------------------------------

# Extract interaction coefficients (scalegrowth * trait) from PP_multinomial_growth
summary_PP_growth <- summary(PP_multinomial_growth)$solutions
coef_data_PP_growth <- as.data.frame(summary_PP_growth) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scalegrowth", trait)) %>%
  mutate(Trait = gsub("trait|:scalegrowth", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Define plotting order of PP categories
# (preserves the order in which categories appear in the model output)
coef_data_PP_growth$Trait <- factor(coef_data_PP_growth$Trait, levels = unique(coef_data_PP_growth$Trait))

# ------------------------------------------------------------
# Subset: breeds with causes of death + reproductive investment data
# (64 breeds after removing missing values)
#
# PURPOSE:
# - Build a dataset including mortality outcomes (OS and PP),
#   body weight (scaleweight), and reproductive investment (scalebwkgxls)
# - Convert proportional mortality into breed-level death counts
# - Match breeds with available genetic similarity matrices
# - Construct SVD-based random-effect covariates (phylo and haplo)
# - Fit multinomial models including both weight and reproductive investment
# ------------------------------------------------------------

# Subset: only breeds with causes of death + reproductive investment data (64 breeds)
dfbwkgxls <- subset(data, select = c(Garamszegi2020,
                                     n,
                                     scaleweight,      # Body weight (standardized)
                                     scalebwkgxls,    # Reproductive investment (standardized)
                                     OScardiovascular,
                                     OSdermatologic,
                                     OSendocrine,
                                     OSgastrointestinal,
                                     OShematopoietic,
                                     OShepatic,
                                     OSmusculoskeletal,
                                     OSneurologic,
                                     OSophthalmologic,
                                     OSrespiratory,
                                     OSurogenital,
                                     OSunclear,
                                     PPcongenital,
                                     PPdegenenerative,
                                     PPinfectious,
                                     PPinflammatory,
                                     PPmetabolic,
                                     PPneoplastic,
                                     PPtoxic,
                                     PPtraumatic,
                                     PPvascular,
                                     PPunclear
))

# Remove rows with missing values to ensure complete cases for modeling
dfbwkgxls <- na.omit(dfbwkgxls)

# Ensure breed identifiers are character strings and set as row names
dfbwkgxls$Garamszegi2020 <- as.character(dfbwkgxls$Garamszegi2020)
row.names(dfbwkgxls) <- dfbwkgxls$Garamszegi2020

# Inspect correlation between standardized body weight and reproductive investment
# (used to assess collinearity between predictors)
cor(dfbwkgxls$scaleweight, dfbwkgxls$scalebwkgxls)

# ------------------------------------------------------------
# Convert OS proportional mortality into breed-level death counts
# count = round(proportion * n)
# ------------------------------------------------------------

# Number of dogs per breed for each cause of death (OS)
dfbwkgxls$OScv_N   <- round(dfbwkgxls$OScardiovascular * dfbwkgxls$n)
dfbwkgxls$OSderm_N <- round(dfbwkgxls$OSdermatologic * dfbwkgxls$n)
dfbwkgxls$OSendo_N <- round(dfbwkgxls$OSendocrine * dfbwkgxls$n)
dfbwkgxls$OSgi_N   <- round(dfbwkgxls$OSgastrointestinal * dfbwkgxls$n)
dfbwkgxls$OShem_N  <- round(dfbwkgxls$OShematopoietic * dfbwkgxls$n)
dfbwkgxls$OShep_N  <- round(dfbwkgxls$OShepatic * dfbwkgxls$n)
dfbwkgxls$OSms_N   <- round(dfbwkgxls$OSmusculoskeletal * dfbwkgxls$n)
dfbwkgxls$OSneuro_N<- round(dfbwkgxls$OSneurologic * dfbwkgxls$n)
dfbwkgxls$OSophth_N<- round(dfbwkgxls$OSophthalmologic * dfbwkgxls$n)
dfbwkgxls$OSresp_N <- round(dfbwkgxls$OSrespiratory * dfbwkgxls$n)
dfbwkgxls$OSuro_N  <- round(dfbwkgxls$OSurogenital * dfbwkgxls$n)
dfbwkgxls$OSuncl_N <- round(dfbwkgxls$OSunclear * dfbwkgxls$n)

# ------------------------------------------------------------
# Convert PP proportional mortality into breed-level death counts
# count = round(proportion * n)
# ------------------------------------------------------------

# Number of dogs per breed for each disease (PP)
dfbwkgxls$PPcongen_N <- round(dfbwkgxls$PPcongenital * dfbwkgxls$n)
dfbwkgxls$PPdegen_N  <- round(dfbwkgxls$PPdegenenerative * dfbwkgxls$n)
dfbwkgxls$PPinfect_N <- round(dfbwkgxls$PPinfectious * dfbwkgxls$n)
dfbwkgxls$PPinflam_N <- round(dfbwkgxls$PPinflammatory * dfbwkgxls$n)
dfbwkgxls$PPmetab_N  <- round(dfbwkgxls$PPmetabolic * dfbwkgxls$n)
dfbwkgxls$PPneopl_N  <- round(dfbwkgxls$PPneoplastic * dfbwkgxls$n)
dfbwkgxls$PPtoxic_N  <- round(dfbwkgxls$PPtoxic * dfbwkgxls$n)
dfbwkgxls$PPtraum_N  <- round(dfbwkgxls$PPtraumatic * dfbwkgxls$n)
dfbwkgxls$PPvasc_N   <- round(dfbwkgxls$PPvascular * dfbwkgxls$n)
dfbwkgxls$PPuncl_N   <- round(dfbwkgxls$PPunclear * dfbwkgxls$n)

# ------------------------------------------------------------
# Match breeds with SNP-based similarity matrix (common ancestry)
# ------------------------------------------------------------

# Genetic similarity matrix (common ancestry)
snps_dfbwkgxls_intersect <- intersect(dfbwkgxls$Garamszegi2020, snps$Garamszegi2020)
snps_dfbwkgxls_common <- snps[snps$Garamszegi2020 %in% snps_dfbwkgxls_intersect, ]
dfbwkgxls_filtered <- dfbwkgxls[dfbwkgxls$Garamszegi2020 %in% snps_dfbwkgxls_intersect, ]

# Transpose genetic data so breeds are indexed consistently
snps_dfbwkgxls_common_transposed <- as.data.frame(t(snps_dfbwkgxls_common))
colnames(snps_dfbwkgxls_common_transposed) <- snps_dfbwkgxls_common_transposed[1,]
snps_dfbwkgxls_common_transposed <- snps_dfbwkgxls_common_transposed[-1,]

snps2 <- rownames_to_column(snps_dfbwkgxls_common_transposed, var = "Garamszegi2020")
row.names(snps2) <- snps2$Garamszegi2020
ancestry_bwkgxls <- subset(snps2, select = -Garamszegi2020)

# ------------------------------------------------------------
# Match breeds with haplotype-based similarity matrix (haplotype sharing)
# ------------------------------------------------------------

# Gene flow matrix
haplotypes_dfbwkgxls_intersect <- intersect(dfbwkgxls$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_dfbwkgxls_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_dfbwkgxls_intersect, ]
dfbwkgxls_filtered <- dfbwkgxls[dfbwkgxls$Garamszegi2020 %in% haplotypes_dfbwkgxls_intersect, ]

# Transpose haplotype data so breeds are indexed consistently
haplotypes_dfbwkgxls_common_transposed <- as.data.frame(t(haplotypes_dfbwkgxls_common))
colnames(haplotypes_dfbwkgxls_common_transposed) <- haplotypes_dfbwkgxls_common_transposed[1,]
haplotypes_dfbwkgxls_common_transposed <- haplotypes_dfbwkgxls_common_transposed[-1,]

haplotypes2 <- rownames_to_column(haplotypes_dfbwkgxls_common_transposed, var = "Garamszegi2020")
row.names(haplotypes2) <- haplotypes2$Garamszegi2020
geneflow_bwkgxls <- subset(haplotypes2, select = -Garamszegi2020)

# ------------------------------------------------------------
# Include genetic data as random factors via SVD transformation
# ------------------------------------------------------------

# Convert matrices to numeric and compute SVD for common ancestry
ancestry_bwkgxls <- data.frame(sapply(ancestry_bwkgxls, as.numeric))
Snpsvd1 <- svd(as.matrix(ancestry_bwkgxls))
dfbwkgxls$phylo <- Snpsvd1$v %*% (t(Snpsvd1$u) * sqrt(Snpsvd1$d))

# Convert matrices to numeric and compute SVD for haplotype sharing
geneflow_bwkgxls <- data.frame(sapply(geneflow_bwkgxls, as.numeric))
Hapsvd2 <- svd(as.matrix(geneflow_bwkgxls))
dfbwkgxls$haplo <- Hapsvd2$v %*% (t(Hapsvd2$u) * sqrt(Hapsvd2$d))

# ------------------------------------------------------------
# Multinomial model including body weight and reproductive investment (OS)
# ------------------------------------------------------------

# Pool selected OS categories into a single "others" category
dfbwkgxls$OSothers_N <- rowSums(dfbwkgxls[, c(27, 31, 34, 37)])

# Fit OS multinomial model with category-specific effects of reproductive investment
# and body weight
OS_multinomial_bwkgxls <- MCMCglmm(cbind(OScv_N,
                                         OSendo_N,
                                         OSgi_N, 
                                         OShem_N, 
                                         OSms_N, 
                                         OSneuro_N,
                                         OSresp_N, 
                                         OSuro_N,
                                         OSothers_N) ~ trait -1 + trait:scalebwkgxls + trait:scaleweight,
                                   random = ~idh(trait):phylo + idh(trait):haplo, 
                                   rcov = ~idh(trait):units, 
                                   data = dfbwkgxls, 
                                   family = "multinomial9",
                                   prior = prior5,
                                   thin = 400,
                                   burnin = 100000,
                                   nitt = 1800000,
                                   verbose = TRUE)

# plot(OS_multinomial_bwkgxls)

# Posterior summaries and MCMC diagnostics
summary(OS_multinomial_bwkgxls)
autocorr.diag(OS_multinomial_bwkgxls$VCV)
autocorr.diag(OS_multinomial_bwkgxls$Sol)
plot.acfs(OS_multinomial_bwkgxls$VCV)

# ------------------------------------------------------------
# Extract category-specific reproductive investment effects (OS)
# ------------------------------------------------------------

# Extract interaction coefficients (scalebwkgxls * trait) from OS_multinomial_bwkgxls
summary_OS_bwkgxls <- summary(OS_multinomial_bwkgxls)$solutions

coef_data_OS_bwkgxls <- as.data.frame(summary_OS_bwkgxls) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scalebwkgxls", trait)) %>%
  mutate(Trait = gsub("trait|:scalebwkgxls", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# ------------------------------------------------------------
# Multinomial model including body weight and reproductive investment (PP)
# ------------------------------------------------------------

# Pool selected PP categories into a single "others" category
dfbwkgxls$PPothers_N <- rowSums(dfbwkgxls[, c(44, 46, 47)])

# Fit PP multinomial model with category-specific effects of reproductive investment
# and body weight
PP_multinomial_bwkgxls = MCMCglmm(cbind(PPcongen_N, 
                                        PPdegen_N, 
                                        PPinflam_N, 
                                        PPinfect_N,
                                        PPmetab_N, 
                                        PPneopl_N,
                                        PPtraum_N,
                                        PPothers_N) ~ trait -1 + trait:scalebwkgxls + trait:scaleweight, 
                                  random = ~idh(trait):phylo + idh(trait):haplo, 
                                  rcov = ~idh(trait):units, 
                                  data = dfbwkgxls, 
                                  family = "multinomial8",
                                  prior = prior3,
                                  thin = 400,
                                  burnin = 100000,
                                  nitt = 1800000,
                                  verbose = TRUE)

# plot(PP_multinomial_bwkgxls)

# Posterior summaries and MCMC diagnostics
summary(PP_multinomial_bwkgxls)
autocorr.diag(PP_multinomial_bwkgxls$VCV)
autocorr.diag(PP_multinomial_bwkgxls$Sol)
plot.acfs(PP_multinomial_bwkgxls$VCV)

# ------------------------------------------------------------
# Extract category-specific reproductive investment effects (PP)
# ------------------------------------------------------------

# Extract interaction coefficients (scalebwkgxls * trait) from PP_multinomial_bwkgxls
summary_PP_bwkgxls <- summary(PP_multinomial_bwkgxls)$solutions

coef_data_PP_bwkgxls <- as.data.frame(summary_PP_bwkgxls) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scalebwkgxls", trait)) %>%
  mutate(Trait = gsub("trait|:scalebwkgxls", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Define plotting order of PP categories
# (preserves the order in which categories appear in the model output)
coef_data_PP_bwkgxls$Trait <- factor(coef_data_PP_bwkgxls$Trait, levels = unique(coef_data_PP_bwkgxls$Trait))

###############################################################################
# ------------------------------------------------------------
# Subset: breeds with causes of death + activity data
# (49 breeds after removing missing values)
#
# PURPOSE:
# - Build a dataset including mortality outcomes (OS and PP),
#   body weight (scaleweight), and activity (scaleactivityCareau2010)
# - Convert proportional mortality into breed-level death counts
# - Match breeds with available genetic similarity matrices
# - Construct SVD-based random-effect covariates (phylo and haplo)
# - Fit multinomial models including both weight and activity as predictors
# ------------------------------------------------------------

# Subset: only breeds with causes of death + activity data (49 breeds)  
dfactivity <- subset(data, select=c(Garamszegi2020,
                                    n,
                                    scaleweight,                  # Body weight (standardized)
                                    scaleactivityCareau2010,      # Activity (standardized)
                                    OScardiovascular,
                                    OSdermatologic,
                                    OSendocrine,
                                    OSgastrointestinal,
                                    OShematopoietic,
                                    OShepatic,
                                    OSmusculoskeletal,
                                    OSneurologic,
                                    OSophthalmologic,
                                    OSrespiratory,
                                    OSurogenital,
                                    OSunclear,
                                    PPcongenital,
                                    PPdegenenerative,
                                    PPinfectious,
                                    PPinflammatory,
                                    PPmetabolic,
                                    PPneoplastic,
                                    PPtoxic,
                                    PPtraumatic,
                                    PPvascular,
                                    PPunclear
))

# Remove rows with missing values to ensure complete cases for modeling
dfactivity <- na.omit(dfactivity)

# Ensure breed identifiers are character strings and set as row names
dfactivity$Garamszegi2020 <- as.character(dfactivity$Garamszegi2020)
row.names(dfactivity) <- dfactivity$Garamszegi2020

# Check correlation between standardized body weight and standardized activity
# (used to assess potential collinearity between predictors)
cor(dfactivity$scaleweight, dfactivity$scaleactivityCareau2010)

# ------------------------------------------------------------
# Convert OS proportional mortality into breed-level death counts
# count = round(proportion * n)
# ------------------------------------------------------------

# Number of dogs per breed for each cause of death (OS)
dfactivity$OScv_N   <- round(dfactivity$OScardiovascular * dfactivity$n)
dfactivity$OSderm_N <- round(dfactivity$OSdermatologic * dfactivity$n)
dfactivity$OSendo_N <- round(dfactivity$OSendocrine * dfactivity$n)
dfactivity$OSgi_N   <- round(dfactivity$OSgastrointestinal * dfactivity$n)
dfactivity$OShem_N  <- round(dfactivity$OShematopoietic * dfactivity$n)
dfactivity$OShep_N  <- round(dfactivity$OShepatic * dfactivity$n)
dfactivity$OSms_N   <- round(dfactivity$OSmusculoskeletal * dfactivity$n)
dfactivity$OSneuro_N<- round(dfactivity$OSneurologic * dfactivity$n)
dfactivity$OSophth_N<- round(dfactivity$OSophthalmologic * dfactivity$n)
dfactivity$OSresp_N <- round(dfactivity$OSrespiratory * dfactivity$n)
dfactivity$OSuro_N  <- round(dfactivity$OSurogenital * dfactivity$n)
dfactivity$OSuncl_N <- round(dfactivity$OSunclear * dfactivity$n)

# ------------------------------------------------------------
# Convert PP proportional mortality into breed-level death counts
# count = round(proportion * n)
# ------------------------------------------------------------

# Number of dogs per breed for each disease (PP)
dfactivity$PPcongen_N <- round(dfactivity$PPcongenital * dfactivity$n)
dfactivity$PPdegen_N  <- round(dfactivity$PPdegenenerative * dfactivity$n)
dfactivity$PPinfect_N <- round(dfactivity$PPinfectious * dfactivity$n)
dfactivity$PPinflam_N <- round(dfactivity$PPinflammatory * dfactivity$n)
dfactivity$PPmetab_N  <- round(dfactivity$PPmetabolic * dfactivity$n)
dfactivity$PPneopl_N  <- round(dfactivity$PPneoplastic * dfactivity$n)
dfactivity$PPtoxic_N  <- round(dfactivity$PPtoxic * dfactivity$n)
dfactivity$PPtraum_N  <- round(dfactivity$PPtraumatic * dfactivity$n)
dfactivity$PPvasc_N   <- round(dfactivity$PPvascular * dfactivity$n)
dfactivity$PPuncl_N   <- round(dfactivity$PPunclear * dfactivity$n)

# ------------------------------------------------------------
# Match breeds with SNP-based similarity matrix (common ancestry)
# ------------------------------------------------------------

# Genetic similarity matrix (common ancestry)
snps_dfactivity_intersect <- intersect(dfactivity$Garamszegi2020, snps$Garamszegi2020)
snps_dfactivity_common <- snps[snps$Garamszegi2020 %in% snps_dfactivity_intersect, ]
dfactivity_filtered <- dfactivity[dfactivity$Garamszegi2020 %in% snps_dfactivity_intersect, ]

# Transpose genetic data so breeds are indexed consistently
snps_dfactivity_common_transposed <- as.data.frame(t(snps_dfactivity_common))
colnames(snps_dfactivity_common_transposed) <- snps_dfactivity_common_transposed[1,]
snps_dfactivity_common_transposed <- snps_dfactivity_common_transposed[-1,]

snps2 <- rownames_to_column(snps_dfactivity_common_transposed, var = "Garamszegi2020")
row.names(snps2) <- snps2$Garamszegi2020
ancestry_activity <- subset(snps2, select = -Garamszegi2020)

# ------------------------------------------------------------
# Match breeds with haplotype-based similarity matrix (haplotype sharing)
# ------------------------------------------------------------

# Gene flow matrix (haplotype sharing)
haplotypes_dfactivity_intersect <- intersect(dfactivity$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_dfactivity_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_dfactivity_intersect, ]
dfactivity_filtered <- dfactivity[dfactivity$Garamszegi2020 %in% haplotypes_dfactivity_intersect, ]

# Transpose haplotype data so breeds are indexed consistently
haplotypes_dfactivity_common_transposed <- as.data.frame(t(haplotypes_dfactivity_common))
colnames(haplotypes_dfactivity_common_transposed) <- haplotypes_dfactivity_common_transposed[1,]
haplotypes_dfactivity_common_transposed <- haplotypes_dfactivity_common_transposed[-1,]

haplotypes2 <- rownames_to_column(haplotypes_dfactivity_common_transposed, var = "Garamszegi2020")
row.names(haplotypes2) <- haplotypes2$Garamszegi2020
geneflow_activity <- subset(haplotypes2, select = -Garamszegi2020)

# ------------------------------------------------------------
# Include genetic data as random factors via SVD transformation
# ------------------------------------------------------------

# Convert matrices to numeric and compute SVD for common ancestry
ancestry_activity <- data.frame(sapply(ancestry_activity, as.numeric))
Snpsvd1 <- svd(as.matrix(ancestry_activity))
dfactivity$phylo <- Snpsvd1$v %*% (t(Snpsvd1$u) * sqrt(Snpsvd1$d))

# Convert matrices to numeric and compute SVD for haplotype sharing
geneflow_activity <- data.frame(sapply(geneflow_activity, as.numeric))
Hapsvd2 <- svd(as.matrix(geneflow_activity))
dfactivity$haplo <- Hapsvd2$v %*% (t(Hapsvd2$u) * sqrt(Hapsvd2$d))

# ------------------------------------------------------------
# Multinomial model including body weight and activity (OS)
# ------------------------------------------------------------

# Pool selected OS categories into a single "others" category
dfactivity$OSothers_N <- rowSums(dfactivity[, c(27, 31, 34, 37)])

# Fit OS multinomial model with category-specific effects of activity and body weight:
# - Fixed effects: trait-specific intercepts plus trait-specific slopes for
#   activity (scaleactivityCareau2010) and weight (scaleweight)
# - Random effects: genetic influence via phylo (common ancestry) and haplo (haplotype sharing)
OS_multinomial_activity <- MCMCglmm(cbind(OScv_N,
                                          OSendo_N,
                                          OSgi_N, 
                                          OShem_N, 
                                          OSms_N, 
                                          OSneuro_N,
                                          OSresp_N, 
                                          OSuro_N,
                                          OSothers_N) ~ trait -1 + trait:scaleactivityCareau2010 + trait:scaleweight,
                                    random = ~idh(trait):phylo + idh(trait):haplo, 
                                    rcov = ~idh(trait):units, 
                                    data = dfactivity, 
                                    family = "multinomial9",
                                    prior = prior5,
                                    thin = 400,
                                    burnin = 100000,
                                    nitt = 1800000,
                                    verbose = TRUE)

# Posterior summaries and MCMC diagnostics (OS model)
summary(OS_multinomial_activity)
autocorr.diag(OS_multinomial_activity$VCV)
autocorr.diag(OS_multinomial_activity$Sol)
plot.acfs(OS_multinomial_activity$VCV)

# ------------------------------------------------------------
# Extract category-specific activity effects (OS)
# ------------------------------------------------------------

# Extract interaction coefficients (scaleactivityCareau2010 * trait) from OS_multinomial_activity
summary_OS_activity <- summary(OS_multinomial_activity)$solutions

coef_data_OS_activity <- as.data.frame(summary_OS_activity) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scaleactivityCareau2010", trait)) %>%
  mutate(Trait = gsub("trait|:scaleactivityCareau2010", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Define plotting order of OS categories
# (preserves the order in which categories appear in the model output)
coef_data_OS_activity$Trait <- factor(coef_data_OS_activity$Trait, levels = unique(coef_data_OS_activity$Trait))

###############################################################################
# ------------------------------------------------------------
# Multinomial model including body weight and activity (PP)
# ------------------------------------------------------------

# Pool selected PP categories into a single "others" category
dfactivity$PPothers_N <- rowSums(dfactivity[, c(44, 46, 47)])

# Fit PP multinomial model with category-specific effects of activity and body weight:
# - Fixed effects: trait-specific intercepts plus trait-specific slopes for
#   activity (scaleactivityCareau2010) and weight (scaleweight)
# - Random effects: genetic influence via phylo (common ancestry) and haplo (haplotype sharing)
PP_multinomial_activity <- MCMCglmm(cbind(PPcongen_N, 
                                          PPdegen_N,
                                          PPinfect_N,
                                          PPinflam_N, 
                                          PPmetab_N, 
                                          PPneopl_N,
                                          PPtraum_N,
                                          PPothers_N) ~ trait -1 + trait:scaleactivityCareau2010 + trait:scaleweight, 
                                    random = ~idh(trait):phylo + idh(trait):haplo, 
                                    rcov = ~idh(trait):units, 
                                    data = dfactivity, 
                                    family = "multinomial8",
                                    prior = prior3,
                                    thin = 400,
                                    burnin = 100000,
                                    nitt = 1800000,
                                    verbose = TRUE)

# Posterior summary (PP model)
summary(PP_multinomial_activity)

# ------------------------------------------------------------
# Extract category-specific activity effects (PP)
# ------------------------------------------------------------

# Extract interaction coefficients (scaleactivityCareau2010 * trait) from PP_multinomial_activity
summary_PP_activity <- summary(PP_multinomial_activity)$solutions

coef_data_PP_activity <- as.data.frame(summary_PP_activity) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scaleactivityCareau2010", trait)) %>%
  mutate(Trait = gsub("trait|:scaleactivityCareau2010", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Define plotting order of PP categories
# (preserves the order in which categories appear in the model output)
coef_data_PP_activity$Trait <- factor(coef_data_PP_activity$Trait, levels = unique(coef_data_PP_activity$Trait))

###############################################################################
# ------------------------------------------------------------
# Subset: breeds with causes of death + aggressiveness data
# (49 breeds after removing missing values)
#
# PURPOSE:
# - Build a dataset including mortality outcomes (OS and PP)
#   and breed-level aggressiveness scores (Careau et al. 2010)
# - Convert proportional mortality into breed-level death counts
# - Match breeds with available genetic similarity matrices
# - Construct SVD-based random-effect covariates (phylo and haplo)
# - Fit multinomial models with category-specific aggressiveness effects
# ------------------------------------------------------------

# Subset: only breeds with causes of death + aggressiveness data (49 breeds)  
dfaggressivenessCareau2010 <- subset(data, select=c(Garamszegi2020,
                                                    n,
                                                    scaleaggressivenessCareau2010,      # Standardized aggressiveness
                                                    OScardiovascular,
                                                    OSdermatologic,
                                                    OSendocrine,
                                                    OSgastrointestinal,
                                                    OShematopoietic,
                                                    OShepatic,
                                                    OSmusculoskeletal,
                                                    OSneurologic,
                                                    OSophthalmologic,
                                                    OSrespiratory,
                                                    OSurogenital,
                                                    OSunclear,
                                                    PPcongenital,
                                                    PPdegenenerative,
                                                    PPinfectious,
                                                    PPinflammatory,
                                                    PPmetabolic,
                                                    PPneoplastic,
                                                    PPtoxic,
                                                    PPtraumatic,
                                                    PPvascular,
                                                    PPunclear
))

# Remove rows with missing values to ensure complete cases
dfaggressivenessCareau2010 <- na.omit(dfaggressivenessCareau2010)

# Ensure breed identifiers are character strings and set as row names
dfaggressivenessCareau2010$Garamszegi2020 <- as.character(dfaggressivenessCareau2010$Garamszegi2020)
row.names(dfaggressivenessCareau2010) <- dfaggressivenessCareau2010$Garamszegi2020

# ------------------------------------------------------------
# Convert OS proportional mortality into breed-level death counts
# count = round(proportion * n)
# ------------------------------------------------------------

# Number of dogs per breed for each cause of death (OS)
dfaggressivenessCareau2010$OScv_N   <- round(dfaggressivenessCareau2010$OScardiovascular * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSderm_N <- round(dfaggressivenessCareau2010$OSdermatologic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSendo_N <- round(dfaggressivenessCareau2010$OSendocrine * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSgi_N   <- round(dfaggressivenessCareau2010$OSgastrointestinal * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OShem_N  <- round(dfaggressivenessCareau2010$OShematopoietic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OShep_N  <- round(dfaggressivenessCareau2010$OShepatic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSms_N   <- round(dfaggressivenessCareau2010$OSmusculoskeletal * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSneuro_N<- round(dfaggressivenessCareau2010$OSneurologic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSophth_N<- round(dfaggressivenessCareau2010$OSophthalmologic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSresp_N <- round(dfaggressivenessCareau2010$OSrespiratory * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSuro_N  <- round(dfaggressivenessCareau2010$OSurogenital * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$OSuncl_N <- round(dfaggressivenessCareau2010$OSunclear * dfaggressivenessCareau2010$n)

# ------------------------------------------------------------
# Convert PP proportional mortality into breed-level death counts
# count = round(proportion * n)
# ------------------------------------------------------------

# Number of dogs per breed for each disease (PP)
dfaggressivenessCareau2010$PPcongen_N <- round(dfaggressivenessCareau2010$PPcongenital * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPdegen_N  <- round(dfaggressivenessCareau2010$PPdegenenerative * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPinfect_N <- round(dfaggressivenessCareau2010$PPinfectious * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPinflam_N <- round(dfaggressivenessCareau2010$PPinflammatory * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPmetab_N  <- round(dfaggressivenessCareau2010$PPmetabolic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPneopl_N  <- round(dfaggressivenessCareau2010$PPneoplastic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPtoxic_N  <- round(dfaggressivenessCareau2010$PPtoxic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPtraum_N  <- round(dfaggressivenessCareau2010$PPtraumatic * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPvasc_N   <- round(dfaggressivenessCareau2010$PPvascular * dfaggressivenessCareau2010$n)
dfaggressivenessCareau2010$PPuncl_N   <- round(dfaggressivenessCareau2010$PPunclear * dfaggressivenessCareau2010$n)

# ------------------------------------------------------------
# Match breeds with SNP-based similarity matrix (common ancestry)
# ------------------------------------------------------------

# Genetic similarity matrix based on shared SNPs (common ancestry)
snps_dfaggressiveness_intersect <- intersect(dfaggressivenessCareau2010$Garamszegi2020, snps$Garamszegi2020)
snps_dfaggressiveness_common <- snps[snps$Garamszegi2020 %in% snps_dfaggressiveness_intersect, ]
dfaggressiveness_filtered <- dfaggressivenessCareau2010[dfaggressivenessCareau2010$Garamszegi2020 %in% snps_dfaggressiveness_intersect, ]

# Transpose genetic matrix to align breeds as rows
snps_dfaggressiveness_common_transposed <- as.data.frame(t(snps_dfaggressiveness_common))
colnames(snps_dfaggressiveness_common_transposed) <- snps_dfaggressiveness_common_transposed[1,]
snps_dfaggressiveness_common_transposed <- snps_dfaggressiveness_common_transposed[-1,]

snps2 <- rownames_to_column(snps_dfaggressiveness_common_transposed, var = "Garamszegi2020")
snps_dfaggressiveness_intersect <- intersect(dfaggressivenessCareau2010$Garamszegi2020, snps2$Garamszegi2020)
snps_dfaggressiveness_common <- snps2[snps2$Garamszegi2020 %in% snps_dfaggressiveness_intersect, ]
row.names(snps_dfaggressiveness_common) <- snps_dfaggressiveness_common$Garamszegi2020
ancestry_aggressiveness <- subset(snps_dfaggressiveness_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# Match breeds with haplotype-based similarity matrix (gene flow)
# ------------------------------------------------------------

# Gene flow matrix based on shared haplotypes
haplotypes_dfaggressiveness_intersect <- intersect(dfaggressivenessCareau2010$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_dfaggressiveness_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_dfaggressiveness_intersect, ]
dfaggressiveness_filtered <- dfaggressivenessCareau2010[dfaggressivenessCareau2010$Garamszegi2020 %in% haplotypes_dfaggressiveness_intersect, ]

# Transpose haplotype matrix to align breeds as rows
haplotypes_dfaggressiveness_common_transposed <- as.data.frame(t(haplotypes_dfaggressiveness_common))
colnames(haplotypes_dfaggressiveness_common_transposed) <- haplotypes_dfaggressiveness_common_transposed[1,]
haplotypes_dfaggressiveness_common_transposed <- haplotypes_dfaggressiveness_common_transposed[-1,]

haplotypes2 <- rownames_to_column(haplotypes_dfaggressiveness_common_transposed, var = "Garamszegi2020")
haplotypes_dfaggressiveness_intersect <- intersect(dfaggressivenessCareau2010$Garamszegi2020, haplotypes2$Garamszegi2020)
haplotypes_dfaggressiveness_common <- haplotypes2[haplotypes2$Garamszegi2020 %in% haplotypes_dfaggressiveness_intersect, ]
row.names(haplotypes_dfaggressiveness_common) <- haplotypes_dfaggressiveness_common$Garamszegi2020
geneflow_aggressiveness <- subset(haplotypes_dfaggressiveness_common, select = -Garamszegi2020)

# ------------------------------------------------------------
# Include genetic information as random effects via SVD
# ------------------------------------------------------------

# Convert common ancestry matrix to numeric and apply SVD
ancestry_aggressiveness <- data.frame(sapply(ancestry_aggressiveness, function(x) as.numeric(as.character(x))))
Snpsvd1 <- svd(as.matrix(ancestry_aggressiveness))
Snpsvd1 <- Snpsvd1$v %*% (t(Snpsvd1$u) * sqrt(Snpsvd1$d))
dfaggressivenessCareau2010$phylo <- Snpsvd1

# Convert haplotype-sharing matrix to numeric and apply SVD
geneflow_aggressiveness <- data.frame(sapply(geneflow_aggressiveness, function(x) as.numeric(as.character(x))))
Hapsvd2 <- svd(as.matrix(geneflow_aggressiveness))
Hapsvd2 <- Hapsvd2$v %*% (t(Hapsvd2$u) * sqrt(Hapsvd2$d))
dfaggressivenessCareau2010$haplo <- Hapsvd2

# Ensure row names and random-effect covariates are correctly assigned
row.names(dfaggressivenessCareau2010) <- dfaggressivenessCareau2010$Garamszegi2020
dfaggressivenessCareau2010$phylo <- Snpsvd1
dfaggressivenessCareau2010$haplo <- Hapsvd2

# ------------------------------------------------------------
# Multinomial model including aggressiveness (OS)
# ------------------------------------------------------------

# Pool selected OS categories into an "others" category
dfaggressivenessCareau2010$OSothers_N <- rowSums(dfaggressivenessCareau2010[, c(27, 31, 34, 37)])

# Fit OS multinomial model with category-specific aggressiveness effects
OS_multinomial_aggressiveness <- MCMCglmm(cbind(OScv_N,
                                                OSendo_N,
                                                OSgi_N, 
                                                OShem_N, 
                                                OSms_N, 
                                                OSneuro_N,
                                                OSresp_N, 
                                                OSuro_N,
                                                OSothers_N) ~ trait-1+trait:scaleaggressivenessCareau2010,
                                          random = ~idh(trait):phylo + idh(trait):haplo, 
                                          rcov = ~idh(trait):units, 
                                          data = dfaggressivenessCareau2010, 
                                          family = "multinomial9",
                                          prior = prior5,
                                          thin = 400,
                                          burnin = 100000,
                                          nitt = 1800000,
                                          verbose = TRUE)

# Posterior summary of the OS aggressiveness model
summary(OS_multinomial_aggressiveness)

# ------------------------------------------------------------
# Multinomial model including aggressiveness (PP)
# ------------------------------------------------------------

# Pool selected PP categories into an "others" category
dfaggressivenessCareau2010$PPothers_N <- rowSums(dfaggressivenessCareau2010[, c(44, 46, 47)])

# Fit PP multinomial model with category-specific aggressiveness effects
PP_multinomial_aggressiveness <- MCMCglmm(cbind(PPinfect_N, 
                                                PPcongen_N, 
                                                PPdegen_N, 
                                                PPinflam_N, 
                                                PPmetab_N, 
                                                PPneopl_N,
                                                PPtraum_N,
                                                PPothers_N) ~ trait-1+trait:scaleaggressivenessCareau2010, 
                                          random = ~idh(trait):phylo + idh(trait):haplo, 
                                          rcov = ~idh(trait):units, 
                                          data = dfaggressivenessCareau2010, 
                                          family = "multinomial8",
                                          prior = prior3,
                                          thin = 400,
                                          burnin = 100000,
                                          nitt = 1800000,
                                          verbose = TRUE)

# Posterior summary of the PP aggressiveness model
summary(PP_multinomial_aggressiveness)

###############################################################################
# ------------------------------------------------------------
# Subset: breeds with causes of death + trainability data
#
# PURPOSE:
# - Build a dataset including mortality outcomes (OS and PP)
#   and breed-level trainability scores (Careau et al. 2010)
# - Convert proportional mortality into breed-level death counts
# - Match breeds with available genetic similarity matrices
# - Construct SVD-based random-effect covariates (phylo and haplo)
# - Fit multinomial models with category-specific trainability
#   effects while controlling for body weight
# ------------------------------------------------------------

# Subset: only breeds with causes of death + trainability data
dftrainability <- subset(data, select=c(Garamszegi2020,
                                        n,
                                        scaleweight,               # Standardized adult body weight
                                        scaletrainabilityCareau2010,  # Standardized trainability score
                                        OScardiovascular,
                                        OSdermatologic,
                                        OSendocrine,
                                        OSgastrointestinal,
                                        OShematopoietic,
                                        OShepatic,
                                        OSmusculoskeletal,
                                        OSneurologic,
                                        OSophthalmologic,
                                        OSrespiratory,
                                        OSurogenital,
                                        OSunclear,
                                        PPcongenital,
                                        PPdegenenerative,
                                        PPinfectious,
                                        PPinflammatory,
                                        PPmetabolic,
                                        PPneoplastic,
                                        PPtoxic,
                                        PPtraumatic,
                                        PPvascular,
                                        PPunclear
))

# Remove rows with missing values to ensure complete cases
dftrainability <- na.omit(dftrainability)

# Ensure breed identifiers are character strings and assign them as row names
dftrainability$Garamszegi2020 <- as.character(dftrainability$Garamszegi2020)
row.names(dftrainability) <- dftrainability$Garamszegi2020

# ------------------------------------------------------------
# Convert OS proportional mortality into breed-level death counts
# count = round(proportion * total sample size)
# ------------------------------------------------------------

# Number of dogs per breed for each organ system (OS)
dftrainability$OScv_N   <- round(dftrainability$OScardiovascular * dftrainability$n)
dftrainability$OSderm_N <- round(dftrainability$OSdermatologic * dftrainability$n)
dftrainability$OSendo_N <- round(dftrainability$OSendocrine * dftrainability$n)
dftrainability$OSgi_N   <- round(dftrainability$OSgastrointestinal * dftrainability$n)
dftrainability$OShem_N  <- round(dftrainability$OShematopoietic * dftrainability$n)
dftrainability$OShep_N  <- round(dftrainability$OShepatic * dftrainability$n)
dftrainability$OSms_N   <- round(dftrainability$OSmusculoskeletal * dftrainability$n)
dftrainability$OSneuro_N<- round(dftrainability$OSneurologic * dftrainability$n)
dftrainability$OSophth_N<- round(dftrainability$OSophthalmologic * dftrainability$n)
dftrainability$OSresp_N <- round(dftrainability$OSrespiratory * dftrainability$n)
dftrainability$OSuro_N  <- round(dftrainability$OSurogenital * dftrainability$n)
dftrainability$OSuncl_N <- round(dftrainability$OSunclear * dftrainability$n)

# ------------------------------------------------------------
# Convert PP proportional mortality into breed-level death counts
# ------------------------------------------------------------

# Number of dogs per breed for each pathophysiological process (PP)
dftrainability$PPcongen_N <- round(dftrainability$PPcongenital * dftrainability$n)
dftrainability$PPdegen_N  <- round(dftrainability$PPdegenenerative * dftrainability$n)
dftrainability$PPinfect_N <- round(dftrainability$PPinfectious * dftrainability$n)
dftrainability$PPinflam_N <- round(dftrainability$PPinflammatory * dftrainability$n)
dftrainability$PPmetab_N  <- round(dftrainability$PPmetabolic * dftrainability$n)
dftrainability$PPneopl_N  <- round(dftrainability$PPneoplastic * dftrainability$n)
dftrainability$PPtoxic_N  <- round(dftrainability$PPtoxic * dftrainability$n)
dftrainability$PPtraum_N  <- round(dftrainability$PPtraumatic * dftrainability$n)
dftrainability$PPvasc_N   <- round(dftrainability$PPvascular * dftrainability$n)
dftrainability$PPuncl_N   <- round(dftrainability$PPunclear * dftrainability$n)

# ------------------------------------------------------------
# Match breeds with SNP-based similarity matrix (common ancestry)
# ------------------------------------------------------------

# Identify overlapping breeds between mortality data and SNP matrix
snps_dftrainability_intersect <- intersect(dftrainability$Garamszegi2020, snps$Garamszegi2020)
snps_dftrainability_common <- snps[snps$Garamszegi2020 %in% snps_dftrainability_intersect, ]
dftrainability_filtered <- dftrainability[dftrainability$Garamszegi2020 %in% snps_dftrainability_intersect, ]

# Transpose SNP matrix so breeds correspond to rows
snps_dftrainability_common_transposed <- as.data.frame(t(snps_dftrainability_common))
colnames(snps_dftrainability_common_transposed) <- snps_dftrainability_common_transposed[1,]
snps_dftrainability_common_transposed <- snps_dftrainability_common_transposed[-1,]

snps2 <- rownames_to_column(snps_dftrainability_common_transposed, var = "Garamszegi2020")
row.names(snps2) <- snps2$Garamszegi2020
ancestry_trainability <- subset(snps2, select = -Garamszegi2020)

# ------------------------------------------------------------
# Match breeds with haplotype-based similarity matrix (gene flow)
# ------------------------------------------------------------

# Identify overlapping breeds between mortality data and haplotype matrix
haplotypes_dftrainability_intersect <- intersect(dftrainability$Garamszegi2020, haplotypes$Garamszegi2020)
haplotypes_dftrainability_common <- haplotypes[haplotypes$Garamszegi2020 %in% haplotypes_dftrainability_intersect, ]
dftrainability_filtered <- dftrainability[dftrainability$Garamszegi2020 %in% haplotypes_dftrainability_intersect, ]

# Transpose haplotype matrix so breeds correspond to rows
haplotypes_dftrainability_common_transposed <- as.data.frame(t(haplotypes_dftrainability_common))
colnames(haplotypes_dftrainability_common_transposed) <- haplotypes_dftrainability_common_transposed[1,]
haplotypes_dftrainability_common_transposed <- haplotypes_dftrainability_common_transposed[-1,]

haplotypes2 <- rownames_to_column(haplotypes_dftrainability_common_transposed, var = "Garamszegi2020")
row.names(haplotypes2) <- haplotypes2$Garamszegi2020
geneflow_trainability <- subset(haplotypes2, select = -Garamszegi2020)

# ------------------------------------------------------------
# Include genetic information as random effects using SVD
# ------------------------------------------------------------

# Convert common ancestry matrix to numeric and apply SVD
ancestry_trainability <- data.frame(sapply(ancestry_trainability, as.numeric))
Snpsvd1 <- svd(as.matrix(ancestry_trainability))
dftrainability$phylo <- Snpsvd1$v %*% (t(Snpsvd1$u) * sqrt(Snpsvd1$d))

# Convert haplotype-sharing matrix to numeric and apply SVD
geneflow_trainability <- data.frame(sapply(geneflow_trainability, as.numeric))
Hapsvd2 <- svd(as.matrix(geneflow_trainability))
dftrainability$haplo <- Hapsvd2$v %*% (t(Hapsvd2$u) * sqrt(Hapsvd2$d))

# ------------------------------------------------------------
# Multinomial model including trainability and body weight (OS)
# ------------------------------------------------------------

# Pool selected OS categories into an "others" category
dftrainability$OSothers_N <- rowSums(dftrainability[, c(27, 31, 34, 37)])

# Fit OS multinomial model with category-specific trainability effects,
# controlling for body weight and genetic non-independence among breeds
OS_multinomial_trainability <- MCMCglmm(cbind(OScv_N,
                                              OSendo_N,
                                              OSgi_N, 
                                              OShem_N, 
                                              OSms_N, 
                                              OSneuro_N,
                                              OSresp_N, 
                                              OSuro_N,
                                              OSothers_N) ~ trait -1 + trait:scaletrainabilityCareau2010 + trait:scaleweight,
                                        random = ~idh(trait):phylo + idh(trait):haplo, 
                                        rcov = ~idh(trait):units, 
                                        data = dftrainability, 
                                        family = "multinomial9",
                                        prior = prior5,
                                        thin = 400,
                                        burnin = 100000,
                                        nitt = 1800000,
                                        verbose = TRUE)

# Posterior summary of the OS trainability model
summary(OS_multinomial_trainability)

# ------------------------------------------------------------
# Multinomial model including trainability and body weight (PP)
# ------------------------------------------------------------

# Pool selected PP categories into an "others" category
dftrainability$PPothers_N <- rowSums(dftrainability[, c(44, 46, 47)])

# Fit PP multinomial model with category-specific trainability effects,
# controlling for body weight and genetic structure
PP_multinomial_trainability <- MCMCglmm(cbind(PPcongen_N, 
                                              PPdegen_N,
                                              PPinfect_N,
                                              PPinflam_N, 
                                              PPmetab_N, 
                                              PPneopl_N,
                                              PPtraum_N,
                                              PPothers_N) ~ trait -1 + trait:scaletrainabilityCareau2010 + trait:scaleweight, 
                                        random = ~idh(trait):phylo + idh(trait):haplo, 
                                        rcov = ~idh(trait):units, 
                                        data = dftrainability, 
                                        family = "multinomial8",
                                        prior = prior3,
                                        thin = 400,
                                        burnin = 100000,
                                        nitt = 1800000,
                                        verbose = TRUE)

# Posterior summary of the PP trainability model
summary(PP_multinomial_trainability)

# ------------------------------------------------------------
# Extract posterior estimates for trainability × trait interaction (PP)
# ------------------------------------------------------------

# Extract posterior samples for the trainability interaction terms
summary_PP_trainability <- summary(PP_multinomial_trainability)$solutions

coef_data_PP_trainability <- as.data.frame(summary_PP_trainability) %>% 
  rownames_to_column(var = "trait") %>%
  filter(grepl("scaletrainabilityCareau2010", trait)) %>%
  mutate(Trait = gsub("trait|:scaletrainabilityCareau2010", "", trait),
         Estimate = post.mean,   
         Lower = `l-95% CI`,    
         Upper = `u-95% CI`)

# Ensure causes of death are plotted in the same order as they appear in the model
coef_data_PP_trainability$Trait <- factor(coef_data_PP_trainability$Trait, levels = unique(coef_data_PP_trainability$Trait))

#################
################# FIGURES

# Create a list to store the results (one small data.frame per VCV column/component)
OS_multinomial_results_list <- lapply(1:ncol(OS_multinomial$VCV), function(i) {
  # Compute the posterior proportion for component i at each MCMC iteration:
  # VCV[, i] / sum(VCV row) = proportion of total variance explained by that component
  OS_multinomial_variable_proportion <- OS_multinomial$VCV[, i] / rowSums(OS_multinomial$VCV)
  
  # Compute the 95% HPD interval of the posterior proportions for this component
  OS_multinomial_hpd <- HPDinterval(OS_multinomial_variable_proportion)
  
  # Compute the posterior median of the proportions for this component
  OS_multinomial_median_val <- median(OS_multinomial_variable_proportion)
  
  # Store results for this component in a one-row data.frame:
  # - Variable name (current VCV column name)
  # - Lower/Upper HPD bounds
  # - Posterior median
  data.frame(
    Variable = colnames(OS_multinomial$VCV)[i],
    Lower = OS_multinomial_hpd[1, "lower"],
    Upper = OS_multinomial_hpd[1, "upper"],
    Median = OS_multinomial_median_val
  )
})

# Combine the list of one-row data.frames into a single results table
OS_multinomial_results <- do.call(rbind, OS_multinomial_results_list)

# Print the results table (each row = one variance component; columns = HPD + median)
print(OS_multinomial_results)

# Define the desired human-readable labels for each variance component, in the exact order
# matching the columns of OS_multinomial$VCV (common ancestry, hybridization, then units)
correct_variable_names <- c(
  "Cardiovascular (common ancestry)", "Endocrine (common ancestry)", "Gastrointestinal (common ancestry)", "Hematopoietic (common ancestry)", "Musculoskeletal (common ancestry)", "Neurologic (common ancestry)",
  "Respiratory (common ancestry)", "Urogenital (common ancestry)", "Cardiovascular (hybridization)", "Endocrine (hybridization)", "Gastrointestinal (hybridization)", "Hematopoietic (hybridization)",
  "Musculoskeletal (hybridization)", "Neurologic (hybridization)", "Respiratory (hybridization)", "Urogenital (hybridization)",
  "Cardiovascular (units)", "Endocrine (units)", "Gastrointestinal (units)", "Hematopoietic (units)", "Musculoskeletal (units)",
  "Neurologic (units)", "Respiratory (units)", "Urogenital (units)"
)

# Replace the "Variable" column in the results table with the human-readable names above
# (This assumes length(correct_variable_names) == nrow(OS_multinomial_results))
OS_multinomial_results$Variable <- correct_variable_names

# Print the relabeled results table to confirm the new names were applied correctly
print(OS_multinomial_results)

# Optional: export the final table to CSV for use in Excel/manuscript tables
# write.csv(OS_multinomial_results, "OS_multinomial_results_final.csv", row.names = FALSE)

# Create a filtered version of the results table excluding residual ("units") components,
# keeping only the two genetic effects (common ancestry and hybridization)
OS_multinomial_results_filtered <- OS_multinomial_results[!grepl("units", OS_multinomial_results$Variable), ]

# Print the filtered table to confirm only non-units components remain
print(OS_multinomial_results_filtered)

# Plot HPD intervals + posterior median for ALL components (including units)
# X axis = component name; Y axis = posterior median; error bars = 95% HPD interval
# NOTE: despite the comment, this plot uses Median only (no Mean is plotted)
ggplot(OS_multinomial_results, aes(x = Variable)) +
  geom_point(aes(y = Median), color = "red", shape = 18, size = 4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
  labs(x = "Variable", y = "Value", 
       title = "HPD Intervals with Mean and Median") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Convert Variable into a factor ordered by reverse alphabetical order
# (This affects the order of categories displayed on the plot axis)
OS_multinomial_results_filtered$Variable <- factor(OS_multinomial_results_filtered$Variable, 
                                                   levels = rev(sort(unique(OS_multinomial_results_filtered$Variable))))

# Plot HPD intervals + posterior median for NON-units components only (common ancestry + hybridization)
# Uses larger points and wider error bars for readability
ggplot(OS_multinomial_results_filtered, aes(x = Variable)) +
  geom_point(aes(y = Median), color = "red", shape = 18, size = 6) +  # Posterior medians
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5, color = "black") +  # 95% HPD intervals
  labs(x = " ", y = "Value", 
       title = "HPD Intervals by Disease") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Rotate and enlarge x labels
    axis.text.y = element_text(size = 12),                         # Enlarge y labels
    axis.title.x = element_text(size = 14),                        # X-axis title size
    axis.title.y = element_text(size = 14),                        # Y-axis title size
    plot.title = element_text(size = 16, face = "bold"),           # Plot title styling
    plot.margin = margin(10, 10, 10, 20)                           # Add extra plot margin space
  )

# Same plot as above, but with axes flipped to make labels easier to read
# After coord_flip(): y = component name; x = median/HPD values
ggplot(OS_multinomial_results_filtered, aes(x = Variable, y = Median)) +
  geom_point(color = "red", shape = 18, size = 5) +  # Posterior medians
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5, color = "black") +  # 95% HPD intervals
  labs(x = " ", y = "HPD interval", 
       title = " ") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),   # Keep x labels horizontal
    axis.text.y = element_text(size = 12),                         # Enlarge y labels
    axis.title.x = element_text(size = 14),                        # X-axis title size
    axis.title.y = element_text(size = 14),                        # Y-axis title size
    plot.title = element_text(size = 16, face = "bold"),           # Plot title styling
    plot.margin = margin(10, 10, 10, 20)                           # Add extra plot margin space
  ) +
  coord_flip()  # Flip axes (categories on y, values on x)

# Compute posterior proportions for ALL OS variance components at each MCMC iteration
# Result: a data.frame with same number of columns as OS_multinomial$VCV, each column is a proportion
VCV_prop_OS <- as.data.frame(OS_multinomial$VCV / rowSums(OS_multinomial$VCV))

# Assign human-readable column names to the proportion table (must match VCV column order)
colnames(VCV_prop_OS) <- c(
  "Cardiovascular (common ancestry)", "Endocrine (common ancestry)", "Gastrointestinal (common ancestry)",
  "Hematopoietic (common ancestry)", "Musculoskeletal (common ancestry)", "Neurologic (common ancestry)",
  "Respiratory (common ancestry)", "Urogenital (common ancestry)",
  "Cardiovascular (hybridization)", "Endocrine (hybridization)", "Gastrointestinal (hybridization)",
  "Hematopoietic (hybridization)", "Musculoskeletal (hybridization)", "Neurologic (hybridization)",
  "Respiratory (hybridization)", "Urogenital (hybridization)",
  "Cardiovascular (units)", "Endocrine (units)", "Gastrointestinal (units)", "Hematopoietic (units)",
  "Musculoskeletal (units)", "Neurologic (units)", "Respiratory (units)", "Urogenital (units)"
)

# Add an explicit iteration index so we can reshape to long format safely
VCV_prop_OS$Iteration <- 1:nrow(VCV_prop_OS)

# Reshape posterior proportions into long format:
# one row per (Iteration × Component), with Proportion as the value
VCV_long_OS <- melt(VCV_prop_OS, id.vars = "Iteration", variable.name = "Component", value.name = "Proportion")

# Drop residual "units" components to keep only the two genetic effects
VCV_long_OS <- subset(VCV_long_OS, !grepl("units", VCV_long_OS$Component))

# Define a manual display order so each OS category appears as:
# common ancestry first, then hybridization, and so on across categories
ordered_components <- c(
  "Cardiovascular (common ancestry)", "Cardiovascular (hybridization)",
  "Endocrine (common ancestry)", "Endocrine (hybridization)",
  "Gastrointestinal (common ancestry)", "Gastrointestinal (hybridization)",
  "Hematopoietic (common ancestry)", "Hematopoietic (hybridization)",
  "Musculoskeletal (common ancestry)", "Musculoskeletal (hybridization)",
  "Neurologic (common ancestry)", "Neurologic (hybridization)",
  "Respiratory (common ancestry)", "Respiratory (hybridization)",
  "Urogenital (common ancestry)", "Urogenital (hybridization)"
)

# Apply the manual ordering to the Component factor (controls facet/label order)
VCV_long_OS$Component <- factor(VCV_long_OS$Component, levels = ordered_components)

# Compute the 95% HPD interval for each component’s posterior proportion distribution
hpd_df <- VCV_long_OS %>%
  group_by(Component) %>%
  summarise(
    Lower = HPDinterval(as.mcmc(Proportion))[1],
    Upper = HPDinterval(as.mcmc(Proportion))[2],
    .groups = "drop"
  )

# Compute the posterior median proportion for each component
median_df <- VCV_long_OS %>%
  group_by(Component) %>%
  summarise(
    Median = median(Proportion),
    .groups = "drop"
  )

# Plot posterior density distributions for each component (one facet per component)
# Adds dashed HPD bounds and a solid median line for each component
ggplot(VCV_long_OS, aes(x = Proportion)) +
  geom_density(fill = "lightblue", alpha = 0.7) +
  geom_vline(data = hpd_df, aes(xintercept = Lower), color = "black", linetype = "dashed") +
  geom_vline(data = hpd_df, aes(xintercept = Upper), color = "black", linetype = "dashed") +
  geom_vline(data = median_df, aes(xintercept = Median), color = "red", linetype = "solid") +
  facet_grid(Component ~ ., scales = "fixed", switch = "y") +
  scale_y_continuous(position = "right") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Posterior distributions and HPD intervals (OS)",
    x = "Proportion of variance explained", y = NULL
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 9, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 7),  # Reduce density-axis text size (visual cleanup)
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.spacing = unit(0.5, "lines")
  )

## ========= Helpers OS (uses objects already created: VCV_long_OS, hpd_df, median_df) =========

# Compute a single global maximum density value across ALL components,
# using a fixed x-range (from–to). This is used to force the same y-limit
# across separate per-category plots so they are directly comparable.
get_ymax_global_os <- function(df, from = 0, to = 0.35) {
  comps <- unique(df$Component)
  max(sapply(comps, function(cmp) {
    dvals <- df$Proportion[df$Component == cmp]
    if (length(na.omit(dvals)) < 2) return(0)
    dx <- density(dvals, from = from, to = to, na.rm = TRUE)
    max(dx$y, na.rm = TRUE)
  }), na.rm = TRUE)
}

# Calculate the global y-maximum for OS components (within x-range 0–0.35)
ymax_all_OS <- get_ymax_global_os(VCV_long_OS, from = 0, to = 0.35)

# Function to plot a single OS category only (two facets: common ancestry vs hybridization)
# with fixed x-limits and fixed y-limits for cross-category comparability.
plot_OS_category <- function(category,
                             x_limits = c(0, 0.35),
                             y_max = ymax_all_OS) {
  # Validate that 'category' is one of the expected OS base names
  valid <- c("Cardiovascular","Endocrine","Gastrointestinal","Hematopoietic",
             "Musculoskeletal","Neurologic","Respiratory","Urogenital")
  if (!category %in% valid) {
    stop(paste0("`category` debe ser uno de: ", paste(valid, collapse = ", ")))
  }
  
  # Build the exact component labels used in VCV_long_OS
  cat_common <- paste0(category, " (common ancestry)")
  cat_haplo  <- paste0(category, " (hybridization)")
  
  # Subset the long posterior draws, HPD table, and median table for these two components
  sub_df   <- dplyr::filter(VCV_long_OS, Component %in% c(cat_common, cat_haplo))
  sub_hpd  <- dplyr::filter(hpd_df,       Component %in% c(cat_common, cat_haplo))
  sub_med  <- dplyr::filter(median_df,    Component %in% c(cat_common, cat_haplo))
  
  # Force consistent facet order (common ancestry on top, hybridization below)
  sub_df$Component  <- factor(sub_df$Component,  levels = c(cat_common, cat_haplo))
  sub_hpd$Component <- factor(sub_hpd$Component, levels = c(cat_common, cat_haplo))
  sub_med$Component <- factor(sub_med$Component, levels = c(cat_common, cat_haplo))
  
  # Generate the per-category density plot with HPD bounds and median lines
  ggplot(sub_df, aes(x = Proportion)) +
    geom_density(fill = "lightblue", alpha = 0.7) +
    geom_vline(data = sub_hpd, aes(xintercept = Lower), color = "black", linetype = "dashed") +
    geom_vline(data = sub_hpd, aes(xintercept = Upper), color = "black", linetype = "dashed") +
    geom_vline(data = sub_med, aes(xintercept = Median), color = "red", linetype = "solid") +
    facet_grid(Component ~ ., scales = "fixed", switch = "y") +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(position = "right", limits = c(0, y_max)) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0("Posterior distributions and HPD intervals (OS) — ", category, " only"),
      x = "Proportion of variance explained", y = NULL
    ) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(size = 9, angle = 0, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.spacing = unit(0.5, "lines")
    )
}

## ========= Example calls (each returns a ggplot object) =========
plot_OS_category("Musculoskeletal")
plot_OS_category("Gastrointestinal")
plot_OS_category("Hematopoietic")

# Create a list to store the results (one small data.frame per VCV column/component)
PP_multinomial_results_list <- lapply(1:ncol(PP_multinomial$VCV), function(i) {
  # Compute the posterior proportion for component i at each MCMC iteration:
  # VCV[, i] / sum(VCV row) = proportion of total variance explained by that component
  PP_multinomial_variable_proportion <- PP_multinomial$VCV[, i] / rowSums(PP_multinomial$VCV)
  
  # Compute the 95% HPD interval of the posterior proportions for this component
  PP_multinomial_hpd <- HPDinterval(PP_multinomial_variable_proportion)
  
  # Compute the posterior median of the proportions for this component
  PP_multinomial_median_val <- median(PP_multinomial_variable_proportion)
  
  # Store results for this component in a one-row data.frame:
  # - Variable name (current VCV column name)
  # - Lower/Upper HPD bounds
  # - Posterior median
  data.frame(
    Variable = colnames(PP_multinomial$VCV)[i],
    Lower = PP_multinomial_hpd[1, "lower"],
    Upper = PP_multinomial_hpd[1, "upper"],
    Median = PP_multinomial_median_val
  )
})

# Combine the list of one-row data.frames into a single results table
PP_multinomial_results <- do.call(rbind, PP_multinomial_results_list)

# Print the results table (each row = one variance component; columns = HPD + median)
print(PP_multinomial_results)

# Define the desired human-readable labels for each variance component, in the exact order
# matching the columns of PP_multinomial$VCV (common ancestry, hybridization, then units)
correct_variable_names <- c(
  "Congenital (common ancestry)", "Degenerative (common ancestry)", "Infectious (common ancestry)", "Inflammatory (common ancestry)", "Metabolic (common ancestry)",
  "Neoplastic (common ancestry)", "Traumatic (common ancestry)", 
  "Congenital (hybridization)", "Degenerative (hybridization)", "Infectious (hybridization)", "Inflammatory (hybridization)", "Metabolic (hybridization)",
  "Neoplastic (hybridization)", "Traumatic (hybridization)", 
  "Congenital (units)", "Degenerative (units)", "Infectious (units)", "Inflammatory (units)", "Metabolic (units)",
  "Neoplastic (units)", "Traumatic (units)"
)

# Replace the "Variable" column in the results table with the human-readable names above
# (This assumes length(correct_variable_names) == nrow(PP_multinomial_results))
PP_multinomial_results$Variable <- correct_variable_names

# Print the relabeled results table to confirm the new names were applied correctly
print(PP_multinomial_results)

# Optional: export the final table to CSV for use in Excel/manuscript tables
# write.csv(PP_multinomial_results, "PP_multinomial_results_final.csv", row.names = FALSE)

# Create a filtered version of the results table excluding residual ("units") components,
# keeping only the two genetic effects (common ancestry and hybridization)
PP_multinomial_results_filtered <- PP_multinomial_results[!grepl("units", PP_multinomial_results$Variable), ]

# Print the filtered table to confirm only non-units components remain
print(PP_multinomial_results_filtered)

# This first block makes a simple point + HPD error bar plot for ALL components
# contained in PP_multinomial_results (assumed to already have Median, Lower, Upper, etc.).
# The x-axis shows each variance component (as a categorical label),
# the point is the posterior median, and the error bar is the 95% HPD interval.

# If PP_multinomial_results contains the columns 'Variable', 'Mean', 'Median', 'Lower', and 'Upper'
ggplot(PP_multinomial_results, aes(x = Variable)) +
  # Plot the posterior median as a red diamond-like point (shape 18)
  geom_point(aes(y = Median), color = "red", shape = 18, size = 5) +
  # Add vertical error bars representing the HPD interval [Lower, Upper]
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "black") +
  # Add axis labels and a title
  labs(x = "Variable", y = "Value", 
       title = "HPD Intervals with Mean and Median") +
  # Use a minimal theme
  theme_minimal() +
  # Rotate x-axis labels to avoid overlap (many variables/components)
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Reorder variables in reverse alphabetical order so the plot shows them
# in a controlled (reproducible) order rather than default factor ordering.
# Ordenar las variables en orden alfabético invertido
PP_multinomial_results_filtered$Variable <- factor(PP_multinomial_results_filtered$Variable, 
                                                   levels = rev(sort(unique(PP_multinomial_results_filtered$Variable))))

# Create the same median + HPD plot but now using the filtered results data.frame,
# which typically excludes "units" components and keeps only genetic components.
# Crear la figura
ggplot(PP_multinomial_results_filtered, aes(x = Variable)) +
  # Plot posterior median points (larger size than the first plot)
  geom_point(aes(y = Median), color = "red", shape = 18, size = 6) +  # Puntos de la mediana
  # Add HPD intervals with a wider bar width for readability
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5, color = "black") +  # Intervalos HPD
  # Set axis labels and title (blank x label used to visually de-emphasize it)
  labs(x = " ", y = "Value", 
       title = "HPD Intervals by Disease") +  # Títulos de los ejes y la gráfica
  theme_minimal() +
  # Customize text sizes and margins for publication-style formatting
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),  # Increase x-axis label size + rotate
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Title styling
    plot.margin = margin(10, 10, 10, 20)  # Add extra margin space (esp. left)
  )

# Create the same plot again but flipping coordinates so variables are on the y-axis
# (often easier to read long labels) and the HPD intervals run horizontally.
# Crear la figura con ejes invertidos
ggplot(PP_multinomial_results_filtered, aes(x = Variable, y = Median)) +
  # Plot posterior medians
  geom_point(color = "red", shape = 18, size = 5) +  # Puntos de la mediana
  # Add HPD intervals [Lower, Upper]
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5, color = "black") +  # Intervalos HPD
  # Provide axis labels; title intentionally left blank
  labs(x = " ", y = "HPD interval", 
       title = " ") +  # Títulos de los ejes y la gráfica
  theme_minimal() +
  # Formatting for readable labels after coordinate flip
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 12),  # Keep x-axis text unrotated
    axis.text.y = element_text(size = 12),  # Increase y-axis label size
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    plot.title = element_text(size = 16, face = "bold"),  # Title styling
    plot.margin = margin(10, 10, 10, 20)  # Add extra margin space
  ) +
  coord_flip()  # Flip axes so Variable labels are on the y-axis

# Compute posterior proportions of variance components by iteration:
# each VCV column divided by the row-wise sum gives the proportion of total variance
# explained by each component at that MCMC iteration.
# Calcular proporciones del modelo
VCV_prop_PP <- as.data.frame(PP_multinomial$VCV / rowSums(PP_multinomial$VCV))

# Assign human-readable names to each variance component column.
# These labels match the interpretation:
# - common ancestry = SNP similarity random effect
# - hybridization = haplotype-sharing random effect
# - units = residual variance component
# Asignar nombres legibles a las columnas (homologado con "hybridization")
colnames(VCV_prop_PP) <- c(
  "Congenital (common ancestry)", "Degenerative (common ancestry)", "Infectious (common ancestry)",
  "Inflammatory (common ancestry)", "Metabolic (common ancestry)", "Neoplastic (common ancestry)",
  "Traumatic (common ancestry)", 
  "Congenital (hybridization)", "Degenerative (hybridization)",
  "Infectious (hybridization)", "Inflammatory (hybridization)", "Metabolic (hybridization)",
  "Neoplastic (hybridization)", "Traumatic (hybridization)", 
  "Congenital (units)", "Degenerative (units)", "Infectious (units)", "Inflammatory (units)",
  "Metabolic (units)", "Neoplastic (units)", "Traumatic (units)"
)

# Add an iteration index so we can reshape to long format while keeping track
# of the MCMC iteration number.
# Añadir columna de iteración
VCV_prop_PP$Iteration <- 1:nrow(VCV_prop_PP)

# Convert the wide matrix of proportions into long format:
# each row becomes one (Iteration, Component, Proportion) observation,
# which is convenient for ggplot + grouped summaries.
# Convertir a formato largo
VCV_long_PP <- melt(VCV_prop_PP, id.vars = "Iteration", variable.name = "Component", value.name = "Proportion")

# Remove residual "units" components, keeping only the two genetic components
# (common ancestry and hybridization) for plotting and summarizing.
# Filtrar solo common ancestry y gene flow
VCV_long_PP <- subset(VCV_long_PP, !grepl("units", VCV_long_PP$Component))

# Define an explicit ordering for the facets:
# for each disease category, show common ancestry first, then hybridization,
# and repeat in the desired disease order.
# Definir orden manual deseado
ordered_components_PP <- c(
  "Congenital (common ancestry)", "Congenital (hybridization)",
  "Degenerative (common ancestry)", "Degenerative (hybridization)",
  "Infectious (common ancestry)", "Infectious (hybridization)",
  "Inflammatory (common ancestry)", "Inflammatory (hybridization)",
  "Metabolic (common ancestry)", "Metabolic (hybridization)",
  "Neoplastic (common ancestry)", "Neoplastic (hybridization)",
  "Traumatic (common ancestry)", "Traumatic (hybridization)"
)

# Apply the facet ordering by turning Component into an ordered factor.
VCV_long_PP$Component <- factor(VCV_long_PP$Component, levels = ordered_components_PP)

# Compute 95% HPD intervals for the posterior distribution of Proportion
# within each Component (densities are based on all MCMC iterations).
# Calcular HPD
hpd_df_PP <- VCV_long_PP %>%
  group_by(Component) %>%
  summarise(
    Lower = HPDinterval(as.mcmc(Proportion))[1],
    Upper = HPDinterval(as.mcmc(Proportion))[2],
    .groups = "drop"
  )

# Compute posterior medians for each Component, to be plotted as a solid vertical line.
# Calcular mediana
median_df_PP <- VCV_long_PP %>%
  group_by(Component) %>%
  summarise(
    Median = median(Proportion),
    .groups = "drop"
  )

# Plot posterior density distributions of the variance proportions:
# - one facet (row) per Component
# - dashed vertical lines = HPD bounds
# - solid vertical line = median
# - fixed x-axis limits to standardize visual comparison across facets
# Generar figura
ggplot(VCV_long_PP, aes(x = Proportion)) +
  # Density of posterior samples for each Component
  geom_density(fill = "lightyellow", alpha = 0.7) +
  # Lower/upper HPD boundaries (dashed)
  geom_vline(data = hpd_df_PP, aes(xintercept = Lower), color = "skyblue4", linetype = "dashed") +
  geom_vline(data = hpd_df_PP, aes(xintercept = Upper), color = "skyblue4", linetype = "dashed") +
  # Posterior median (solid red)
  geom_vline(data = median_df_PP, aes(xintercept = Median), color = "red", linetype = "solid") +
  # One row per Component; labels on the left; densities drawn in a single column
  facet_grid(Component ~ ., scales = "fixed", switch = "y") +
  # Fix x-axis range for comparability across all facets (PP)
  scale_x_continuous(limits = c(0, 0.35)) +  # ← here the x-range is explicitly constrained
  # Put y-axis ticks on the right (matches your vertical density layout style)
  scale_y_continuous(position = "right") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Posterior distributions and HPD intervals (PP)",
    x = "Proportion of variance explained", y = NULL
  ) +
  # Styling: move facet labels outside and tune text sizes/spacings
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(size = 9, angle = 0, hjust = 1),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.spacing = unit(0.5, "lines")
  )

# ------------------------------------------------------------
# Utility 1: compute a single global maximum density value across components,
# using the same x-range you plot (0 to 0.35). This is used to enforce the same
# y-axis (density) scale across multiple separate plots.
# --- 1) global y-maximum using the same X-range (0–0.35) ---
get_ymax_global_pp <- function(df, from = 0, to = 0.35) {
  comps <- unique(df$Component)
  max(sapply(comps, function(cmp) {
    dvals <- df$Proportion[df$Component == cmp]
    if (length(na.omit(dvals)) < 2) return(0)
    dx <- density(dvals, from = from, to = to, na.rm = TRUE)
    max(dx$y, na.rm = TRUE)
  }), na.rm = TRUE)
}

# Compute the global y-maximum once, to reuse in per-category plots
ymax_all_PP <- get_ymax_global_pp(VCV_long_PP, from = 0, to = 0.35)

# ------------------------------------------------------------
# Utility 2: function to plot ONE PP disease category at a time,
# showing two facets (common ancestry and hybridization) with:
# - density curves
# - dashed HPD lines
# - solid median line
# and with fixed x and fixed y scales for direct comparability across categories.
# --- 2) function to plot a specific PP category ---
plot_PP_category <- function(category,
                             x_limits = c(0, 0.35),
                             y_max = ymax_all_PP) {
  # Build exact component names expected in VCV_long_PP
  cat_common <- paste0(category, " (common ancestry)")
  cat_haplo  <- paste0(category, " (hybridization)")
  
  # Filter posterior samples and summary stats (HPD + median) for just this category
  sub_df   <- dplyr::filter(VCV_long_PP, Component %in% c(cat_common, cat_haplo))
  sub_hpd  <- dplyr::filter(hpd_df_PP,   Component %in% c(cat_common, cat_haplo))
  sub_med  <- dplyr::filter(median_df_PP,Component %in% c(cat_common, cat_haplo))
  
  # Enforce consistent facet order: common ancestry on top, hybridization below
  sub_df$Component <- factor(sub_df$Component, levels = c(cat_common, cat_haplo))
  sub_hpd$Component <- factor(sub_hpd$Component, levels = c(cat_common, cat_haplo))
  sub_med$Component <- factor(sub_med$Component, levels = c(cat_common, cat_haplo))
  
  # Produce the category-specific density plot with fixed x and fixed y
  ggplot(sub_df, aes(x = Proportion)) +
    geom_density(fill = "lightyellow", alpha = 0.7) +
    geom_vline(data = sub_hpd, aes(xintercept = Lower), color = "skyblue4", linetype = "dashed") +
    geom_vline(data = sub_hpd, aes(xintercept = Upper), color = "skyblue4", linetype = "dashed") +
    geom_vline(data = sub_med, aes(xintercept = Median), color = "red", linetype = "solid") +
    facet_grid(Component ~ ., scales = "fixed", switch = "y") +
    scale_x_continuous(limits = x_limits) +
    # Fix y-limit using the precomputed global max density
    scale_y_continuous(position = "right", limits = c(0, y_max)) +
    theme_minimal(base_size = 11) +
    labs(
      title = paste0("Posterior distributions and HPD intervals (PP) — ", category, " only"),
      x = "Proportion of variance explained", y = NULL
    ) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(size = 12, angle = 0, hjust = 1),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.spacing = unit(0.5, "lines")
    )
}

# --- 3) Example calls: generate category-specific PP plots ---
plot_PP_category("Neoplastic")
plot_PP_category("Congenital")
plot_PP_category("Metabolic")

# ------------------------------------------------------------
# Function: overlay the two genetic effects (common ancestry vs hybridization)
# for ONE category, in ONE panel, with fixed axis limits:
# - OS uses xlim from OS global proportion range
# - PP uses xlim = c(0, 0.35) to match your PP figure
# It also computes a global maximum density (across all components) so the y scale
# stays consistent with the facet version (scales = "fixed").
# Notes in the code mention needed objects and packages.
# Necesita: dplyr, ggplot2, y que existan
# VCV_long_OS, hpd_df, median_df  (OS model)
# VCV_long_PP, hpd_df_PP, median_df_PP (PP model)

plot_dual_effects <- function(model = c("OS", "PP"),
                              category = "Musculoskeletal",
                              xlim_override = NULL,
                              ylim_override = NULL) {
  model <- match.arg(model)
  
  # Select the correct long-format posterior object and summary tables (HPD + median)
  # depending on whether we are plotting OS or PP results.
  if (model == "OS") {
    df_long  <- VCV_long_OS
    hpd_all  <- hpd_df
    med_all  <- median_df
    # Global x-limits based on OS posterior proportions (for consistent scaling)
    xlim_global <- range(df_long$Proportion, na.rm = TRUE)
  } else {
    df_long  <- VCV_long_PP
    hpd_all  <- hpd_df_PP
    med_all  <- median_df_PP
    # PP x-limits fixed to match the earlier PP figure
    xlim_global <- c(0, 0.35)
  }
  
  # Build the exact component labels for this category
  comp_common <- paste0(category, " (common ancestry)")
  comp_hybrid <- paste0(category, " (hybridization)")
  
  # Keep only the two components to overlay
  df_two <- df_long %>%
    filter(Component %in% c(comp_common, comp_hybrid))
  
  # If the category name does not match expected labels, stop with an informative message
  if (nrow(df_two) == 0) {
    stop("No se encontraron datos para la categoría solicitada. ",
         "Asegúrate de usar el nombre base exactamente (p. ej., ",
         "'Cardiovascular', 'Endocrine', 'Gastrointestinal', ",
         "'Hematopoietic', 'Musculoskeletal', 'Neurologic', ",
         "'Respiratory', 'Urogenital', 'Congenital', 'Degenerative', ",
         "'Infectious', 'Inflammatory', 'Metabolic', 'Neoplastic', 'Traumatic').")
  }
  
  # HPD bounds and medians restricted to the two selected components
  hpd_two <- hpd_all %>% filter(Component %in% c(comp_common, comp_hybrid))
  med_two <- med_all %>% filter(Component %in% c(comp_common, comp_hybrid))
  
  # Compute a global maximum density across ALL components in the model,
  # to keep y-axis scaling comparable to the facet plot with fixed scales.
  dens_max <- df_long %>%
    split(.$Component) %>%
    lapply(function(dd) {
      x <- dd$Proportion
      x <- x[is.finite(x)]
      if (length(unique(x)) > 1) max(stats::density(x)$y) else 0
    }) %>%
    unlist() %>%
    max(na.rm = TRUE)
  
  # Allow optional manual overrides for x and y limits
  xlim_use <- if (is.null(xlim_override)) xlim_global else xlim_override
  ylim_use <- if (is.null(ylim_override)) c(0, dens_max) else ylim_override
  
  # Define colors for the two effects (fill + outline)
  pal <- c(
    setNames(if (model == "OS") "#2D8079" else "#1C404C", comp_common),   # common ancestry
    setNames(if (model == "OS") "#DD9D46" else "#CE5042", comp_hybrid)    # hybridization
  )
  
  # Overlay plot:
  # - two density curves (one per effect)
  # - dashed HPD bounds
  # - solid median
  ggplot(df_two, aes(x = Proportion, fill = Component, color = Component)) +
    geom_density(alpha = 0.72, adjust = 1) +
    # Lower and upper HPD bounds
    geom_vline(data = hpd_two, aes(xintercept = Lower, color = Component),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    geom_vline(data = hpd_two, aes(xintercept = Upper, color = Component),
               linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    # Median
    geom_vline(data = med_two, aes(xintercept = Median, color = Component),
               linetype = "solid", linewidth = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    coord_cartesian(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    labs(
      title = paste0("Posterior distributions and HPD (", model, " — ", category, ")"),
      x = "Proportion of variance explained",
      y = "Density",
      fill = "Effect", color = "Effect"
    ) +
    theme_minimal(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top"
    )
}

# Example calls: generate overlay plots for specific OS categories
# Ejemplos de uso:
# 1) OS, Musculoskeletal
plot_dual_effects(model = "OS", category = "Musculoskeletal")

# 2) OS, Gastrointestinal
plot_dual_effects(model = "OS", category = "Gastrointestinal")

# 3) OS, Hematopoietic
plot_dual_effects(model = "OS", category = "Hematopoietic")

# Example calls: generate overlay plots for specific PP categories (fixed xlim = c(0, 0.35))
# 4) PP, Neoplastic (respecting xlim = c(0, 0.35) from your PP plot)
plot_dual_effects(model = "PP", category = "Neoplastic")

# 5) PP, Congenital (respecting xlim = c(0, 0.35) from your PP plot)
plot_dual_effects(model = "PP", category = "Congenital")

# 6) PP, Metabolic (respecting xlim = c(0, 0.35) from your PP plot)
plot_dual_effects(model = "PP", category = "Metabolic")

# Create an output directory (if it does not already exist) for saving overlay plots
dir.create("figs_dual", showWarnings = FALSE, recursive = TRUE)

# Wrapper function to generate an overlay plot (via plot_dual_effects)
# and save it as a PNG with standardized dimensions and resolution.
save_dual_plot_simple <- function(model, category,
                                  out_dir = "figs_dual",
                                  width_in = 10, height_in = 8,
                                  dpi = 300) {
  p <- plot_dual_effects(model = model, category = category)
  
  # Utility to create file-safe names (replace non-alphanumerics with underscores)
  slug <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  base <- paste0(slug(model), "_", slug(category))
  
  # Save PNG at 300 dpi with a white background
  ggsave(file.path(out_dir, paste0(base, ".png")),
         plot = p, width = width_in, height = height_in,
         units = "in", dpi = dpi, bg = "white", device = "png")
  
  # Optional formats (commented out): TIFF (LZW) and PDF vector output
  # # TIFF 300 DPI (with LZW compression if available)
  # ggsave(file.path(out_dir, paste0(base, ".tiff")),
  #        plot = p, width = width_in, height = height_in,
  #        units = "in", dpi = dpi, bg = "white", device = "tiff",
  #        compression = "lzw")
  # 
  # # PDF vector output (dpi not relevant)
  # ggsave(file.path(out_dir, paste0(base, ".pdf")),
  #        plot = p, width = width_in, height = height_in,
  #        units = "in", device = "pdf", bg = "white")
}

# Print current working directory so you can confirm where outputs will be written
getwd()

###############################################################################
###############################################################################
## Posterior proportions with overlaid effects
## (common ancestry vs hybridization) by category, with fixed scales
## ==========

# Utility function: compute the global maximum density across components,
# used to fix the y-axis limits so all facets share the same density scale.
# Utilidad: calcula el máximo de densidad global (para fijar el mismo ylim)
dens_max_global <- function(df_long) {
  df_long %>%
    split(.$Component) %>%
    lapply(function(dd) {
      x <- dd$Proportion
      x <- x[is.finite(x)]
      if (length(unique(x)) > 1) max(stats::density(x)$y) else 0
    }) %>%
    unlist() %>%
    max(na.rm = TRUE)
}

# Generic constructor for a facet plot with overlaid densities (two effects per category),
# using fixed x and fixed y scales across facets.
# Constructor genérico para la figura facetada con traslape
plot_dual_facet <- function(model = c("OS", "PP"),
                            xlim_override = NULL,
                            ylim_override = NULL) {
  model <- match.arg(model)
  
  # Choose the correct long-format posterior object and summary tables
  # depending on OS vs PP.
  if (model == "OS") {
    df_long  <- VCV_long_OS                # must exist (generated earlier in your pipeline)
    hpd_all  <- hpd_df
    med_all  <- median_df
    # OS: fixed scales, x-range taken from the global OS posterior proportions
    xlim_global <- range(df_long$Proportion, na.rm = TRUE)
  } else {
    df_long  <- VCV_long_PP                # analogous object for PP
    hpd_all  <- hpd_df_PP
    med_all  <- median_df_PP
    # PP: fixed x-range to match your earlier PP plot
    xlim_global <- c(0, 0.35)
  }
  
  # Remove residual variance components so the figure only shows genetic effects
  # (common ancestry vs hybridization).
  df_long <- df_long %>% filter(!str_detect(Component, "\\(units\\)"))
  
  # Parse component labels to extract:
  # - Base: the disease category name
  # - Effect: which genetic component it corresponds to
  df_long <- df_long %>%
    mutate(
      Base   = str_replace(Component, " \\(.*\\)$", ""),
      Effect = case_when(
        str_detect(Component, "common ancestry") ~ "Common ancestry",
        str_detect(Component, "hybridization")   ~ "Hybridization",
        TRUE                                     ~ NA_character_
      )
    ) %>%
    filter(!is.na(Effect))
  
  # Set a controlled facet order for categories so the panel reads consistently
  if (model == "OS") {
    base_order <- c("Cardiovascular","Endocrine","Gastrointestinal","Hematopoietic",
                    "Musculoskeletal","Neurologic","Respiratory","Urogenital")
  } else {
    base_order <- c("Congenital","Degenerative","Infectious","Inflammatory",
                    "Metabolic","Neoplastic","Traumatic")
  }
  df_long$Base <- factor(df_long$Base, levels = base_order)
  
  # Convert HPD and median tables (which are indexed by Component) into
  # Base × Effect summaries so we can draw correct vertical lines per facet.
  hpd_use <- hpd_all %>%
    mutate(
      Base   = str_replace(Component, " \\(.*\\)$", ""),
      Effect = case_when(
        str_detect(Component, "common ancestry") ~ "Common ancestry",
        str_detect(Component, "hybridization")   ~ "Hybridization",
        TRUE                                     ~ NA_character_
      )
    ) %>%
    filter(!is.na(Effect)) %>%
    select(Base, Effect, Lower, Upper)
  
  med_use <- med_all %>%
    mutate(
      Base   = str_replace(Component, " \\(.*\\)$", ""),
      Effect = case_when(
        str_detect(Component, "common ancestry") ~ "Common ancestry",
        str_detect(Component, "hybridization")   ~ "Hybridization",
        TRUE                                     ~ NA_character_
      )
    ) %>%
    filter(!is.na(Effect)) %>%
    select(Base, Effect, Median)
  
  # Compute global y-limit (maximum density) so every facet shares the same y scale
  dens_max <- dens_max_global(df_long %>%
                                mutate(Component = paste0(Base, " (", Effect, ")")))
  xlim_use <- if (is.null(xlim_override)) xlim_global else xlim_override
  ylim_use <- if (is.null(ylim_override)) c(0, dens_max) else ylim_override
  
  # Two-color palette for the two effects, using transparency to visualize overlap
  pal_fill <- c("Common ancestry" = "#2D8079", "Hybridization" = "gray")
  pal_col  <- c("Common ancestry" = "#2D8079", "Hybridization" = "gray")
  
  # Build the facet plot:
  # - within each Base (category), overlay density curves for the two effects
  # - add HPD bounds (dashed) and median (solid)
  ggplot(df_long, aes(x = Proportion, fill = Effect, color = Effect)) +
    geom_density(alpha = 0.85, adjust = 1) +
    # HPD bounds per effect
    geom_vline(data = hpd_use, aes(xintercept = Lower, color = Effect),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    geom_vline(data = hpd_use, aes(xintercept = Upper, color = Effect),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    # Medians per effect
    geom_vline(data = med_use, aes(xintercept = Median, color = Effect),
               linetype = "solid", linewidth = 0.6, show.legend = FALSE) +
    facet_grid(Base ~ ., scales = "fixed", switch = "y") +
    scale_fill_manual(values = pal_fill) +
    scale_color_manual(values = pal_col) +
    coord_cartesian(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    scale_y_continuous(position = "right") +
    labs(
      title = paste0("Posterior distributions and HPD (", model, ")"),
      x = "Proportion of variance explained",
      y = NULL,
      fill = "Effect", color = "Effect"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(size = 9, angle = 0, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.spacing = unit(0.5, "lines"),
      legend.position = "top"
    )
}

## Example usage:
# Build and print the OS facet overlay figure
# Ejemplos:
# Figure 1 (OS) with overlap:
fig1_OS <- plot_dual_facet(model = "OS")
print(fig1_OS)
# ggsave("Figure1_OS_dualEffects.png", fig1_OS, width = 6, height = 8, dpi = 300)

# Build and print the PP facet overlay figure
# Figure 2 (PP) with overlap:
fig2_PP <- plot_dual_facet(model = "PP")
print(fig2_PP)
# ggsave("Figure2_PP_dualEffects.png", fig2_PP, width = 6, height = 7, dpi = 300)

################V EXTRA FIGURES

## ============================================================
## OS_multinomial_weight
## ============================================================

# This section builds a plotting table (posterior mean + 95% HPD)
# for the trait-specific effect of body weight on OS causes of death,
# then produces a horizontal forest-style plot.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the weight effect across OS causes of death
table <- data.frame(
  Cause    = c("Cardiovascular", "Endocrine", "Gastrointestinal",
               "Hematopoietic", "Musculoskeletal", "Neurologic",
               "Respiratory", "Urogenital"),
  
  Post.Mean = c(0.031, -0.225, 0.146, 0.111, 0.349, -0.061, -0.003, -0.117),
  l.95      = c(-0.057, -0.340, 0.076, 0.005, 0.271, -0.149, -0.087, -0.193),
  u.95      = c(0.121, -0.115, 0.208, 0.212, 0.428, 0.029, 0.086, -0.034)
)

# Order causes by the posterior mean effect size and label whether the 95% HPD
# interval excludes zero (interpreted as a credible non-zero effect)
table <- table %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Create a horizontal forest plot:
# - dashed line at 0 (null effect)
# - horizontal error bars = 95% HPD
# - points = posterior means, filled by whether HPD excludes 0
p1 <- ggplot(table5, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of body weight (posterior mean, 95% HPD)",
       y = "Cause of death (OS)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PP_multinomial_weight
## ============================================================
# Same plotting workflow as above, but for PP causes of death:
# build a table with posterior mean + 95% HPD for the weight effect,
# then produce a forest-style plot.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the weight effect across PP causes of death
table_PP_weight <- data.frame(
  Cause = c("Infectious", "Degenerative", "Congenital",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(0.101, -0.145, 0.090,
                -0.083, -0.021, 0.258,
                0.057),
  
  l.95 = c(0.009, -0.256, -0.040,
           -0.215, -0.109, 0.144,
           -0.037),
  
  u.95 = c(0.194, -0.027, 0.218,
           0.065, 0.062, 0.382,
           0.154)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_weight <- table_PP_weight %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Forest plot for PP: posterior mean points + HPD intervals; dashed line at 0
p2 <- ggplot(table_PP_weight, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of body weight (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## OS_multinomial_growth
## ============================================================
# Same workflow, but now plotting the trait-specific effect of growth
# on OS causes of death.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the growth effect across OS causes of death
table_OS_growth <- data.frame(
  Cause = c("Cardiovascular", "Endocrine", "Gastrointestinal",
            "Hematopoietic", "Musculoskeletal", "Neurologic",
            "Respiratory", "Urogenital"),
  
  Post.Mean = c(-0.029, -0.168, 0.089,
                0.132, 0.268, -0.095,
                -0.056, -0.089),
  
  l.95 = c(-0.135, -0.291, 0.008,
           0.021, 0.176, -0.201,
           -0.144, -0.187),
  
  u.95 = c(0.070, -0.040, 0.175,
           0.257, 0.367, 0.006,
           0.031, -0.009)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_OS_growth <- table_OS_growth %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Forest plot for OS growth effects
p3 <- ggplot(table_OS_growth, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of growth (posterior mean, 95% HPD)",
       y = "Cause of death (OS)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PP_multinomial_growth
## ============================================================
# Same workflow, but for the trait-specific effect of growth on PP causes of death.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the growth effect across PP causes of death
table_PP_growth <- data.frame(
  Cause = c("Infectious", "Congenital", "Degenerative",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(0.070, 0.071, -0.120,
                -0.023, -0.007, 0.271,
                0.040),
  
  l.95 = c(-0.030, -0.063, -0.251,
           -0.156, -0.104, 0.132,
           -0.071),
  
  u.95 = c(0.187, 0.207, 0.004,
           0.118, 0.092, 0.410,
           0.149)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_growth <- table_PP_growth %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Forest plot for PP growth effects
p4 <- ggplot(table_PP_growth, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of growth (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## 5. PANEL 2×2 — base grid only (no external packages)
## ============================================================
# Combine the four ggplot objects into a single 2x2 panel layout
# using the base 'grid' system.

library(grid)   # base R package used to arrange plots on a grid layout

# Initialize a new plotting page and create a 2x2 layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

# Helper function to place plots into specific grid cells
vplayout <- function(row, col) viewport(layout.pos.row = row,
                                        layout.pos.col = col)

# Print each plot into its assigned position
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
print(p3, vp = vplayout(2, 1))
print(p4, vp = vplayout(2, 2))


######################################
######################################

## -----------------------------------------
## Table: PP × reproductive investment
## -----------------------------------------
# Create a plotting table for the PP model including reproductive investment:
# posterior mean and 95% HPD interval for each PP cause.

table_PP_repinv <- data.frame(
  Cause = c("Congenital", "Degenerative", "Inflammatory",
            "Infectious", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.477, -0.427, -0.372,
                -0.052, -0.231, 0.040,
                -0.034),
  
  l.95 = c(-0.785, -0.749, -0.663,
           -0.249, -0.474, -0.067,
           -0.240),
  
  u.95 = c(-0.167, -0.107, -0.070,
           0.137, 0.011, 0.159,
           0.186)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_repinv <- table_PP_repinv %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: reproductive investment effect (PP)
p_PP_repinv <- ggplot(table_PP_repinv, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  
  labs(x = "Effect of reproductive investment (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Display the PP × reproductive investment plot
p_PP_repinv

## -----------------------------------------
## Table: PP × trainability
## -----------------------------------------
# Create a plotting table for the PP model including trainability:
# posterior mean and 95% HPD interval for each PP cause.

table_PP_train <- data.frame(
  Cause = c("Congenital", "Degenerative", "Infectious",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.151, 0.046, -0.116,
                0.091, -0.036, 0.031,
                -0.028),
  
  l.95 = c(-0.380, -0.149, -0.222,
           -0.055, -0.182, -0.033,
           -0.154),
  
  u.95 = c(0.093, 0.223, -0.002,
           0.269, 0.107, 0.094,
           0.092)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_train <- table_PP_train %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: trainability effect (PP)
p_PP_train <- ggplot(table_PP_train, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  
  labs(x = "Effect of trainability (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Display the PP × trainability plot
p_PP_train

## ============================================================
## PP × REPRODUCTIVE INVESTMENT
## ============================================================
# Repeat the PP × reproductive investment plotting workflow:
# define the table, flag HPD significance, and build a forest plot.

table_PP_repinv <- data.frame(
  Cause = c("Congenital", "Degenerative", "Inflammatory",
            "Infectious", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.477, -0.427, -0.372,
                -0.052, -0.231, 0.040,
                -0.034),
  
  l.95 = c(-0.785, -0.749, -0.663,
           -0.249, -0.474, -0.067,
           -0.240),
  
  u.95 = c(-0.167, -0.107, -0.070,
           0.137, 0.011, 0.159,
           0.186)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_repinv <- table_PP_repinv %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: reproductive investment effect (PP)
p_PP_repinv <- ggplot(table_PP_repinv, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of reproductive investment (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PP × TRAINABILITY
## ============================================================
# Repeat the PP × trainability plotting workflow:
# define the table, flag HPD significance, and build a forest plot.

table_PP_train <- data.frame(
  Cause = c("Congenital", "Degenerative", "Infectious",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.151, 0.046, -0.116,
                0.091, -0.036, 0.031,
                -0.028),
  
  l.95 = c(-0.380, -0.149, -0.222,
           -0.055, -0.182, -0.033,
           -0.154),
  
  u.95 = c(0.093, 0.223, -0.002,
           0.269, 0.107, 0.094,
           0.092)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_train <- table_PP_train %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: trainability effect (PP)
p_PP_train <- ggplot(table_PP_train, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of trainability (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PANEL 1×2 (base grid only; no extra packages)
## ============================================================
# Arrange the PP reproductive investment plot and PP trainability plot side-by-side
# in a 1x2 panel using the base grid system.

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

# Helper function to place plots into specific grid cells
vplayout <- function(row, col) viewport(layout.pos.row = row,
                                        layout.pos.col = col)

# Print each plot into its assigned position
print(p_PP_repinv, vp = vplayout(1, 1))
print(p_PP_train,  vp = vplayout(1, 2))


# ------------------------------------------------------------
# Diagnostics / descriptive plots: relationship between weight and growth
# ------------------------------------------------------------

# Correlation test between standardized weight and standardized growth
cor.test(data$scaleweight, data$scalegrowth)

# Scatter plot of standardized weight vs standardized growth with semi-transparent points
plot(data$scaleweight, data$scalegrowth,
     pch = 19, col = rgb(0,0,1,0.4),
     xlab = "Scale weight (adult size A)",
     ylab = "Scale growth (G)")

# Add a linear regression line of growth ~ weight
abline(lm(scalegrowth ~ scaleweight, data = data), col="red")

# Repeat scatter plot (same variables), using the same point styling
plot(data$scaleweight, data$scalegrowth,
     pch = 19,
     col = rgb(0,0,1,0.4),
     xlab = "Scale weight (adult size A)",
     ylab = "Scale growth (G)")

# Add the regression line with increased line width for visibility
abline(lm(scalegrowth ~ scaleweight, data = data), col = "red", lwd = 2)

# Add a text annotation with the correlation value on the plot
text(x = min(data$scaleweight, na.rm = TRUE),
     y = max(data$scalegrowth, na.rm = TRUE),
     labels = paste0("r = ", round(0.8867644, 3)),
     pos = 4, cex = 1.3, font = 2,
     col = "black")





























## ============================================================
## OS_multinomial_weight
## ============================================================
# This section builds a plotting table (posterior mean + 95% HPD)
# for the trait-specific effect of body weight on OS causes of death,
# then produces a horizontal forest-style plot.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the weight effect across OS causes of death
table <- data.frame(
  Cause    = c("Cardiovascular", "Endocrine", "Gastrointestinal",
               "Hematopoietic", "Musculoskeletal", "Neurologic",
               "Respiratory", "Urogenital"),
  
  Post.Mean = c(0.031, -0.225, 0.146, 0.111, 0.349, -0.061, -0.003, -0.117),
  l.95      = c(-0.057, -0.340, 0.076, 0.005, 0.271, -0.149, -0.087, -0.193),
  u.95      = c(0.121, -0.115, 0.208, 0.212, 0.428, 0.029, 0.086, -0.034)
)

# Order causes by the posterior mean effect size and label whether the 95% HPD
# interval excludes zero (interpreted as a credible non-zero effect)
table <- table %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Create a horizontal forest plot:
# - dashed line at 0 (null effect)
# - horizontal error bars = 95% HPD
# - points = posterior means, filled by whether HPD excludes 0
p1 <- ggplot(table5, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of body weight (posterior mean, 95% HPD)",
       y = "Cause of death (OS)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PP_multinomial_weight
## ============================================================
# Same plotting workflow as above, but for PP causes of death:
# build a table with posterior mean + 95% HPD for the weight effect,
# then produce a forest-style plot.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the weight effect across PP causes of death
table_PP_weight <- data.frame(
  Cause = c("Infectious", "Degenerative", "Congenital",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(0.101, -0.145, 0.090,
                -0.083, -0.021, 0.258,
                0.057),
  
  l.95 = c(0.009, -0.256, -0.040,
           -0.215, -0.109, 0.144,
           -0.037),
  
  u.95 = c(0.194, -0.027, 0.218,
           0.065, 0.062, 0.382,
           0.154)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_weight <- table_PP_weight %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Forest plot for PP: posterior mean points + HPD intervals; dashed line at 0
p2 <- ggplot(table_PP_weight, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of body weight (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## OS_multinomial_growth
## ============================================================
# Same workflow, but now plotting the trait-specific effect of growth
# on OS causes of death.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the growth effect across OS causes of death
table_OS_growth <- data.frame(
  Cause = c("Cardiovascular", "Endocrine", "Gastrointestinal",
            "Hematopoietic", "Musculoskeletal", "Neurologic",
            "Respiratory", "Urogenital"),
  
  Post.Mean = c(-0.029, -0.168, 0.089,
                0.132, 0.268, -0.095,
                -0.056, -0.089),
  
  l.95 = c(-0.135, -0.291, 0.008,
           0.021, 0.176, -0.201,
           -0.144, -0.187),
  
  u.95 = c(0.070, -0.040, 0.175,
           0.257, 0.367, 0.006,
           0.031, -0.009)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_OS_growth <- table_OS_growth %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Forest plot for OS growth effects
p3 <- ggplot(table_OS_growth, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of growth (posterior mean, 95% HPD)",
       y = "Cause of death (OS)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PP_multinomial_growth
## ============================================================
# Same workflow, but for the trait-specific effect of growth on PP causes of death.

# Manually assemble posterior summaries (posterior mean and 95% HPD interval)
# for the growth effect across PP causes of death
table_PP_growth <- data.frame(
  Cause = c("Infectious", "Congenital", "Degenerative",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(0.070, 0.071, -0.120,
                -0.023, -0.007, 0.271,
                0.040),
  
  l.95 = c(-0.030, -0.063, -0.251,
           -0.156, -0.104, 0.132,
           -0.071),
  
  u.95 = c(0.187, 0.207, 0.004,
           0.118, 0.092, 0.410,
           0.149)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_growth <- table_PP_growth %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig = ifelse(l.95 > 0 | u.95 < 0, "HPD excludes 0", "HPD includes 0")
  )

# Forest plot for PP growth effects
p4 <- ggplot(table_PP_growth, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig), height = 0.2) +
  geom_point(aes(fill = sig), shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of growth (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## 5. PANEL 2×2 — base grid only (no external packages)
## ============================================================
# Combine the four ggplot objects into a single 2x2 panel layout
# using the base 'grid' system.

library(grid)   # base R package used to arrange plots on a grid layout

# Initialize a new plotting page and create a 2x2 layout
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))

# Helper function to place plots into specific grid cells
vplayout <- function(row, col) viewport(layout.pos.row = row,
                                        layout.pos.col = col)

# Print each plot into its assigned position
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(1, 2))
print(p3, vp = vplayout(2, 1))
print(p4, vp = vplayout(2, 2))


######################################
######################################

## -----------------------------------------
## Table: PP × reproductive investment
## -----------------------------------------
# Create a plotting table for the PP model including reproductive investment:
# posterior mean and 95% HPD interval for each PP cause.

table_PP_repinv <- data.frame(
  Cause = c("Congenital", "Degenerative", "Inflammatory",
            "Infectious", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.477, -0.427, -0.372,
                -0.052, -0.231, 0.040,
                -0.034),
  
  l.95 = c(-0.785, -0.749, -0.663,
           -0.249, -0.474, -0.067,
           -0.240),
  
  u.95 = c(-0.167, -0.107, -0.070,
           0.137, 0.011, 0.159,
           0.186)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_repinv <- table_PP_repinv %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: reproductive investment effect (PP)
p_PP_repinv <- ggplot(table_PP_repinv, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  
  labs(x = "Effect of reproductive investment (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Display the PP × reproductive investment plot
p_PP_repinv

## -----------------------------------------
## Table: PP × trainability
## -----------------------------------------
# Create a plotting table for the PP model including trainability:
# posterior mean and 95% HPD interval for each PP cause.

table_PP_train <- data.frame(
  Cause = c("Congenital", "Degenerative", "Infectious",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.151, 0.046, -0.116,
                0.091, -0.036, 0.031,
                -0.028),
  
  l.95 = c(-0.380, -0.149, -0.222,
           -0.055, -0.182, -0.033,
           -0.154),
  
  u.95 = c(0.093, 0.223, -0.002,
           0.269, 0.107, 0.094,
           0.092)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_train <- table_PP_train %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: trainability effect (PP)
p_PP_train <- ggplot(table_PP_train, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  
  labs(x = "Effect of trainability (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

# Display the PP × trainability plot
p_PP_train

## ============================================================
## PP × REPRODUCTIVE INVESTMENT
## ============================================================
# Repeat the PP × reproductive investment plotting workflow:
# define the table, flag HPD significance, and build a forest plot.

table_PP_repinv <- data.frame(
  Cause = c("Congenital", "Degenerative", "Inflammatory",
            "Infectious", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.477, -0.427, -0.372,
                -0.052, -0.231, 0.040,
                -0.034),
  
  l.95 = c(-0.785, -0.749, -0.663,
           -0.249, -0.474, -0.067,
           -0.240),
  
  u.95 = c(-0.167, -0.107, -0.070,
           0.137, 0.011, 0.159,
           0.186)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_repinv <- table_PP_repinv %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: reproductive investment effect (PP)
p_PP_repinv <- ggplot(table_PP_repinv, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of reproductive investment (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PP × TRAINABILITY
## ============================================================
# Repeat the PP × trainability plotting workflow:
# define the table, flag HPD significance, and build a forest plot.

table_PP_train <- data.frame(
  Cause = c("Congenital", "Degenerative", "Infectious",
            "Inflammatory", "Metabolic", "Neoplastic",
            "Traumatic"),
  
  Post.Mean = c(-0.151, 0.046, -0.116,
                0.091, -0.036, 0.031,
                -0.028),
  
  l.95 = c(-0.380, -0.149, -0.222,
           -0.055, -0.182, -0.033,
           -0.154),
  
  u.95 = c(0.093, 0.223, -0.002,
           0.269, 0.107, 0.094,
           0.092)
)

# Order causes by posterior mean and flag whether HPD excludes zero
table_PP_train <- table_PP_train %>%
  mutate(
    Cause = factor(Cause, levels = Cause[order(Post.Mean)]),
    sig   = ifelse(l.95 > 0 | u.95 < 0,
                   "HPD excludes 0", "HPD includes 0")
  )

# Forest plot: trainability effect (PP)
p_PP_train <- ggplot(table_PP_train, aes(x = Post.Mean, y = Cause)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(xmin = l.95, xmax = u.95, color = sig),
                 height = 0.2) +
  geom_point(aes(fill = sig),
             shape = 21, size = 3, stroke = 0.4) +
  labs(x = "Effect of trainability (posterior mean, 95% HPD)",
       y = "Cause of death (PP)", color = "", fill = "") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")



## ============================================================
## PANEL 1×2 (base grid only; no extra packages)
## ============================================================
# Arrange the PP reproductive investment plot and PP trainability plot side-by-side
# in a 1x2 panel using the base grid system.

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2)))

# Helper function to place plots into specific grid cells
vplayout <- function(row, col) viewport(layout.pos.row = row,
                                        layout.pos.col = col)

# Print each plot into its assigned position
print(p_PP_repinv, vp = vplayout(1, 1))
print(p_PP_train,  vp = vplayout(1, 2))


# ------------------------------------------------------------
# Diagnostics / descriptive plots: relationship between weight and growth
# ------------------------------------------------------------

# Correlation test between standardized weight and standardized growth
cor.test(data$scaleweight, data$scalegrowth)

# Scatter plot of standardized weight vs standardized growth with semi-transparent points
plot(data$scaleweight, data$scalegrowth,
     pch = 19, col = rgb(0,0,1,0.4),
     xlab = "Scale weight (adult size A)",
     ylab = "Scale growth (G)")

# Add a linear regression line of growth ~ weight
abline(lm(scalegrowth ~ scaleweight, data = data), col="red")

# Repeat scatter plot (same variables), using the same point styling
plot(data$scaleweight, data$scalegrowth,
     pch = 19,
     col = rgb(0,0,1,0.4),
     xlab = "Scale weight (adult size A)",
     ylab = "Scale growth (G)")

# Add the regression line with increased line width for visibility
abline(lm(scalegrowth ~ scaleweight, data = data), col = "red", lwd = 2)

# Add a text annotation with the correlation value on the plot
text(x = min(data$scaleweight, na.rm = TRUE),
     y = max(data$scalegrowth, na.rm = TRUE),
     labels = paste0("r = ", round(0.8867644, 3)),
     pos = 4, cex = 1.3, font = 2,
     col = "black")

