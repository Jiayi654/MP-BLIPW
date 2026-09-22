##############################################################################
# 04_nhanes_analysis_hip_FIXED.R
#
# Final NHANES MP-BLIPW analysis using BMI + HIP CIRCUMFERENCE.
#
# Requires:
#   00_mp_blipw_glm.R
#   nhanes_mpblipw_data_hip.csv
#
# Methods:
#   1. Complete-case OLS
#   2. SIPW
#   3. BLIPW (BMI)
#   4. BLIPW (Hip)
#   5. MP-BLIPW (BMI + Hip)
#
# This version includes explicit rank checks so that coding problems are caught
# before solve(G1) is called.
##############################################################################

# Optional: uncomment if you want to force the local project folder.
# setwd("/Users/jyw/Downloads/mpblipw")


##############################################################################
# 1. Locate and source the core estimator
##############################################################################

core_file <- "00_mp_blipw_glm.R"

if (!file.exists(core_file)) {
  stop(
    "Cannot find ", core_file, " in the current working directory:\n",
    getwd(),
    "\n\nPut 00_mp_blipw_glm.R in the same folder as this script, or run:\n",
    'setwd("/Users/jyw/Downloads/mpblipw")'
  )
}

source(core_file)


##############################################################################
# 2. Load prepared NHANES data
##############################################################################

data_file <- "nhanes_mpblipw_data_hip.csv"

if (!file.exists(data_file)) {
  stop(
    "Cannot find ", data_file, ".\n",
    "Run 03_prepare_nhanes_hip_FIXED.R first."
  )
}

dat <- read.csv(
  data_file,
  na.strings = c("", "NA", "NaN"),
  stringsAsFactors = FALSE
)

required_vars <- c(
  "Y", "X", "W1", "W2", "Age",
  "Male", "NHBlack", "Hispanic", "Smoker", "Diabetes"
)

missing_vars <- setdiff(required_vars, names(dat))

if (length(missing_vars) > 0) {
  stop(
    "Prepared dataset is missing required variables: ",
    paste(missing_vars, collapse = ", ")
  )
}

# The preparation script already restricts to adults with observed HDL.
dat <- dat[!is.na(dat$Y), , drop = FALSE]

N <- nrow(dat)

if (N == 0) {
  stop("The analysis dataset contains zero rows.")
}


##############################################################################
# 3. Define complete-case and auxiliary-phase indicators
##############################################################################

# Complete case = DEXA body-fat measurement observed.
R <- as.integer(!is.na(dat$X))

# Each auxiliary phase must contain the CC sample.
# If a CC participant happens to have the auxiliary variable missing, phase
# membership is forced to 1, matching the original working analysis.
R1 <- pmax(R, as.integer(!is.na(dat$W1)))
R2 <- pmax(R, as.integer(!is.na(dat$W2)))

if (!all(R1[R == 1] == 1)) {
  stop("Phase 1 nesting condition C_CC subset C_1 is violated.")
}

if (!all(R2[R == 1] == 1)) {
  stop("Phase 2 nesting condition C_CC subset C_2 is violated.")
}


##############################################################################
# 4. Descriptive sample information
##############################################################################

cat("\n============================================================\n")
cat("NHANES MP-BLIPW: BMI + HIP CIRCUMFERENCE\n")
cat("============================================================\n")

cat(sprintf("N                          = %d\n", N))
cat(sprintf(
  "X (DEXA %% body fat) obs.  = %d (%.3f%%)\n",
  sum(R),
  100 * mean(R)
))
cat(sprintf(
  "Phase 1 (BMI or CC)        = %d (%.3f%%)\n",
  sum(R1),
  100 * mean(R1)
))
cat(sprintf(
  "Phase 2 (Hip or CC)        = %d (%.3f%%)\n",
  sum(R2),
  100 * mean(R2)
))


##############################################################################
# 5. Prepare model variables
##############################################################################

# Standardize age to improve numerical conditioning.
age_s <- as.numeric(scale(dat$Age))

if (anyNA(age_s) || any(!is.finite(age_s))) {
  stop("Standardized age contains missing or non-finite values.")
}

# X is unused for R=0 observations in the IPW estimating equations.
# Fill those placeholders with zero so the full-length design matrix is valid.
X_use <- dat$X
X_use[is.na(X_use)] <- 0

# W1/W2 values outside their observed phase do not contribute to that phase's
# estimating equations. Median values are used only as harmless placeholders.
W1_use <- dat$W1
W2_use <- dat$W2

if (all(is.na(W1_use))) stop("W1 (BMI) is completely missing.")
if (all(is.na(W2_use))) stop("W2 (hip circumference) is completely missing.")

W1_use[is.na(W1_use)] <- median(W1_use, na.rm = TRUE)
W2_use[is.na(W2_use)] <- median(W2_use, na.rm = TRUE)

# Primary incomplete covariate matrix.
X_mat <- matrix(X_use, ncol = 1)
colnames(X_mat) <- "DEXA_bodyfat"

# Fully observed adjustment covariates in the final primary model.
Z_mat <- cbind(
  Age_SD = age_s,
  Male = dat$Male,
  NHBlack = dat$NHBlack,
  Hispanic = dat$Hispanic,
  Smoker = dat$Smoker
)

if (anyNA(Z_mat) || any(!is.finite(Z_mat))) {
  stop(
    "Z_mat contains missing/non-finite values.\n",
    "Check sex, race, Hispanic, and smoking recoding in 03_prepare_nhanes_hip_FIXED.R."
  )
}

# Working models use intercept + surrogate, matching the working hip analysis.
D1 <- cbind(
  Intercept = 1,
  BMI = W1_use
)

D2 <- cbind(
  Intercept = 1,
  Hip = W2_use
)


##############################################################################
# 6. CHECK PRIMARY DESIGN MATRIX RANK
##############################################################################

ii <- which(R == 1)

if (length(ii) == 0) {
  stop("There are no complete cases with DEXA observed.")
}

D_check <- cbind(
  Intercept = 1,
  DEXA_bodyfat = X_mat[ii, 1],
  Z_mat[ii, , drop = FALSE]
)

cat("\nPrimary model complete-case diagnostics:\n")
cat("  Complete cases =", nrow(D_check), "\n")
cat("  Number of columns =", ncol(D_check), "\n")
cat("  Matrix rank =", qr(D_check)$rank, "\n")

unique_counts <- apply(
  D_check,
  2,
  function(x) length(unique(x))
)

cat("\nUnique values per primary-model column:\n")
print(unique_counts)

if (qr(D_check)$rank < ncol(D_check)) {
  cat("\nBinary covariate distributions among complete cases:\n")
  print(table(Male = dat$Male[ii], useNA = "ifany"))
  print(table(NHBlack = dat$NHBlack[ii], useNA = "ifany"))
  print(table(Hispanic = dat$Hispanic[ii], useNA = "ifany"))
  print(table(Smoker = dat$Smoker[ii], useNA = "ifany"))

  stop(
    "\nPrimary design matrix is rank deficient. ",
    "Do not add a ridge term just to bypass this error. ",
    "Fix the variable coding first."
  )
}


##############################################################################
# 7. Phase inclusion probabilities
##############################################################################

# Match the original working analysis:
# empirical phase inclusion probability for non-CC participants, while CC
# participants are included in every phase with probability 1.
pi1 <- rep(mean(R1), N)
pi2 <- rep(mean(R2), N)

pi1[R == 1] <- 1
pi2[R == 1] <- 1

pi1 <- pmax(pi1, 0.05)
pi2 <- pmax(pi2, 0.05)


##############################################################################
# 8. Complete-case selection model
##############################################################################

# Match the original working hip analysis:
# P(DEXA observed | standardized age, sex, race, diabetes)
#
# Diabetes is included in the selection model, but not the primary outcome
# model. Missing values were set to 0 in the preparation script.
D_sel <- cbind(
  Intercept = 1,
  Age_SD = age_s,
  Male = dat$Male,
  NHBlack = dat$NHBlack,
  Hispanic = dat$Hispanic,
  Diabetes = dat$Diabetes
)

if (anyNA(D_sel) || any(!is.finite(D_sel))) {
  stop("Selection-model design matrix contains missing/non-finite values.")
}

cat("\nSelection-model diagnostics:\n")
cat("  Number of columns =", ncol(D_sel), "\n")
cat("  Matrix rank =", qr(D_sel)$rank, "\n")

if (qr(D_sel)$rank < ncol(D_sel)) {
  stop(
    "Selection-model design matrix is rank deficient. ",
    "Check age/sex/race/diabetes coding before fitting selection probabilities."
  )
}

pi_obj <- fit_selection(
  R = R,
  D_sel = D_sel
)

cat(sprintf(
  "\nSelection model: pi range [%.3f, %.3f], mean = %.3f\n",
  min(pi_obj$pi),
  max(pi_obj$pi),
  mean(pi_obj$pi)
))


##############################################################################
# 9. Complete-case OLS
##############################################################################

cc_data <- data.frame(
  Y = dat$Y,
  X = dat$X,
  Age_SD = age_s,
  Male = dat$Male,
  NHBlack = dat$NHBlack,
  Hispanic = dat$Hispanic,
  Smoker = dat$Smoker,
  R = R
)

cc_fit <- lm(
  Y ~ X + Age_SD + Male + NHBlack + Hispanic + Smoker,
  data = cc_data,
  subset = R == 1
)

cc_coef <- coef(cc_fit)
cc_se <- coef(summary(cc_fit))[, "Std. Error"]

if (anyNA(cc_coef) || anyNA(cc_se)) {
  stop("Complete-case OLS produced aliased coefficients.")
}


##############################################################################
# 10. SIPW
##############################################################################

cat("\nFitting SIPW...\n")

sip <- sipw_fit_glm(
  Y = dat$Y,
  X = X_mat,
  Z = Z_mat,
  R = R,
  family = gaussian(),
  pi_obj = pi_obj
)

cat("SIPW completed.\n")


##############################################################################
# 11. Single-phase BLIPW using BMI
##############################################################################

cat("Fitting BLIPW(BMI)...\n")

bl_bmi <- mp_blipw_glm(
  Y = dat$Y,
  X = X_mat,
  Z = Z_mat,
  W_list = list(D1),
  pi_obj = pi_obj,
  pi_list = list(pi1),
  R = R,
  R_list = list(R1),
  family = gaussian()
)

cat("BLIPW(BMI) completed.\n")


##############################################################################
# 12. Single-phase BLIPW using hip circumference
##############################################################################

cat("Fitting BLIPW(Hip)...\n")

bl_hip <- mp_blipw_glm(
  Y = dat$Y,
  X = X_mat,
  Z = Z_mat,
  W_list = list(D2),
  pi_obj = pi_obj,
  pi_list = list(pi2),
  R = R,
  R_list = list(R2),
  family = gaussian()
)

cat("BLIPW(Hip) completed.\n")


##############################################################################
# 13. MP-BLIPW using BMI + hip circumference
##############################################################################

cat("Fitting MP-BLIPW(BMI + Hip)...\n")

mp <- mp_blipw_glm(
  Y = dat$Y,
  X = X_mat,
  Z = Z_mat,
  W_list = list(D1, D2),
  pi_obj = pi_obj,
  pi_list = list(pi1, pi2),
  R = R,
  R_list = list(R1, R2),
  family = gaussian()
)

cat("MP-BLIPW completed.\n")


##############################################################################
# 14. Build publication-style results
##############################################################################

parameter_names <- c(
  "Intercept",
  "DEXA body fat",
  "Age (SD)",
  "Male sex",
  "Non-Hispanic Black",
  "Hispanic",
  "Ever smoker"
)

expected_p <- length(parameter_names)

if (
  length(cc_coef) != expected_p ||
  length(sip$coef) != expected_p ||
  length(bl_bmi$coef) != expected_p ||
  length(bl_hip$coef) != expected_p ||
  length(mp$coef) != expected_p
) {
  stop("Unexpected number of regression coefficients in one or more models.")
}

z975 <- qnorm(0.975)

make_method_table <- function(method_name, estimate, se, sipw_se) {
  estimate <- as.numeric(estimate)
  se <- as.numeric(se)

  p_value <- 2 * pnorm(-abs(estimate / se))

  if (method_name == "CC-OLS") {
    vr <- rep(NA_real_, length(se))
  } else if (method_name == "SIPW") {
    vr <- rep(0, length(se))
  } else {
    vr <- 100 * (sipw_se^2 - se^2) / sipw_se^2
  }

  data.frame(
    Method = method_name,
    Parameter = parameter_names,
    Estimate = estimate,
    SE = se,
    CI_Lower = estimate - z975 * se,
    CI_Upper = estimate + z975 * se,
    P_value = p_value,
    Variance_Reduction_pct = vr,
    stringsAsFactors = FALSE
  )
}

results_long <- rbind(
  make_method_table("CC-OLS", cc_coef, cc_se, sip$se),
  make_method_table("SIPW", sip$coef, sip$se, sip$se),
  make_method_table("BLIPW (BMI)", bl_bmi$coef, bl_bmi$se, sip$se),
  make_method_table("BLIPW (Hip)", bl_hip$coef, bl_hip$se, sip$se),
  make_method_table("MP-BLIPW", mp$coef, mp$se, sip$se)
)


##############################################################################
# 15. Primary DEXA-body-fat result
##############################################################################

primary_result <- results_long[
  results_long$Parameter == "DEXA body fat",
  ,
  drop = FALSE
]

cat("\n============================================================\n")
cat("PRIMARY ASSOCIATION: DEXA BODY FAT -> HDL-C\n")
cat("============================================================\n")

print(
  transform(
    primary_result,
    Estimate = round(Estimate, 3),
    SE = round(SE, 3),
    CI_Lower = round(CI_Lower, 3),
    CI_Upper = round(CI_Upper, 3),
    P_value = signif(P_value, 3),
    Variance_Reduction_pct = round(Variance_Reduction_pct, 3)
  ),
  row.names = FALSE
)


##############################################################################
# 16. Variance reduction for all parameters
##############################################################################

variance_reduction <- results_long[
  results_long$Method %in%
    c("BLIPW (BMI)", "BLIPW (Hip)", "MP-BLIPW"),
  c("Method", "Parameter", "Variance_Reduction_pct"),
  drop = FALSE
]


##############################################################################
# 17. Sample-structure table
##############################################################################

sample_structure <- data.frame(
  N = N,
  Complete_cases = sum(R),
  Phase1_BMI = sum(R1),
  Phase2_Hip = sum(R2),
  Complete_case_pct = 100 * mean(R),
  Phase1_pct = 100 * mean(R1),
  Phase2_pct = 100 * mean(R2),
  Min_pi_hat = min(pi_obj$pi),
  Max_pi_hat = max(pi_obj$pi)
)

cat("\nSample structure:\n")
print(sample_structure, row.names = FALSE)


##############################################################################
# 18. Save numerical outputs
##############################################################################

dir.create("results", showWarnings = FALSE)

write.csv(
  sample_structure,
  "results/nhanes_sample_structure_hip.csv",
  row.names = FALSE
)

write.csv(
  primary_result,
  "results/nhanes_primary_association_hip.csv",
  row.names = FALSE
)

write.csv(
  variance_reduction,
  "results/nhanes_variance_reduction_hip.csv",
  row.names = FALSE
)

write.csv(
  results_long,
  "results/nhanes_all_coefficients_hip.csv",
  row.names = FALSE
)

results_object <- list(
  sample_structure = sample_structure,
  cc = list(coef = cc_coef, se = cc_se),
  sipw = sip,
  blipw_bmi = bl_bmi,
  blipw_hip = bl_hip,
  mp_blipw = mp,
  primary_result = primary_result,
  all_results = results_long,
  variance_reduction = variance_reduction
)

saveRDS(
  results_object,
  "results/nhanes_mpblipw_results_hip.rds"
)


##############################################################################
# 19. Forest plot for primary coefficient
##############################################################################

dir.create("figures", showWarnings = FALSE)

method_order <- c(
  "CC-OLS",
  "SIPW",
  "BLIPW (BMI)",
  "BLIPW (Hip)",
  "MP-BLIPW"
)

plot_dat <- primary_result[
  match(method_order, primary_result$Method),
  ,
  drop = FALSE
]

png(
  filename = "figures/nhanes_betaX_forestplot_hip.png",
  width = 1200,
  height = 750,
  res = 160
)

ypos <- rev(seq_along(method_order))

x_range <- range(
  c(plot_dat$CI_Lower, plot_dat$CI_Upper),
  finite = TRUE
)

plot(
  plot_dat$Estimate,
  ypos,
  xlim = x_range,
  ylim = c(0.5, length(method_order) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "Regression coefficient for DEXA body fat",
  pch = 19,
  main = "NHANES: DEXA body fat and HDL cholesterol"
)

axis(
  side = 2,
  at = ypos,
  labels = method_order,
  las = 1
)

segments(
  x0 = plot_dat$CI_Lower,
  y0 = ypos,
  x1 = plot_dat$CI_Upper,
  y1 = ypos
)

abline(v = 0, lty = 2)

dev.off()


##############################################################################
# 20. Final checks
##############################################################################

if (qr(D_check)$rank != ncol(D_check)) {
  stop("Unexpected loss of primary-model matrix rank after fitting.")
}

if (any(!is.finite(results_long$Estimate))) {
  stop("Non-finite coefficient estimate detected.")
}

if (any(!is.finite(results_long$SE))) {
  stop("Non-finite standard error detected.")
}

cat("\n============================================================\n")
cat("NHANES HIP ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Results saved in: results/\n")
cat("Figure saved in: figures/\n")
