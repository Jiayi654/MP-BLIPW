##############################################################################
#  nhanes_mpblipw_analysis.R
#
#  Illustration of MP-BLIPW using NHANES 2017-2018
#
#  SCIENTIFIC QUESTION ====
#  
#  What is the association between adiposity — measured precisely by
#  dual-energy X-ray absorptiometry (DEXA) — and HDL cholesterol,
#  after adjusting for age, sex, race/ethnicity, and smoking?
#
#  MISSING DATA STRUCTURE (two-phase design within the examined sub-sample)
#  ========================================================================
#  The NHANES examination sub-sample forms a natural two-phase study:
#
#    Full cohort (N = 4,380):
#      All examined adults ≥20 with HDL-C measured (lab sub-sample).
#      Y  = HDL cholesterol (mg/dL)                [outcome]
#      Z  = Age, sex, smoking, alcohol, SES        [always observed]
#      W1 = BMI (kg/m²) from examination            [Phase-1 surrogate, ~99%]
#      W2 = Hip circumference (cm) from exam      [Phase-2 surrogate, ~97%]
#
#    Complete cases (n = 1,938, 44%):
#      Subset randomly selected for the DEXA body-composition exam.
#      X  = DEXA % total body fat                  [gold-standard exposure]
#
#  MAPPING TO MP-BLIPW
#  ===================
#  The DEXA measurement is missing for 56% of subjects. Within the CC:
#    cor(X, W1) ≈ 0.57   →  BMI is a moderate surrogate for body fat
#    cor(X, W2) ≈ 0.51   →  Hip circumference is an additional surrogate
#  MP-BLIPW uses both surrogates simultaneously via the two-phase design:
#    Phase 1: BMI observed for nearly all examined subjects (CC ⊆ Phase 1)
#    Phase 2: Hip circumference observed for nearly all (CC ⊆ Phase 2)
#
#  METHODS COMPARED
#  ================
#  1. CC-OLS    — complete-case OLS on the 1,938 DEXA sub-sample
#  2. SIPW      — IPW-corrected estimator using all 4,380 subjects
#  3. BLIPW(W1) — single-phase BLIPW with BMI as working-model surrogate
#  4. BLIPW(W2) — single-phase BLIPW with Hip circumference
#  5. MP-BLIPW  — two-phase BLIPW using both BMI and hip simultaneously
#
#  DATA ACCESS
#  ===========
#  Data loaded via the Python 'nhanes' package (pip install nhanes), which
#  provides the NHANES 2017-2018 harmonised dataset without downloading
#  individual XPT files. An alternative is the nhanesA R package or direct
#  download of XPT files from https://wwwn.cdc.gov/nchs/nhanes/
#
#  REQUIRED FILES
#  ==============
#  • mp_blipw_glm.R   — MP-BLIPW implementation (Supplemental S1)
#  Both are sourced from the working directory.
##############################################################################

## ── 0. Packages and source ────────────────────────────────────────────────
# Install once if needed:
#   pip install nhanes          # Python (for data download)
#   install.packages("reticulate")  # if using Python from R
#   install.packages("nhanesA") # alternative: R-native NHANES access

suppressPackageStartupMessages(library(MASS))

setwd("/Users/kathleenzang/Desktop/missing data/2026")
source("mp_blipw_glm_final.R")


## ── 1. Load and prepare NHANES 2017-2018 data ────────────────────────────

# Read the exported CSV
nhanes <- read.csv("nhanes_mpblipw_data.csv", na.strings = c("","NA","NaN"))


#install.packages("nhanesA")
library(nhanesA)


nhanesTables(data_group='DEMO', year=2017)
nhanesTables(data_group='EXAM', year=2017)

nhanesTableVars('DEMO', 'DEMO_J', nchar=1000)
nhanesTableVars('LAB', 'HDL_J', nchar=1000)
nhanesTableVars('EXAM', 'BPX_J', nchar=1000)
nhanesTableVars('EXAM', 'BMX_J', nchar=1000)
nhanesTableVars('EXAM', 'DXX_J', nchar=1000)
# nhanesTableVars('Q', 'PAQ_J', nchar=1000)
nhanesTableVars('Q', 'SMQ_J', nchar=1000)
nhanesTableVars('Q', 'ALQ_J', nchar=1000)
nhanesTableVars('Q', 'DIQ_J', nchar=1000)

nhanesSearch("liver", ystart=2017, ystop=2018, nchar=1000, ignore.case=TRUE)[,1:3]

demo <- nhanes("DEMO_J"); lab  <- nhanes("HDL_J"); bpx <- nhanes("BPX_J")
bmx  <- nhanes("BMX_J");  dxx  <- nhanes("DXX_J");   smq <- nhanes("SMQ_J")
alq  <- nhanes("ALQ_J");  diq  <- nhanes("DIQ_J");   kiq <- nhanes("KIQ_U_J")
mcq <- nhanes("MCQ_J")
dat  <- Reduce(function(a,b) merge(a,b,by="SEQN",all=TRUE),
               list(demo,lab,bpx,bmx,dxx,smq,alq,diq,kiq,mcq))

# LBDHDD Direct HDL-Cholesterol (mg/dL) 
# DXDTOPF Total Percent Fat (DEXA)
# BMXBMI Body Mass Index (kg/m**2)
# BMXWAIST Waist Circumference (cm)
# BMXHIP Hip Circumference (cm)
# RIDAGEYR Age in years of the participant at the time of screening. Individuals 80 and over are topcoded at 80
# RIAGENDR Gender of the participant.
# INDFMPIR Poverty income ratio (PIR) - a ratio of family income to poverty threshold
#ALQ121 During the past 12 months, about how often did {you/SP} drink any type of alcoholic beverage?\r\nPROBE: How many days per week, per month, or per year did {you/SP} drink?
# SMQ020 These next questions are about cigarette smoking and other tobacco use. {Have you/Has SP} smoked at least 100 cigarettes in {your/his/her} entire life?
# DIQ010 The next questions are about specific medical conditions. {Other than during pregnancy, {have you/has SP}/{Have you/Has SP}} ever been told by a doctor or health professional that {you have/{he/she/SP} has} diabetes or sugar diabetes?
# KIQ022 Have you/Has SP} ever been told by a doctor or other health professional that {you/s/he} had weak or failing kidneys?  Do not include kidney stones, bladder infections, or incontinence.
# MCQ160L Has a doctor or other health professional ever told {you/SP}  that {you/s/he} . . .had any kind of liver condition?
# WTMEC2YR Survey weights
# INDFMPIR

selcols <- c("SEQN","LBDHDD","DXDTOPF","BMXBMI","BMXWAIST","BMXHIP","RIDAGEYR","RIAGENDR","SMQ020","ALQ111","ALQ121","ALQ130","ALQ142","DIQ010","KIQ022","MCQ160L","WTMEC2YR", "INDFMPIR")

dat_ <- subset(dat, RIDAGEYR>=20 & !is.na(LBDHDD), select = selcols)
summary(dat_)
dim(dat_)
sum(dat_$WTMEC2YR)

cat(sprintf("Dataset loaded: N = %d rows, %d columns\n\n",
            nrow(nhanes), ncol(nhanes)))

## ── 2. Restrict to analysis sample ───────────────────────────────────────
#  Analysis cohort: all subjects with HDL-C (Y) measured (already done by export).
#  Within this cohort, DEXA % body fat (X) is missing for ~56%.

dat <- within(dat_, {
  
  # Outcome
  Y <- LBDHDD
  
  # Exposure
  X <- DXDTOPF
  
  # Surrogates
  W1 <- BMXBMI
  W2 <- BMXHIP
  
  # Confounders
  Age <- RIDAGEYR
  
  Male <- ifelse(RIAGENDR == "Male", 1,
                 ifelse(RIAGENDR == "Female", 0, NA))
  
  Smoker <- ifelse(SMQ020 == "Yes", 1,
                   ifelse(SMQ020 == "No", 0, NA))
  
  Income <- INDFMPIR
  
  Alcohol <- ifelse(ALQ121 == "Never in the last year", 0,
                    ifelse(ALQ121 == "1 to 2 times in the last year", 1,
                           ifelse(ALQ121 == "3 to 6 times in the last year", 2,
                                  ifelse(ALQ121 == "7 to 11 times in the last year", 3,
                                         ifelse(ALQ121 == "Once a month", 4,
                                                ifelse(ALQ121 == "2 to 3 times a month", 5,
                                                       ifelse(ALQ121 == "Once a week", 6,
                                                              ifelse(ALQ121 == "2 times a week", 7,
                                                                     ifelse(ALQ121 == "3 to 4 times a week", 8,
                                                                            ifelse(ALQ121 == "Nearly every day", 9,
                                                                                   ifelse(ALQ121 == "Every day", 10, NA)))))))))))
  Alcohol_use <- ifelse(Alcohol == 0, 0,
                        ifelse(!is.na(Alcohol), 1, NA))
})

# Missing values for alcohol use and poverty income ratio
cat("Alcohol_use missing:", sum(is.na(dat$Alcohol_use)),
    "(", round(100 * mean(is.na(dat$Alcohol_use)), 1), "%)\n")

cat("Income missing:", sum(is.na(dat$Income)),
    "(", round(100 * mean(is.na(dat$Income)), 1), "%)\n")
dat <- dat[!is.na(dat$Alcohol_use) & !is.na(dat$Income), ]
cat("N after excluding missing alcohol/income:", nrow(dat), "\n")

cat("Missing Alcohol_use:", sum(is.na(dat$Alcohol_use)), "\n")
cat("Missing Income:", sum(is.na(dat$Income)), "\n")


N   <- nrow(dat)

cat("=== Analysis cohort (HDL-C observed) ===\n")
cat(sprintf("N                          = %d\n", N))
cat(sprintf("X (DEXA %% body fat) obs.  = %d  (%.1f%%)\n",
            sum(!is.na(dat$X)), 100*mean(!is.na(dat$X))))
cat(sprintf("W1 (BMI) observed          = %d  (%.1f%%)\n",
            sum(!is.na(dat$W1)), 100*mean(!is.na(dat$W1))))
cat(sprintf("W2 (waist circ.) observed  = %d  (%.1f%%)\n\n",
             sum(!is.na(dat$W2)), 100*mean(!is.na(dat$W2))))

# Complete-case and phase indicators
R  <- as.integer(!is.na(dat$X))            # CC: has DEXA
R1 <- pmax(R, as.integer(!is.na(dat$W1))) # Phase 1: has BMI  (or CC)
R2 <- pmax(R, as.integer(!is.na(dat$W2))) # Phase 2: has hip (or CC)
stopifnot(all(R[R==1] == 1 & R1[R==1] == 1))  # verify C_CC ⊆ C1
stopifnot(all(R[R==1] == 1 & R2[R==1] == 1))  # verify C_CC ⊆ C2

cat(sprintf("CC  (DEXA + HDL-C):        n  = %d  (%.1f%%)\n", sum(R), 100*mean(R)))
cat(sprintf("Phase 1 (BMI + HDL-C):     n1 = %d  (%.1f%%)\n", sum(R1), 100*mean(R1)))
cat(sprintf("Phase 2 (Hip + HDL-C):   n2 = %d  (%.1f%%)\n\n", sum(R2), 100*mean(R2)))


## ── 3. Descriptive statistics ─────────────────────────────────────────────
cat("=== Descriptive statistics ===\n")
desc <- function(v, label) {
  cat(sprintf("  %-30s mean = %6.2f  SD = %5.2f  n = %d\n",
              label, mean(v, na.rm=TRUE), sd(v, na.rm=TRUE),
              sum(!is.na(v))))
}
desc(dat$Y,   "HDL-C (mg/dL)")
desc(dat$X,   "DEXA % body fat")
desc(dat$W1,  "BMI (kg/m²)")
desc(dat$W2,  "Hip circumference (cm)")
desc(dat$Age, "Age (years)")
cat(sprintf("  %-30s  %.1f%% male\n", "Sex", 100*mean(dat$Male, na.rm=TRUE)))
# cat(sprintf("  %-30s  %.1f%% NH Black\n", "Race: NH Black",
#             100*mean(dat$NHBlack, na.rm=TRUE)))
# cat(sprintf("  %-30s  %.1f%% Hispanic\n", "Race: Hispanic",
#             100*mean(dat$Hispanic, na.rm=TRUE)))
cat(sprintf("  %-30s  %.1f%% ever-smoker\n", "Smoking",
            100*mean(dat$Smoker, na.rm=TRUE)))
cat(sprintf("  %-30s  %.1f%% past-year alcohol use\n", "Alcohol use",
            100*mean(dat$Alcohol_use, na.rm=TRUE)))
desc(dat$Income, "Income ($)")
# cat(sprintf("  %-30s  %.1f%% diabetic\n\n", "Diabetes",
#             100*mean(dat$Diabetes, na.rm=TRUE)))

# Surrogate correlations (in CC subjects)
cc_idx <- R == 1
cat("Correlations among CC subjects (n =", sum(cc_idx), "):\n")
cat(sprintf("  cor(X, W1) = cor(DEXA, BMI)   = %.3f\n",
            cor(dat$X[cc_idx], dat$W1[cc_idx], use = "complete.obs")))
cat(sprintf("  cor(X, W2) = cor(DEXA, Hip) = %.3f\n",
            cor(dat$X[cc_idx], dat$W2[cc_idx], use = "complete.obs")))
cat(sprintf("  cor(W1,W2) = cor(BMI, Hip)  = %.3f\n\n",
            cor(dat$W1[cc_idx], dat$W2[cc_idx], use = "complete.obs")))


## ── 4. Prepare matrices ───────────────────────────────────────────────────
# Standardise age (improves numerical conditioning)
age_s  <- as.vector(scale(dat$Age))
inc_s <- as.vector(scale(dat$Income)) #standardise income

# Impute median for the tiny fraction of W1/W2 missing among non-CC subjects
# (CC subjects always have W1 and W2 by construction of nesting)
W1 <- dat$W1; W1[is.na(W1)] <- median(dat$W1, na.rm = TRUE)
W2 <- dat$W2; W2[is.na(W2)] <- median(dat$W2, na.rm = TRUE)
W2[is.na(W2)] <- median(W2, na.rm = TRUE)

# Set X = 0 for non-CC subjects (values unused in IPW estimating equations)
X  <- dat$X;  X[is.na(X)]   <- 0.0

# Primary design matrix columns: intercept + X + Z
# Note: Diabetes is constant within CC (all DEXA-selected subjects were non-diabetic
# in this cycle due to NHANES exclusion criteria), so it is excluded from the
# primary model but included in the selection model.
X_mat <- matrix(X, ncol = 1)
# Z_mat <- cbind(age_s,         # Age (standardised)
#                dat$Male,      # Sex (1 = male)
#                dat$inc_s,
#                dat$Alcohol_use,
#                dat$Smoker # Ever-smoker
#                )

Z_mat <- cbind(
  age_s,
  dat$Male,
  inc_s,
  dat$Alcohol_use,
  dat$Smoker
)

colnames(Z_mat) <- c(
  "Age",
  "Male",
  "Income",
  "Alcohol",
  "Smoker"
)

dim(Z_mat)

# Working model design matrices (intercept + surrogate)
D1 <- cbind(1, W1)   # Phase 1: intercept + BMI
D2 <- cbind(1, W2)   # Phase 2: intercept + hip circumference

# Phase inclusion probabilities
# DEXA sub-cohort assignment is approximately MCAR conditional on age/sex/race
# (NHANES uses stratified random sub-sampling for body composition)
pi1 <- rep(mean(R1), N); pi1[R == 1] <- 1.0
pi2 <- rep(mean(R2), N); pi2[R == 1] <- 1.0


## ── 5. Selection model ───────────────────────────────────────────────────
# P(DEXA observed | age, sex, race, diabetes) fitted by logistic regression.
# Diabetes can vary in the full cohort and is included here.
D_sel  <- cbind(1, age_s, dat$Male, dat$Alcohol_use, inc_s, dat$Smoker)
pi_obj <- fit_selection(R = R, D_sel = D_sel) #MISSING ALC AND INC removed

cat(sprintf("Selection model (logistic): pi range [%.3f, %.3f], mean = %.3f\n\n",
            min(pi_obj$pi), max(pi_obj$pi), mean(pi_obj$pi)))

pi_obj <- fit_selection(R = R, D_sel = D_sel)

cat(sprintf("Selection model (logistic): pi range [%.3f, %.3f], mean = %.3f\n\n",
            min(pi_obj$pi), max(pi_obj$pi), mean(pi_obj$pi)))
## ── 6. Fit all estimators ─────────────────────────────────────────────────
cat("Fitting estimators...\n")

# (a) Complete-case OLS (benchmark; uses only CC subjects with true DEXA X)
cc_lm <- lm(dat$Y ~ dat$X + age_s + dat$Male + dat$Alcohol_use + dat$Income + dat$Smoker, subset = R == 1)
b_cc  <- coef(cc_lm); se_cc <- coef(summary(cc_lm))[, 2]

# (b) SIPW (corrects for selection into CC via IPW)
sip <- sipw_fit_glm(dat$Y, X_mat, Z_mat, R, gaussian(), pi_obj)

# (c) BLIPW: single-phase, Phase 1 only (BMI as surrogate)
bl1 <- mp_blipw_glm(dat$Y, X_mat, Z_mat,
                    W_list  = list(D1),
                    pi_obj  = pi_obj,
                    pi_list = list(pi1),
                    R       = R,
                    R_list  = list(R1),
                    family  = gaussian())

# (d) BLIPW: single-phase, Phase 2 only (hip circumference)
bl2 <- mp_blipw_glm(dat$Y, X_mat, Z_mat,
                    W_list  = list(D2),
                    pi_obj  = pi_obj,
                    pi_list = list(pi2),
                    R       = R,
                    R_list  = list(R2),
                    family  = gaussian())

# (e) MP-BLIPW: both phases simultaneously (BMI + hip circumference)
mp2 <- mp_blipw_glm(dat$Y, X_mat, Z_mat,
                    W_list  = list(D1, D2),
                    pi_obj  = pi_obj,
                    pi_list = list(pi1, pi2),
                    R       = R,
                    R_list  = list(R1, R2),
                    family  = gaussian())

cat("Done.\n\n")


## ── 7. Results ────────────────────────────────────────────────────────────
z95 <- qnorm(0.975)

# Parameter names (primary model: intercept + X + 7 Z)
pnames <- c("Intercept", "DEXA % body fat", 
            "Age (SD)",
            "Male sex", 
            "Poverty income ratio (PIR)",
            "Alcohol use",
            "Ever smoker"
            )

cat("=======================================================================\n")
cat("  NHANES 2017-2018: HDL-C (mg/dL) ~ DEXA % Body Fat + Covariates\n")
cat("  Model: E[HDL | X, Z] = beta_0 + beta_X * X + beta_Z' * Z  (Gaussian)\n")
cat("=======================================================================\n\n")


for (i in seq_along(pnames)) {
  if( i==1 ){
    cat("--- Full coefficient table ---\n")
    cat(sprintf("%-20s %8s %7s | %8s %7s | %8s %7s | %8s %7s\n",
                "Parameter", "CC", "SE", "SIPW", "SE", "BL1(BMI)", "SE", "MP(K=2)", "SE"))
    cat(strrep("-", 87), "\n")
  }
  cat(sprintf("%-20s %8.3f %7.3f | %8.3f %7.3f | %8.3f %7.3f | %8.3f %7.3f\n",
      pnames[i],
      b_cc[i],      se_cc[i],
      sip$coef[i],  sip$se[i],
      bl1$coef[i],  bl1$se[i],
      mp2$coef[i],  mp2$se[i]))
}
cat(strrep("-", 87), "\n\n")

cat("--- Key coefficient: DEXA % body fat -> HDL-C ---\n")
cat(sprintf("  Interpretation: beta_X = change in HDL-C (mg/dL) per 1 pp increase\n"))
cat(sprintf("                  in DEXA total body fat, adjusting for age, sex, race,\n"))
cat(sprintf("                  and smoking.\n\n"))

fmt_row <- "  %-24s  %7.4f  %6.4f  [%7.4f, %7.4f]  %s\n"
cat(sprintf("  %-24s  %7s  %6s  %15s  %s\n",
    "Method", "Est.", "SE", "95% CI", "SE vs SIPW"))
cat(sprintf("  %s\n", strrep("-", 72)))
v0 <- sip$se[2]^2
for (m in list(
    list("CC-OLS (n=1,938)",      b_cc[2],      se_cc[2],      "reference"),
    list("SIPW",                  sip$coef[2],  sip$se[2],     "—"),
    list("BLIPW (K=1, BMI)",      bl1$coef[2],  bl1$se[2],
         sprintf("-%.1f%%", 100*(v0-bl1$se[2]^2)/v0)),
    list("BLIPW (K=1, Hip)",    bl2$coef[2],  bl2$se[2],
         sprintf("-%.1f%%", 100*(v0-bl2$se[2]^2)/v0)),
    list("MP-BLIPW (K=2)",        mp2$coef[2],  mp2$se[2],
         sprintf("-%.1f%%", 100*(v0-mp2$se[2]^2)/v0)))) {
  cat(sprintf(fmt_row, m[[1]], m[[2]], m[[3]],
              m[[2]] - z95*m[[3]], m[[2]] + z95*m[[3]], m[[4]]))
}

cat("\n--- Variance reductions vs SIPW (all parameters) ---\n")
cat(sprintf("  %-20s %10s %12s %12s\n",
    "Parameter", "BL1(BMI)", "BL2(Hip)", "MP(K=2)"))
cat(sprintf("  %s\n", strrep("-", 58)))
for (i in seq_along(pnames)) {
  v_sip <- sip$se[i]^2
  cat(sprintf("  %-20s %9.1f%% %11.1f%% %11.1f%%\n",
      pnames[i],
      100*(v_sip - bl1$se[i]^2) / v_sip,
      100*(v_sip - bl2$se[i]^2) / v_sip,
      100*(v_sip - mp2$se[i]^2) / v_sip))
}

cat("\n--- 95% confidence intervals for beta_X ---\n")
cat(sprintf("  CC-OLS:         [%.3f, %.3f]\n",
    b_cc[2] - z95*se_cc[2], b_cc[2] + z95*se_cc[2]))
cat(sprintf("  SIPW:           [%.3f, %.3f]  (width = %.3f)\n",
    sip$coef[2]-z95*sip$se[2], sip$coef[2]+z95*sip$se[2], 2*z95*sip$se[2]))
cat(sprintf("  BLIPW(K=1,BMI): [%.3f, %.3f]  (width = %.3f)\n",
    bl1$coef[2]-z95*bl1$se[2], bl1$coef[2]+z95*bl1$se[2], 2*z95*bl1$se[2]))
cat(sprintf("  MP-BLIPW(K=2):  [%.3f, %.3f]  (width = %.3f, %+.1f%% vs SIPW)\n",
    mp2$coef[2]-z95*mp2$se[2], mp2$coef[2]+z95*mp2$se[2], 2*z95*mp2$se[2],
    100*(2*z95*mp2$se[2] - 2*z95*sip$se[2])/(2*z95*sip$se[2])))


## ── 8. Save results ───────────────────────────────────────────────────────
results <- list(
  N = N, n_cc = sum(R), n_p1 = sum(R1), n_p2 = sum(R2),
  pi_range = range(pi_obj$pi),
  cc  = list(coef = b_cc,      se = se_cc),
  sip = list(coef = sip$coef,  se = sip$se),
  bl1 = list(coef = bl1$coef,  se = bl1$se),
  bl2 = list(coef = bl2$coef,  se = bl2$se),
  mp2 = list(coef = mp2$coef,  se = mp2$se),
  var_red = data.frame(
    param   = pnames,
    bl1_pct = 100*(sip$se^2 - bl1$se^2)/sip$se^2,
    bl2_pct = 100*(sip$se^2 - bl2$se^2)/sip$se^2,
    mp2_pct = 100*(sip$se^2 - mp2$se^2)/sip$se^2
  )
)
saveRDS(results, "nhanes_mpblipw_results.rds")
cat("\nResults saved to nhanes_mpblipw_results.rds\n")
cat("\nSession info:\n"); print(sessionInfo(), locale = FALSE)

results <- readRDS("/Users/kathleenzang/Desktop/missing data/2026/nhanes_mpblipw_results.rds")
results
