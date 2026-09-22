##############################################################################
# 03_prepare_nhanes_hip_FIXED.R
#
# Prepare NHANES 2017-2018 data for the FINAL hip-circumference MP-BLIPW
# analysis used in the manuscript.
#
# Final analysis:
#   Y  = HDL cholesterol (LBDHDD)
#   X  = DEXA total percent body fat (DXDTOPF)
#   W1 = BMI (BMXBMI)
#   W2 = Hip circumference (BMXHIP)
#   Z  = age, sex, race/ethnicity indicators, smoking
#
# IMPORTANT:
#   - Adults only: RIDAGEYR >= 20
#   - HDL must be observed
#   - WHR is NOT created or used
#   - Race/sex/smoking recoding is written to work with either numeric NHANES
#     codes or labelled/factor output from nhanesA.
##############################################################################

# Optional: uncomment if you want to force the local project folder.
# setwd("/Users/jyw/Downloads/mpblipw")

if (!requireNamespace("nhanesA", quietly = TRUE)) {
  stop(
    "Package 'nhanesA' is required.\n",
    "Install it once with: install.packages('nhanesA')"
  )
}

##############################################################################
# 1. Helpers
##############################################################################

clean_nhanes_table <- function(x, table_name = "NHANES table") {
  x <- as.data.frame(x)

  if (!"SEQN" %in% names(x)) {
    stop(table_name, " does not contain SEQN.")
  }

  # Remove duplicated column names if present.
  x <- x[, !duplicated(names(x)), drop = FALSE]

  # Ensure one row per participant before merging.
  x <- x[!duplicated(x$SEQN), , drop = FALSE]

  x
}


safe_merge_nhanes <- function(tbls, names_vec) {
  if (length(tbls) != length(names_vec)) {
    stop("tbls and names_vec must have the same length.")
  }

  tbls <- Map(clean_nhanes_table, tbls, names_vec)

  Reduce(
    function(a, b) merge(a, b, by = "SEQN", all = TRUE),
    tbls
  )
}


as_clean_character <- function(x) {
  trimws(as.character(x))
}


as_numeric_code <- function(x) {
  # Works when nhanesA returns numeric values or numeric-looking character data.
  suppressWarnings(as.numeric(as_clean_character(x)))
}


recode_male <- function(x) {
  chr <- tolower(as_clean_character(x))
  num <- as_numeric_code(x)

  out <- rep(NA_real_, length(chr))

  out[(!is.na(num) & num == 1) |
        chr %in% c("male", "1") |
        grepl("^1[[:space:]:-]*male$", chr)] <- 1

  out[(!is.na(num) & num == 2) |
        chr %in% c("female", "2") |
        grepl("^2[[:space:]:-]*female$", chr)] <- 0

  out
}


recode_yes_no <- function(x) {
  chr <- tolower(as_clean_character(x))
  num <- as_numeric_code(x)

  out <- rep(NA_real_, length(chr))

  # Yes = 1
  out[(!is.na(num) & num == 1) |
        chr == "yes" |
        chr == "1" |
        grepl("^1[[:space:]:-]*yes", chr)] <- 1

  # No = 0
  out[(!is.na(num) & num == 2) |
        chr == "no" |
        chr == "2" |
        grepl("^2[[:space:]:-]*no", chr)] <- 0

  out
}


recode_race <- function(x) {
  chr <- tolower(as_clean_character(x))
  num <- as_numeric_code(x)

  # NHANES RIDRETH3 standard codes:
  # 1 = Mexican American
  # 2 = Other Hispanic
  # 3 = Non-Hispanic White
  # 4 = Non-Hispanic Black
  # 6 = Non-Hispanic Asian
  # 7 = Other Race / Multi-Racial

  is_black <- (!is.na(num) & num == 4) |
    grepl("non-hispanic black", chr, fixed = TRUE)

  is_hispanic <- (!is.na(num) & num %in% c(1, 2)) |
    grepl("mexican american", chr, fixed = TRUE) |
    grepl("other hispanic", chr, fixed = TRUE)

  is_asian <- (!is.na(num) & num == 6) |
    grepl("non-hispanic asian", chr, fixed = TRUE)

  # Replace NA logical results with FALSE before numeric conversion.
  is_black[is.na(is_black)] <- FALSE
  is_hispanic[is.na(is_hispanic)] <- FALSE
  is_asian[is.na(is_asian)] <- FALSE

  data.frame(
    NHBlack = as.numeric(is_black),
    Hispanic = as.numeric(is_hispanic),
    Asian = as.numeric(is_asian)
  )
}


##############################################################################
# 2. Download required NHANES 2017-2018 components
##############################################################################

message("Downloading NHANES 2017-2018 component files...")

demo <- nhanesA::nhanes("DEMO_J")
hdl  <- nhanesA::nhanes("HDL_J")
bmx  <- nhanesA::nhanes("BMX_J")
dxx  <- nhanesA::nhanes("DXX_J")
smq  <- nhanesA::nhanes("SMQ_J")
diq  <- nhanesA::nhanes("DIQ_J")


##############################################################################
# 3. Merge by SEQN
##############################################################################

raw <- safe_merge_nhanes(
  tbls = list(demo, hdl, bmx, dxx, smq, diq),
  names_vec = c("DEMO_J", "HDL_J", "BMX_J", "DXX_J", "SMQ_J", "DIQ_J")
)


##############################################################################
# 4. Verify required variables
##############################################################################

required_vars <- c(
  "SEQN",
  "LBDHDD",
  "DXDTOPF",
  "BMXBMI",
  "BMXHIP",
  "RIDAGEYR",
  "RIAGENDR",
  "RIDRETH3",
  "SMQ020",
  "DIQ010"
)

missing_vars <- setdiff(required_vars, names(raw))

if (length(missing_vars) > 0) {
  stop(
    "The following required NHANES variables are missing:\n",
    paste(missing_vars, collapse = ", ")
  )
}


##############################################################################
# 5. Restrict to the FINAL analysis cohort
##############################################################################

# Final report analysis cohort:
#   adults >= 20 years old with HDL-C observed.
raw <- raw[
  !is.na(raw$RIDAGEYR) &
    raw$RIDAGEYR >= 20 &
    !is.na(raw$LBDHDD),
  ,
  drop = FALSE
]

if (nrow(raw) == 0) {
  stop("No observations remain after applying age >= 20 and observed HDL.")
}


##############################################################################
# 6. Recode categorical variables
##############################################################################

race <- recode_race(raw$RIDRETH3)

Male <- recode_male(raw$RIAGENDR)
Smoker <- recode_yes_no(raw$SMQ020)
Diabetes <- recode_yes_no(raw$DIQ010)

# Diabetes is only used in the selection model.
# Match the original working analysis: unresolved/missing diabetes responses
# are coded 0 so the selection-model design matrix keeps all cohort rows.
Diabetes[is.na(Diabetes)] <- 0


##############################################################################
# 7. Sanity checks BEFORE saving
##############################################################################

if (anyNA(Male)) {
  bad <- unique(as_clean_character(raw$RIAGENDR[is.na(Male)]))
  stop(
    "Sex recoding produced missing values. Unrecognized RIAGENDR values: ",
    paste(bad, collapse = ", ")
  )
}

if (anyNA(Smoker)) {
  bad <- unique(as_clean_character(raw$SMQ020[is.na(Smoker)]))
  stop(
    "Smoking recoding produced missing values. Unrecognized SMQ020 values: ",
    paste(bad, collapse = ", ")
  )
}

# These checks specifically prevent the singularity encountered previously.
if (length(unique(race$NHBlack)) < 2) {
  stop(
    "NHBlack is constant after recoding. Check RIDRETH3 coding before analysis.\n",
    "Observed RIDRETH3 values:\n",
    paste(unique(as_clean_character(raw$RIDRETH3)), collapse = " | ")
  )
}

if (length(unique(race$Hispanic)) < 2) {
  stop(
    "Hispanic is constant after recoding. Check RIDRETH3 coding before analysis."
  )
}


##############################################################################
# 8. Construct final analysis dataset
##############################################################################

dat <- data.frame(
  SEQN = raw$SEQN,
  Y = as.numeric(raw$LBDHDD),
  X = as.numeric(raw$DXDTOPF),
  W1 = as.numeric(raw$BMXBMI),
  W2 = as.numeric(raw$BMXHIP),
  Age = as.numeric(raw$RIDAGEYR),
  Male = Male,
  NHBlack = race$NHBlack,
  Hispanic = race$Hispanic,
  Asian = race$Asian,
  Smoker = Smoker,
  Diabetes = Diabetes
)


##############################################################################
# 9. Basic data-quality checks
##############################################################################

if (any(!is.finite(dat$Y))) {
  stop("Y contains non-finite values.")
}

if (any(!is.finite(dat$Age))) {
  stop("Age contains non-finite values.")
}

if (all(is.na(dat$X))) {
  stop("DEXA body-fat variable X is completely missing.")
}

if (all(is.na(dat$W1))) {
  stop("BMI variable W1 is completely missing.")
}

if (all(is.na(dat$W2))) {
  stop("Hip circumference variable W2 is completely missing.")
}


##############################################################################
# 10. Display sample structure
##############################################################################

R  <- as.integer(!is.na(dat$X))
R1 <- pmax(R, as.integer(!is.na(dat$W1)))
R2 <- pmax(R, as.integer(!is.na(dat$W2)))

cat("\n============================================================\n")
cat("FINAL NHANES HIP-CIRCUMFERENCE DATASET\n")
cat("============================================================\n")

cat(sprintf("Analysis cohort (age >=20, HDL observed): N = %d\n", nrow(dat)))
cat(sprintf(
  "DEXA body fat observed:                     %d (%.3f%%)\n",
  sum(!is.na(dat$X)),
  100 * mean(!is.na(dat$X))
))
cat(sprintf(
  "BMI directly observed:                      %d (%.3f%%)\n",
  sum(!is.na(dat$W1)),
  100 * mean(!is.na(dat$W1))
))
cat(sprintf(
  "Hip circumference directly observed:        %d (%.3f%%)\n",
  sum(!is.na(dat$W2)),
  100 * mean(!is.na(dat$W2))
))
cat(sprintf(
  "Phase 1 membership (BMI or CC):              %d (%.3f%%)\n",
  sum(R1),
  100 * mean(R1)
))
cat(sprintf(
  "Phase 2 membership (Hip or CC):              %d (%.3f%%)\n",
  sum(R2),
  100 * mean(R2)
))

cat("\nRace/smoking checks:\n")
print(table(NHBlack = dat$NHBlack, useNA = "ifany"))
print(table(Hispanic = dat$Hispanic, useNA = "ifany"))
print(table(Male = dat$Male, useNA = "ifany"))
print(table(Smoker = dat$Smoker, useNA = "ifany"))


##############################################################################
# 11. Save
##############################################################################

output_file <- "nhanes_mpblipw_data_hip.csv"

write.csv(
  dat,
  file = output_file,
  row.names = FALSE,
  na = ""
)

cat("\nSaved:", normalizePath(output_file, mustWork = FALSE), "\n")


##############################################################################
# 12. Compare with the sample structure reported in the final manuscript
##############################################################################

# The final report reported:
#   N = 4937
#   complete cases = 2142
#   Phase 1 (BMI) = 4870
#   Phase 2 (Hip) = 4706
expected <- c(
  N = 4937,
  CompleteCases = 2142,
  Phase1_BMI = 4870,
  Phase2_Hip = 4706
)

observed <- c(
  N = nrow(dat),
  CompleteCases = sum(R),
  Phase1_BMI = sum(R1),
  Phase2_Hip = sum(R2)
)

cat("\nReported-vs-current sample counts:\n")
print(
  data.frame(
    Quantity = names(expected),
    Reported = as.integer(expected),
    Current = as.integer(observed),
    Match = expected == observed
  ),
  row.names = FALSE
)

if (!all(expected == observed)) {
  warning(
    "Current data counts do not exactly match the final report. ",
    "Do NOT silently change the manuscript. First check the installed nhanesA ",
    "version and the downloaded NHANES tables."
  )
}

cat("\n03_prepare_nhanes_hip_FIXED.R completed successfully.\n")
