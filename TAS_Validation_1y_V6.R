# =============================================================================
#  External Validation Pipeline — TAS — 1-YEAR OUTCOME (V6)
# =============================================================================
#  FILES NEEDED IN THE SAME FOLDER AS THIS SCRIPT (from SANT):
#    SANT_transfer_bundle_1yr.rds          (ML model + Platt + threshold;
#                                           from ~/Desktop/LASSO_19.10)
#    PAnCReAS_score_bundle_1yr.rds         (point map, reference unit,
#                                           uncentred S0(1 yr), band cut-offs;
#                                           from the Cox output folder)
#    cox_bootstrap_validation.csv          (optional: SANT corrected C)
#    cox_1yr_harrell_C_optimism_corrected.csv  (optional)
#
#  INSTRUCTIONS FOR THE TASMANIAN DATA TEAM:
#    1. Place this script and the files above in ONE folder
#    2. Set TAS_DATA_PATH below and fill in the VARIABLE MAPPING (Section 3)
#    3. Run: source("TAS_validation_1yr.R")
#    4. Send the ~/Desktop/TAS_Validation_1yr/ folder back to SANT.
# =============================================================================

COHORT <- "TAS"

cat("=============================================================\n")
cat(sprintf("  %s External Validation — 1-year outcome (V6)\n", COHORT))
cat("  ML model: blind validation (SANT Platt + SANT Youden threshold)\n")
cat("  Cox score: PAnCReAS 1-year bands (read from the SANT score bundle)\n")
cat("  Census date: 2022-12-01 (hardcoded, enforced)\n")
cat("=============================================================\n\n")

# =============================================================================
# 0. CONFIGURATION — EDIT THIS SECTION ONLY
# =============================================================================

# ── Full path to the Tasmanian dataset ──────────────────────────────────────
# Use forward slashes or double backslashes on Windows, e.g.:
#   Windows: "C:\\Data\\TAS\\tas_pancreatitis_data.csv"
#   Mac/Linux: "/home/user/data/tas_pancreatitis_data.csv"
TAS_DATA_PATH <- "EDIT_THIS_PATH/tas_data.csv"
DATA_PATH     <- TAS_DATA_PATH

# ── Census date — DO NOT MODIFY ──────────────────────────────────────────────
CENSUS_DATE <- as.Date("2022-12-01")
cat(sprintf("Census date (HARDCODED): %s\n", format(CENSUS_DATE)))

# ── Outcome horizon — DO NOT MODIFY ─────────────────────────────────────────
HORIZON_YRS <- 1        # cancer within 1 year of first AP
FU_CUTOFF   <- CENSUS_DATE - 365   # first AP must be on/before this date
cat(sprintf("Last eligible first-AP date (census - 365 days): %s\n\n", format(FU_CUTOFF)))

# ── Output directories ───────────────────────────────────────────────────────
# out_dir      : release folder - send this one back to SANT
# internal_dir : unsuppressed copies - keep locally, do NOT send
out_dir      <- file.path(path.expand("~"), "Desktop", paste0(COHORT, "_Validation_1yr"))
internal_dir <- paste0(out_dir, "_INTERNAL")
dir.create(out_dir,      showWarnings = FALSE, recursive = TRUE)
dir.create(internal_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output directory (release):  ", out_dir, "\n")
cat("Internal directory (keep):   ", internal_dir, "\n\n")
out <- function(f) file.path(out_dir, f)

EXPECTED_POINT_NAMES <- c("age_50_59", "age_60_69", "age_70plus", "male",
                          "biliary_ap", "alcohol_ap", "unspecified_ap", "other_ap")

# =============================================================================
# 1. PACKAGES
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(survival)
  library(dcurves)
})
cat("Packages loaded.\n\n")

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================

NA_COL <- "__NA__"

get_col <- function(df, varname, na_val = NA_real_) {
  if (length(varname) == 1 && varname == NA_COL) return(rep(na_val, nrow(df)))
  if (!all(varname %in% names(df)))
    stop(sprintf("Column(s) not found in %s dataset: %s\n  Available columns: %s", COHORT,
                 paste(varname[!varname %in% names(df)], collapse = ", "),
                 paste(names(df)[1:min(20, ncol(df))], collapse = ", ")))
  df[[varname]]
}

get_any_col <- function(df, varnames, na_val = 0L) {
  if (length(varnames) == 1 && varnames[1] == NA_COL) return(rep(na_val, nrow(df)))
  if (!all(varnames %in% names(df)))
    stop(sprintf("Column(s) not found in %s dataset: %s", COHORT,
                 paste(varnames[!varnames %in% names(df)], collapse = ", ")))
  as.integer(rowSums(df[varnames] == 1, na.rm = TRUE) > 0)
}

wilson_ci <- function(k, n) {
  z  <- 1.96; p <- k / n
  lo <- (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
  hi <- (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
  c(lo = round(lo, 4), hi = round(hi, 4))
}

evaluate_at_threshold <- function(probs, labels_bin, threshold) {
  pred <- ifelse(probs >= threshold, 1, 0)
  tp <- sum(pred == 1 & labels_bin == 1)
  fp <- sum(pred == 1 & labels_bin == 0)
  fn <- sum(pred == 0 & labels_bin == 1)
  tn <- sum(pred == 0 & labels_bin == 0)
  sens <- tp / max(tp + fn, 1)
  spec <- tn / max(tn + fp, 1)
  ppv  <- ifelse((tp + fp) == 0, NA_real_, tp / (tp + fp))
  npv  <- ifelse((tn + fn) == 0, NA_real_, tn / (tn + fn))
  f1   <- ifelse(is.na(ppv) || (sens + ppv) == 0, NA_real_,
                 2 * sens * ppv / (sens + ppv))
  list(tp=tp, fp=fp, fn=fn, tn=tn,
       sens=sens, spec=spec, ppv=ppv, npv=npv, f1=f1,
       sens_ci=wilson_ci(tp, tp+fn),
       spec_ci=wilson_ci(tn, tn+fp))
}

# ── Cox helper: KM cumulative incidence at time t (+ events by t) ───────────
km_at <- function(d, t) {
  n_ev <- sum(d$event == 1 & d$surv_years <= t)
  if (nrow(d) == 0 || sum(d$event) == 0)
    return(tibble(n_events = n_ev, n_at_risk = NA_integer_, km_pct = 0,
                  km_lo95 = NA_real_, km_hi95 = NA_real_))
  s <- summary(survfit(Surv(surv_years, event) ~ 1, data = d), times = t, extend = TRUE)
  tibble(n_events  = n_ev,
         n_at_risk = s$n.risk,
         km_pct  = (1 - s$surv)  * 100,
         km_lo95 = (1 - s$upper) * 100,   # upper/lower flip for 1 - S(t)
         km_hi95 = (1 - s$lower) * 100)
}

# ── Small-cell suppression (custodian rule: counts 1-9 -> "<10") ─────────────
# Release copies go to out_dir; unsuppressed copies to internal_dir.
# Linked cells that would let a suppressed count be back-calculated are "S".
SDC_MIN <- 10
sdc_small <- function(x) {
  x <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
  !is.na(x) & x > 0 & x < SDC_MIN
}
sdc_lab <- function(x) ifelse(sdc_small(x), "<10", format(x, big.mark = ",", trim = TRUE))

# rules: named list, count column -> linked columns to set to "S".
# drop_cols are removed from the release copy AFTER the flags are computed,
# so a rule can be keyed on a column that is not released (e.g. FP).
sdc <- function(df, rules, label = "", drop_cols = character(0)) {
  out_df <- df %>% mutate(across(everything(), as.character))
  n_supp <- 0
  for (cc in names(rules)) {
    f <- sdc_small(df[[cc]])
    if (!any(f)) next
    out_df[[cc]][f] <- "<10"
    for (lc in intersect(rules[[cc]], names(out_df))) out_df[[lc]][f] <- "S"
    n_supp <- n_supp + sum(f)
  }
  if (n_supp > 0)
    cat(sprintf("  [SDC] %s: %d cell(s) suppressed - check secondary suppression by hand.\n",
                label, n_supp))
  out_df %>% select(-any_of(drop_cols))
}

write_release <- function(df, rules, file, drop_cols = character(0)) {
  write.csv(df, file.path(internal_dir, file), row.names = FALSE, fileEncoding = "UTF-8")
  rel <- sdc(df, rules, file, drop_cols)
  write.csv(rel, out(file), row.names = FALSE, fileEncoding = "UTF-8")
  invisible(rel)
}

# Calibration-plot warning: observed rate x bin N reveals small counts
warn_small_bins <- function(n_events, what) {
  if (any(sdc_small(n_events)))
    warning(sprintf("%s: %d bin/point(s) based on 1-9 cancers - the observed rate ",
                    what, sum(sdc_small(n_events))),
            "can be read off the plot. Check with the data custodian before release.")
}

# Quantile bins that tolerate tied predictions (duplicate cut-points dropped)
quantile_bins <- function(p, n_bins = 10) {
  brks <- unique(quantile(p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(brks) < 2) return(factor(rep("all", length(p))))
  cut(p, breaks = brks, include.lowest = TRUE)
}

# =============================================================================
# 3. VARIABLE MAPPING
#    !! EDIT THIS SECTION TO MATCH YOUR TASMANIAN DATASET COLUMN NAMES !!
#
#    For each variable, replace the string in quotes with the exact column
#    name as it appears in your CSV header row.
#    If a column does not exist in your data, leave it as NA_COL <- "__NA__"
#    and it will be treated as all-zero.
#
#    DATE FORMAT: set DATE_FORMAT to match how dates appear in your data.
#      "%Y-%m-%d"  = 2019-05-14  (Preferred)
#      "%d/%m/%Y"  = 14/05/2019
#      "%m/%d/%Y"  = 05/14/2019
#
#    SEX_MALE_VALUE: set to whatever value codes MALE in your sex column.
#      e.g. 1, "M", "Male", "MALE"
# =============================================================================

# ── Acute pancreatitis ───────────────────────────────────────────────────────
VAR_AP_FLAG          <- "EDIT_ME"   # Binary flag: 1 = confirmed AP patient
VAR_AP_DATE          <- "EDIT_ME"   # Date of first AP episode

# ── Cancer outcome ───────────────────────────────────────────────────────────
VAR_CANCER_DATE      <- "EDIT_ME"   # Date of pancreatic cancer diagnosis (NA if none)
VAR_CANCER_FLAG      <- "EDIT_ME"   # Binary flag: 1 = pancreatic cancer diagnosed
VAR_CANCER_DIFFWEEKS <- "EDIT_ME"   # Weeks from AP to cancer (leave as NA_COL if not available
# — will be calculated automatically from dates)

# ── Chronic pancreatitis ──────────────────────────────
VAR_CHRONIC_FLAG     <- "EDIT_ME"   # Binary flag: 1 = chronic pancreatitis

# ── Death ────────────────────────────────────────────────────────────────────
VAR_DEATH_FLAG       <- "EDIT_ME"   # Binary flag: 1 = deceased
VAR_DEATH_DATE       <- "EDIT_ME"   # Date of death (NA if alive)

# ── Patient characteristics ──────────────────────────────────────────────────
VAR_AGE              <- "EDIT_ME"   # Age at first AP episode (continuous, years)
VAR_SEX              <- "EDIT_ME"   # Sex variable
SEX_MALE_VALUE       <- 1           # Value that codes MALE (edit if "M", "Male", etc.)

# ── AP aetiology ─────────────────────────────────────────────────────────────
VAR_BILIARY          <- "EDIT_ME"   # Binary flag: 1 = biliary AP (ICD K85.1)
VAR_ALCOHOL          <- "EDIT_ME"   # Binary flag: 1 = alcohol AP (ICD K85.2)
VAR_UNSPECIFIED      <- "EDIT_ME"   # Binary flag: 1 = unspecified AP (ICD K85.9)
VAR_OTHER_AP         <- c("EDIT_ME") # One or more columns for "other" AP aetiologies

# ── Date format ───────────────────────────────────────────────────────────────
DATE_FORMAT          <- "%Y-%m-%d"  # Change to match your date format

# =============================================================================
# 4. LOAD SANT OBJECTS
# =============================================================================

cat("=== STEP 4: Loading SANT objects ===\n\n")

# ── 4A. PAnCReAS score bundle (Cox) ─────────────────────────────────────────
if (!file.exists("PAnCReAS_score_bundle_1yr.rds"))
  stop("PAnCReAS_score_bundle_1yr.rds not found.\n",
       "Copy it from the SANT Cox output folder (saved by PAnCReAS_Cox_1-year_FINAL.R, Step 16).")
score_bundle <- readRDS("PAnCReAS_score_bundle_1yr.rds")
if (!identical(as.numeric(score_bundle$horizon_yrs), 1))
  stop("Score bundle is not the 1-year model.")

point_map <- score_bundle$point_map
if (!setequal(names(point_map), EXPECTED_POINT_NAMES))
  stop("Point map names do not match the 1-year model.\n  Expected: ",
       paste(EXPECTED_POINT_NAMES, collapse = ", "),
       "\n  Found:    ", paste(names(point_map), collapse = ", "))
if (any(is.na(point_map))) stop("Point map contains NA values.")

ref_unit <- score_bundle$ref_unit
if (length(ref_unit) != 1 || !is.finite(ref_unit))
  stop("Could not read a single reference unit from the score bundle.")

# Baseline survival must be UNCENTRED (= a score-0 patient), otherwise
# predicted risk, calibration and DCA are all wrong. Checked at SANT (Step 10A).
if (!isFALSE(score_bundle$S0_centered))
  stop("Score bundle S0 is not flagged as uncentred (S0_centered = FALSE). ",
       "Re-run the SANT Cox script to regenerate the bundle.")
s0_1yr <- score_bundle$S0_1yr

# Band cut-offs and score range come from SANT - single source of truth
BAND_LOW_MAX  <- score_bundle$band_low_max
BAND_HIGH_MIN <- score_bundle$band_high_min
SCORE_MIN     <- score_bundle$score_min     # lowest achievable score (no floor needed)
SCORE_MAX     <- score_bundle$score_max
b_sant        <- score_bundle$coefficients[EXPECTED_POINT_NAMES]
if (any(!is.finite(b_sant))) stop("Score bundle does not contain the 8 SANT log-HRs.")

BAND_LABELS <- c(
  Low      = sprintf("Low (\u2264%d)", BAND_LOW_MAX),
  Moderate = sprintf("Moderate (%d\u2013%d)", BAND_LOW_MAX + 1, BAND_HIGH_MIN - 1),
  High     = sprintf("High (\u2265%d)", BAND_HIGH_MIN)
)
BAND_COLOURS <- setNames(c("#0D0887", "#CC4778", "#F0A227"), unname(BAND_LABELS))

assign_band <- function(score) {
  factor(unname(BAND_LABELS[dplyr::case_when(
    score <= BAND_LOW_MAX  ~ "Low",
    score >= BAND_HIGH_MIN ~ "High",
    TRUE                   ~ "Moderate"
  )]), levels = unname(BAND_LABELS))
}

# SANT optimism-corrected 1-year Harrell's C (bootstrap, B = 500)
corrected_C <- corrected_C_score <- NA_real_
if (file.exists("cox_bootstrap_validation.csv")) {
  boot_df     <- read.csv("cox_bootstrap_validation.csv")
  corrected_C <- as.numeric(boot_df$Value[grepl("optimism-corrected", boot_df$Metric,
                                                ignore.case = TRUE)])[1]
} else warning("cox_bootstrap_validation.csv not found - SANT corrected C (Cox model) unavailable.")
if (file.exists("cox_1yr_harrell_C_optimism_corrected.csv")) {
  c10c <- read.csv("cox_1yr_harrell_C_optimism_corrected.csv")
  corrected_C_score <- c10c$Corrected_C[grepl("integer", c10c$Predictor)][1]
  if (is.na(corrected_C)) corrected_C <- c10c$Corrected_C[grepl("linear", c10c$Predictor)][1]
} else warning("cox_1yr_harrell_C_optimism_corrected.csv not found - SANT corrected C (integer score) unavailable.")

cat("PAnCReAS scoring card (SANT derivation, 1-year model):\n")
cat(sprintf("  Age 50-59        %+d\n", point_map[["age_50_59"]]))
cat(sprintf("  Age 60-69        %+d\n", point_map[["age_60_69"]]))
cat(sprintf("  Age >= 70        %+d\n", point_map[["age_70plus"]]))
cat(sprintf("  Male             %+d\n", point_map[["male"]]))
cat(sprintf("  Biliary AP       %+d\n", point_map[["biliary_ap"]]))
cat(sprintf("  Alcohol AP       %+d\n", point_map[["alcohol_ap"]]))
cat(sprintf("  Unspecified AP   %+d\n", point_map[["unspecified_ap"]]))
cat(sprintf("  Other AP         %+d\n", point_map[["other_ap"]]))
cat(sprintf("  Bands: %s | %s | %s\n", BAND_LABELS[["Low"]], BAND_LABELS[["Moderate"]],
            BAND_LABELS[["High"]]))
cat(sprintf("  Achievable score range %d to %d (no floor applied)\n", SCORE_MIN, SCORE_MAX))
cat(sprintf("  Reference unit   %.4f | S0(1 yr) = %.6f (uncentred)\n", ref_unit, s0_1yr))
cat(sprintf("  SANT optimism-corrected 1-yr C: Cox model %.3f | integer score %.3f\n\n",
            corrected_C, corrected_C_score))

# ── 4B. ML transfer bundle (LASSO or XGBoost) ───────────────────────────────
if (!file.exists("SANT_transfer_bundle_1yr.rds"))
  stop("SANT_transfer_bundle_1yr.rds not found - this script needs the 1-YEAR bundle.")
bundle <- readRDS("SANT_transfer_bundle_1yr.rds")
if (!identical(as.numeric(bundle$outcome_horizon_yrs), 1))
  stop("Bundle is not the 1-year model (outcome_horizon_yrs != 1).")
ML_NAME         <- bundle$best_model_name
platt_sant      <- bundle$platt_model_final
feature_names   <- bundle$feature_names
sant_prevalence <- bundle$sant_prevalence
THRESHOLD_A     <- bundle$best_threshold          # SANT Youden, calibrated scale
sant_auroc_ml   <- if (!is.null(bundle$sant_final_auroc)) bundle$sant_final_auroc else NA_real_
model_final     <- bundle$model_final

# Predict the raw (pre-Platt) probability. Columns are forced into the
# bundle's feature order first.
ml_predict_raw <- function(X_mat) {
  miss <- setdiff(feature_names, colnames(X_mat))
  if (length(miss)) stop("Predictor column(s) missing for the ML model: ", paste(miss, collapse = ", "))
  X_mat <- X_mat[, feature_names, drop = FALSE]
  if (ML_NAME == "LASSO") {
    lam <- if (!is.null(bundle$lasso_lambda)) bundle$lasso_lambda else model_final$lambda[1]
    as.numeric(predict(model_final, newx = X_mat, type = "response", s = lam)[, 1])
  } else if (ML_NAME == "XGBoost") {
    as.numeric(predict(model_final, xgboost::xgb.DMatrix(X_mat)))
  } else if (ML_NAME == "Random Forest") {
    as.numeric(predict(model_final, as.data.frame(X_mat), type = "prob")[, "Cancer"])
  } else stop("Unknown model type in bundle: ", ML_NAME)
}

pkg <- c(LASSO = "glmnet", XGBoost = "xgboost", `Random Forest` = "randomForest")[[ML_NAME]]
if (!requireNamespace(pkg, quietly = TRUE))
  stop(sprintf("Package '%s' is needed for the %s model: install.packages(\"%s\")", pkg, ML_NAME, pkg))
suppressPackageStartupMessages(library(pkg, character.only = TRUE))

# Load check: one trivial prediction. A failure here is reported as-is
# (plain predict(); the old xgboost::predict() call does not exist and gave
# a misleading "version mismatch" message).
if (ML_NAME == "XGBoost" && file.exists("model_final.ubj")) {
  model_final <- xgboost::xgb.load("model_final.ubj")
  cat("XGBoost model loaded from: model_final.ubj (portable format)\n")
}
test_mat <- matrix(0, nrow = 1, ncol = length(feature_names),
                   dimnames = list(NULL, feature_names))
tryCatch(invisible(ml_predict_raw(test_mat)), error = function(e) {
  stop(sprintf("The %s model in the bundle could not make a prediction.\n  Error: %s\n%s",
               ML_NAME, conditionMessage(e),
               if (ML_NAME == "XGBoost")
                 "  If this is an XGBoost version mismatch, ask SANT for model_final.ubj and place it here.\n"
               else ""))
})
cat(sprintf("ML model loaded: %s (%s) | features: %s\n", ML_NAME, class(model_final)[1],
            paste(feature_names, collapse = ", ")))
cat(sprintf("SANT prevalence: %.2f%% | threshold (SANT Youden, from bundle): %.4f\n\n",
            sant_prevalence * 100, THRESHOLD_A))

# =============================================================================
# 5. LOAD TASMANIAN DATA
# =============================================================================

cat(sprintf("=== STEP 5: Loading %s data ===\n\n", COHORT))

if (!file.exists(DATA_PATH)) {
  stop(sprintf(
    "Tasmanian data file not found at:\n  %s\nPlease update TAS_DATA_PATH in Section 0.",
    DATA_PATH
  ))
}

raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
cat(sprintf("Raw %s rows: %d\n", COHORT, nrow(raw)))
cat(sprintf("Columns detected (%d): %s\n\n",
            ncol(raw),
            paste(names(raw)[1:min(15, ncol(raw))], collapse = ", ")))

# ── Pre-processing: construct death flag if not present ──────────────────────
# If your death flag is already coded 0/1, skip this block or comment it out.
# If death is indicated by the presence of a date, use this:
if (VAR_DEATH_FLAG != NA_COL && !VAR_DEATH_FLAG %in% names(raw)) {
  raw$nv_death <- ifelse(is.na(raw[[VAR_DEATH_DATE]]) |
                           raw[[VAR_DEATH_DATE]] == "", 0, 1)
  VAR_DEATH_FLAG <- "nv_death"
  cat("NOTE: Death flag constructed from death date column.\n")
}

# ── Pre-processing: construct cancer flag if not present ─────────────────────
# If cancer flag is indicated by presence of a cancer date, use this:
if (VAR_CANCER_FLAG != NA_COL && !VAR_CANCER_FLAG %in% names(raw)) {
  raw$nv_cancer <- ifelse(is.na(raw[[VAR_CANCER_DATE]]) |
                            raw[[VAR_CANCER_DATE]] == "", 0, 1)
  VAR_CANCER_FLAG <- "nv_cancer"
  cat("NOTE: Cancer flag constructed from cancer date column.\n")
}

# ── Pre-processing: sex recode ───────────────────────────────────────────────
# Check and recode sex to numeric if needed.
# Edit this block to match how sex is coded in your data.
cat("Sex variable distribution (before recode):\n")
print(table(raw[[VAR_SEX]], useNA = "always"))

# EXAMPLE: If sex is "Male"/"Female" strings, recode to 1/2:
# raw[[VAR_SEX]] <- ifelse(raw[[VAR_SEX]] == "Male",   1,
#                   ifelse(raw[[VAR_SEX]] == "Female", 2, NA))

# EXAMPLE: If sex is "M"/"F":
# raw[[VAR_SEX]] <- ifelse(raw[[VAR_SEX]] == "M", 1,
#                   ifelse(raw[[VAR_SEX]] == "F", 2, NA))

# After recoding, confirm SEX_MALE_VALUE matches the male code in your data.
cat("Sex variable distribution (after recode — confirm 1=Male):\n")
print(table(raw[[VAR_SEX]], useNA = "always"))
cat("\n")

# ── Pre-processing: "other AP" composite flag ─────────────────────────────────
# If your data has multiple "other AP aetiology" columns, combine them here.
# Edit VAR_OTHER_AP in Section 3 to list all relevant column names.
# If only one column, VAR_OTHER_AP <- c("column_name") is fine.

# =============================================================================
# 6. BUILD ANALYTICAL DATASET (with flowchart)
# =============================================================================
# Exclusions, applied in this order (same rules as the SANT derivation, plus
# the chronic-pancreatitis exclusion used in the validation cohorts):
#   1. missing/invalid first-AP date
#   2. chronic pancreatitis
#   3. pancreatic cancer diagnosed BEFORE first AP
#   4. pancreatic cancer flagged but no diagnosis date (cannot be timed)
#   5. first AP after census - 365 days (< 1 year possible follow-up, so a
#      "no cancer within 1 year" label could not be verified)
#   6. missing predictor data
# =============================================================================

cat("=== STEP 6: Building analytical dataset ===\n\n")

d0 <- raw %>%
  filter(get_col(., VAR_AP_FLAG) == 1) %>%
  mutate(
    nv_k85_firstage = get_col(., VAR_AGE),
    ap_date         = as.Date(get_col(., VAR_AP_DATE),     format = DATE_FORMAT),
    cancer_date     = as.Date(get_col(., VAR_CANCER_DATE), format = DATE_FORMAT),
    death_date      = as.Date(get_col(., VAR_DEATH_DATE),  format = DATE_FORMAT),
    cancer_flag     = get_col(., VAR_CANCER_FLAG),
    chronic_flag    = get_col(., VAR_CHRONIC_FLAG),
    death_flag      = get_col(., VAR_DEATH_FLAG),
    
    # Weeks from AP to cancer (calculated from dates if column not available)
    diffweeks = case_when(
      VAR_CANCER_DIFFWEEKS != NA_COL ~
        as.numeric(get_col(., VAR_CANCER_DIFFWEEKS)),
      !is.na(cancer_date) & !is.na(ap_date) ~
        as.numeric(cancer_date - ap_date) / 7,
      TRUE ~ NA_real_
    ),
    
    # Age categories (identical to SANT derivation)
    age_cat = case_when(
      nv_k85_firstage < 50                          ~ 0,
      nv_k85_firstage >= 50 & nv_k85_firstage < 60 ~ 1,
      nv_k85_firstage >= 60 & nv_k85_firstage < 70 ~ 2,
      nv_k85_firstage >= 70                         ~ 3
    ),
    
    # Predictors (also used by the ML model)
    male           = as.integer(get_col(., VAR_SEX)         == SEX_MALE_VALUE),
    biliary_ap     = as.integer(get_col(., VAR_BILIARY)      == 1),
    alcohol_ap     = as.integer(get_col(., VAR_ALCOHOL)      == 1),
    unspecified_ap = as.integer(get_col(., VAR_UNSPECIFIED)  == 1),
    other_ap       = get_any_col(., VAR_OTHER_AP, na_val = 0L),
    
    # Exclusion flags
    ca_flag   = cancer_flag %in% 1,
    cp_flag   = chronic_flag %in% 1,
    prior_ca  = ca_flag & !is.na(diffweeks) & diffweeks < 0,
    ca_nodate = ca_flag & !prior_ca & (is.na(cancer_date) | is.na(diffweeks))
  )

d1 <- d0 %>% filter(!is.na(ap_date))
d2 <- d1 %>% filter(!cp_flag)
d3 <- d2 %>% filter(!prior_ca)
d4 <- d3 %>% filter(!ca_nodate)
d5 <- d4 %>% filter(ap_date <= FU_CUTOFF)
d6 <- d5 %>% drop_na(age_cat, male, biliary_ap, alcohol_ap, unspecified_ap)

dat_val <- d6 %>%
  mutate(
    # Cancer AFTER AP and ON/BEFORE census (chronic pancreatitis already excluded)
    event_full = as.integer(ca_flag & !is.na(cancer_date) & diffweeks >= 0 &
                              cancer_date <= CENSUS_DATE),
    
    # Follow-up end date: cancer date > death date (on/before census) > census date
    end_date = case_when(
      event_full == 1                                                  ~ cancer_date,
      death_flag == 1 & !is.na(death_date) & death_date <= CENSUS_DATE ~ death_date,
      TRUE                                                             ~ CENSUS_DATE   # HARDCODED
    ),
    
    # Survival time (floor at 1 day for same-admission cancer diagnoses)
    surv_years_full = pmax(as.numeric(end_date - ap_date) / 365.25, 1 / 365.25),
    
    # -- 1-YEAR OUTCOME (model outcome) --------------------------------------
    event      = as.integer(event_full == 1 & surv_years_full <= HORIZON_YRS),
    surv_years = pmin(surv_years_full, HORIZON_YRS)
  )

flowchart <- tibble(
  Step = c("AP patients (AP flag = 1)",
           "Excluded: missing or invalid first-AP date",
           "Excluded: chronic pancreatitis",
           "Excluded: pancreatic cancer diagnosed before first AP",
           "Excluded: pancreatic cancer recorded with no diagnosis date",
           sprintf("Excluded: first AP after %s (< 1 year possible follow-up)", format(FU_CUTOFF)),
           "Excluded: missing predictor data",
           "Analytic cohort",
           "  Pancreatic cancer within 1 year (outcome)",
           "  Pancreatic cancer at any time to census (reference)"),
  N = c(nrow(d0), nrow(d0) - nrow(d1), nrow(d1) - nrow(d2), nrow(d2) - nrow(d3),
        nrow(d3) - nrow(d4), nrow(d4) - nrow(d5), nrow(d5) - nrow(d6), nrow(dat_val),
        sum(dat_val$event), sum(dat_val$event_full))
)
cat(sprintf("=== FLOWCHART (%s) ===\n", COHORT))
print(as.data.frame(flowchart), row.names = FALSE)
write_release(flowchart, list(N = character(0)), sprintf("Flowchart_%s_1yr.csv", COHORT))
cat(sprintf("Saved: Flowchart_%s_1yr.csv\n", COHORT))
n_ca_after_census <- sum(dat_val$ca_flag & !is.na(dat_val$cancer_date) &
                           dat_val$cancer_date > CENSUS_DATE)
cat(sprintf("Cancers diagnosed after the census date (censored, not events): %d\n\n",
            n_ca_after_census))

# Binarise aetiology flags if any are > 1 (safety check)
for (v in c("biliary_ap","alcohol_ap","unspecified_ap","other_ap")) {
  if (max(dat_val[[v]], na.rm = TRUE) > 1) {
    dat_val[[v]] <- as.integer(dat_val[[v]] > 0)
    cat(sprintf("AUTO-CORRECTION: %s binarised (values > 1 detected).\n", v))
  }
}

# ── PAnCReAS score, bands, predicted risk (after binarisation) ──────────────
dat_val <- dat_val %>%
  mutate(
    age_50_59  = as.integer(age_cat == 1),
    age_60_69  = as.integer(age_cat == 2),
    age_70plus = as.integer(age_cat == 3),
    cox_score =
      age_50_59      * point_map[["age_50_59"]]      +
      age_60_69      * point_map[["age_60_69"]]      +
      age_70plus     * point_map[["age_70plus"]]     +
      male           * point_map[["male"]]           +
      biliary_ap     * point_map[["biliary_ap"]]     +
      alcohol_ap     * point_map[["alcohol_ap"]]     +
      unspecified_ap * point_map[["unspecified_ap"]] +
      other_ap       * point_map[["other_ap"]],      # no floor (see SANT Step 6)
    score_band = assign_band(cox_score),
    risk_band  = score_band,
    # SANT-predicted absolute risk (1-year): uncentred S0 so score 0 = baseline
    pred_1yr = 1 - s0_1yr ^ exp(cox_score * ref_unit),
    # Clinical rule comparator (non-gallstone/alcohol AP OR age >= 50)
    clinical_rule = as.integer((biliary_ap == 0 & alcohol_ap == 0) | age_cat >= 1)
  )
if (min(dat_val$cox_score) < SCORE_MIN || max(dat_val$cox_score) > SCORE_MAX)
  stop("Scores fall outside the SANT achievable range - check the predictor coding.")

y_val <- dat_val$event

cat("=== DATASET SUMMARY ===\n")
cat(sprintf("  Total N (analytic cohort): %d\n",   nrow(dat_val)))
cat(sprintf("  Cancer events, 1 year:     %d (%.2f%%)  <- outcome\n", sum(y_val), mean(y_val)*100))
cat(sprintf("  Cancer events, any time:   %d (%.2f%%)\n", sum(dat_val$event_full), mean(dat_val$event_full)*100))
cat(sprintf("  No cancer within 1 year:   %d\n",   sum(y_val == 0)))
cat(sprintf("  Median follow-up (uncapped): %.1f years\n", median(dat_val$surv_years_full)))
cat(sprintf("  Died within 1 year without cancer: %d\n",
            sum(dat_val$event == 0 & dat_val$surv_years_full < HORIZON_YRS)))
cat(sprintf("  Same-day cancer diagnoses (1 day): %d\n",
            sum(dat_val$event == 1 & dat_val$surv_years <= 0.01)))
cat(sprintf("  Census date applied:       %s\n\n",  format(CENSUS_DATE)))

cat("Predictor distribution (% = 1):\n")
for (v in c("age_50_59","age_60_69","age_70plus","male",
            "biliary_ap","alcohol_ap","unspecified_ap","other_ap"))
  cat(sprintf("  %-16s %.1f%%\n", v, mean(dat_val[[v]], na.rm=TRUE)*100))
cat("Aetiology codes per patient (biliary + alcohol + unspecified + other):\n")
print(table(with(dat_val, biliary_ap + alcohol_ap + unspecified_ap + other_ap)))
cat(sprintf("\nScore range: %d to %d (achievable %d to %d)\n\n",
            min(dat_val$cox_score), max(dat_val$cox_score), SCORE_MIN, SCORE_MAX))

# =============================================================================
# 7. MISSING DATA AUDIT (analytic cohort)
# =============================================================================

cat("=== STEP 7: Missing data audit ===\n\n")

miss_df <- dat_val %>%
  select(event, nv_k85_firstage, male, biliary_ap,
         alcohol_ap, unspecified_ap, other_ap) %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to="variable", values_to="n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(dat_val) * 100, 2))
print(miss_df)
cat(sprintf("(%d patients were excluded for missing predictor data - see flowchart.)\n\n",
            nrow(d5) - nrow(d6)))

# =============================================================================
# 8. ML MODEL: RAW PREDICTIONS
# =============================================================================

cat(sprintf("=== STEP 8: %s raw predictions ===\n\n", ML_NAME))

X_val_cont <- dat_val %>%
  select(nv_k85_firstage, male, biliary_ap,
         alcohol_ap, unspecified_ap, other_ap)

if (any(is.na(X_val_cont$nv_k85_firstage))) {
  cat("WARNING: Continuous age contains NAs — using age_cat approximation.\n")
  X_val_cont$nv_k85_firstage <- dat_val$age_cat * 10 + 45
}

X_val_mat <- model.matrix(~ . - 1, data = X_val_cont)
X_val_mat <- X_val_mat[, feature_names, drop = FALSE]   # force bundle column order

prob_raw  <- ml_predict_raw(X_val_mat)

cat(sprintf("Raw %s prob — median: %.4f | max: %.4f\n\n",
            ML_NAME, median(prob_raw), max(prob_raw)))

roc_raw      <- pROC::roc(y_val, prob_raw, quiet = TRUE)
auroc_raw    <- as.numeric(pROC::auc(roc_raw))
auroc_raw_ci <- as.numeric(pROC::ci.auc(roc_raw, method = "delong"))
cat(sprintf("Raw %s AUROC (leakage-free): %.4f (95%% CI: %.3f-%.3f)\n",
            ML_NAME, auroc_raw, auroc_raw_ci[1], auroc_raw_ci[3]))
cat("This is the primary leakage-free discrimination estimate.\n\n")

# =============================================================================
# 9. ANALYSIS A — BLIND EXTERNAL VALIDATION
# =============================================================================

cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("  %s — BLIND EXTERNAL VALIDATION (1-year outcome)\n", toupper(ML_NAME)))
cat(sprintf("  SANT Platt model applied directly. No %s labels used.\n", COHORT))
cat(sprintf("  Threshold: SANT Youden optimum (%.4f)\n", THRESHOLD_A))
cat("════════════════════════════════════════════════════════════\n\n")

prob_A <- as.numeric(predict(platt_sant,
                             newdata = data.frame(prob = prob_raw),
                             type    = "response"))

cat(sprintf("SANT Platt applied: median = %.4f | max = %.4f\n",
            median(prob_A), max(prob_A)))

if (max(prob_A) < THRESHOLD_A) {
  cat(sprintf(
    "\nNOTE: SANT Platt compresses %s probabilities (max = %.4f) below\n",
    COHORT, max(prob_A)))
  cat(sprintf("the SANT threshold (%.4f). Sensitivity at SANT threshold = 0%%.\n",
              THRESHOLD_A))
  cat("This is a valid finding: discrimination may be preserved (AUROC)\n")
  cat("while calibration does not transfer. Reported as-is.\n\n")
}

roc_A      <- pROC::roc(y_val, prob_A, quiet = TRUE)
auroc_A    <- as.numeric(pROC::auc(roc_A))
auroc_A_ci <- as.numeric(pROC::ci.auc(roc_A, method = "delong"))
perf_A     <- evaluate_at_threshold(prob_A, y_val, THRESHOLD_A)

pr_A       <- PRROC::pr.curve(scores.class0 = prob_A[y_val==1],
                              scores.class1 = prob_A[y_val==0], curve = FALSE)
brier_A    <- mean((prob_A - y_val)^2)
brier_null <- mean(y_val) * (1 - mean(y_val))
brier_skA  <- 1 - (brier_A / brier_null)

cat(sprintf("AUROC (SANT Platt):    %.4f (95%% CI: %.3f-%.3f)\n",
            auroc_A, auroc_A_ci[1], auroc_A_ci[3]))
if (!is.na(sant_auroc_ml)) cat(sprintf("AUROC drop from SANT:  %.4f\n", sant_auroc_ml - auroc_A))
cat(sprintf("AUPRC:                 %.4f\n", pr_A$auc.integral))
cat(sprintf("Brier:                 %.4f (skill: %.4f)\n", brier_A, brier_skA))
cat(sprintf("Sensitivity at %.4f: %.4f (%d / %d cancers caught)\n",
            THRESHOLD_A, perf_A$sens, perf_A$tp, perf_A$tp + perf_A$fn))
cat(sprintf("Specificity at %.4f: %.4f\n",    THRESHOLD_A, perf_A$spec))
cat(sprintf("PPV:                   %.4f\n\n", ifelse(is.na(perf_A$ppv), NaN, perf_A$ppv)))

results_A <- tibble(
  Analysis         = "A — Blind external validation",
  Model            = ML_NAME,
  Platt_model      = sprintf("SANT-derived (no %s labels used)", COHORT),
  Threshold        = THRESHOLD_A,
  Threshold_source = "SANT Youden optimum on calibrated 1-year probability (from bundle)",
  Census_date      = format(CENSUS_DATE),
  N                = nrow(dat_val),
  Cancer_events    = sum(y_val),
  AUROC            = round(auroc_A,         4),
  AUROC_lo         = round(auroc_A_ci[1],   4),
  AUROC_hi         = round(auroc_A_ci[3],   4),
  AUROC_raw        = round(auroc_raw,       4),
  SANT_AUROC       = sant_auroc_ml,
  AUROC_drop       = round(sant_auroc_ml - auroc_A, 4),
  AUPRC            = round(pr_A$auc.integral, 4),
  Sensitivity      = round(perf_A$sens, 4),
  Specificity      = round(perf_A$spec, 4),
  PPV              = round(perf_A$ppv,  4),
  NPV              = round(perf_A$npv,  4),
  Brier            = round(brier_A,     4),
  Brier_skill      = round(brier_skA,   4),
  TP = perf_A$tp, FP = perf_A$fp, FN = perf_A$fn, TN = perf_A$tn,
  Calibration_note = ifelse(
    max(prob_A) < THRESHOLD_A,
    sprintf("SANT Platt compressed %s probs below threshold — sensitivity = 0%% is a valid finding", COHORT),
    "SANT Platt adequate"
  )
)
write_release(results_A,
              list(TP = c("Sensitivity", "PPV", "NPV", "FN"),
                   FN = c("Sensitivity", "PPV", "NPV", "TP"),
                   FP = c("Specificity", "PPV", "TN"),
                   TN = c("Specificity", "NPV", "FP")),
              "AnalysisA_blind_validation.csv")
cat("Saved: AnalysisA_blind_validation.csv\n\n")

# ── Table S6 — blind validation column (manuscript layout, <10 rule applied) ─
tp <- perf_A$tp; fp <- perf_A$fp; fn <- perf_A$fn; tn <- perf_A$tn
s6 <- function(cells, txt) if (any(sdc_small(cells))) "S" else txt
nni_A <- ifelse(is.na(perf_A$ppv) || perf_A$ppv == 0, NA, 1 / perf_A$ppv)
table_s6 <- tibble(
  Metric = c("Model", "N (AP)", "Cancer events within 1 year, n (%)",
             "Cancer events, any time, n (%)",
             "Platt calibration model", "Threshold applied", "AUROC (95% CI)",
             "AUPRC", "Sensitivity", "Specificity", "NNI", "PPV", "NPV",
             "Cancers caught / total", "Patients flagged, n (%)",
             "Brier score (calibrated)", "Brier skill score"),
  Value = c(
    ML_NAME,
    format(nrow(dat_val), big.mark = ","),
    if (sdc_small(sum(y_val))) "<10" else sprintf("%d (%.2f%%)", sum(y_val), 100 * mean(y_val)),
    if (sdc_small(sum(dat_val$event_full))) "<10" else
      sprintf("%d (%.2f%%)", sum(dat_val$event_full), 100 * mean(dat_val$event_full)),
    "SANT-derived",
    sprintf("%.4f", THRESHOLD_A),
    sprintf("%.3f (%.3f–%.3f)", auroc_A_ci[2], auroc_A_ci[1], auroc_A_ci[3]),
    sprintf("%.4f", pr_A$auc.integral),
    s6(c(tp, fn), sprintf("%.1f%%", 100 * perf_A$sens)),
    s6(c(tn, fp), sprintf("%.1f%%", 100 * perf_A$spec)),
    s6(c(tp, fp, fn), ifelse(is.na(nni_A), "—", sprintf("%.1f", nni_A))),
    s6(c(tp, fp, fn), sprintf("%.2f%%", 100 * perf_A$ppv)),
    s6(c(tn, fn, tp), sprintf("%.2f%%", 100 * perf_A$npv)),
    sprintf("%s / %d", if (sdc_small(tp)) "<10" else if (sdc_small(fn)) "S" else tp, tp + fn),
    if (sdc_small(tp + fp)) "<10" else if (sdc_small(fp)) "S" else
      sprintf("%s (%.1f%%)", format(tp + fp, big.mark = ","), 100 * (tp + fp) / nrow(dat_val)),
    sprintf("%.4f", brier_A),
    sprintf("%.4f", brier_skA))
)
names(table_s6)[2] <- sprintf("%s Blind Validation", COHORT)
cat(sprintf("=== TABLE S6 — %s Blind Validation column ===\n", COHORT))
print(as.data.frame(table_s6), row.names = FALSE)
write.csv(table_s6, out(sprintf("Table_S6_%s_1yr.csv", COHORT)), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Saved: Table_S6_%s_1yr.csv\n\n", COHORT))

# =============================================================================
# 11. PAnCReAS COX SCORE (1-YEAR MODEL) — RISK STRATIFICATION,
#     DISCRIMINATION, CALIBRATION
# =============================================================================

cat("=== STEP 11: PAnCReAS Cox score (1-year model) ===\n\n")

score_dist <- dat_val %>%
  group_by(cox_score) %>%
  summarise(n = n(), n_cancer = sum(event),
            cancer_rate_pct = round(100 * mean(event), 2), .groups = "drop") %>%
  mutate(band = assign_band(cox_score))
print(score_dist, n = Inf)
write_release(score_dist,
              list(n = c("n_cancer", "cancer_rate_pct"),
                   n_cancer = c("cancer_rate_pct")),
              sprintf("Cox_score_distribution_%s.csv", COHORT))

band_counts <- dat_val %>%
  group_by(risk_band, .drop = FALSE) %>%
  summarise(N = n(), Pct_of_cohort = round(100 * n() / nrow(dat_val), 1),
            N_cancer = sum(event),
            Rate_1yr_pct = round(100 * mean(event), 2),
            N_cancer_any_time = sum(event_full), .groups = "drop")
cat("\nBand sizes and 1-year cancer rates:\n")
print(band_counts)
write_release(band_counts,
              list(N = c("Pct_of_cohort"),
                   N_cancer = c("Rate_1yr_pct"),
                   N_cancer_any_time = character(0)),
              sprintf("Cox_band_counts_%s.csv", COHORT))

# -- Direction check (report only - never flip) -------------------------------
dir_fit <- coxph(Surv(surv_years, event) ~ cox_score, data = dat_val)
hr_per_pt <- exp(coef(dir_fit))
cat(sprintf("\nDirection check: HR per score point = %.3f (should be > 1)\n", hr_per_pt))
rates <- band_counts$Rate_1yr_pct
if (hr_per_pt < 1 || isTRUE(rates[3] < rates[1])) {
  cat(sprintf("WARNING: score is inversely associated with cancer in %s.\n", COHORT),
      "        Check SEX_MALE_VALUE and the aetiology/age mapping BEFORE anything else.\n",
      "        Do NOT reverse the score - report the result as observed.\n")
}
cat("\n")

# ── 11A. Discrimination ──────────────────────────────────────────────────────
cat("=== STEP 11A: Discrimination ===\n\n")
roc_cox <- pROC::roc(dat_val$event, dat_val$cox_score, quiet = TRUE, direction = "<")
auc_ci  <- as.numeric(pROC::ci.auc(roc_cox, method = "delong"))

# Harrell's C (time-to-event; larger score = shorter time -> reverse = TRUE)
cc  <- concordance(Surv(surv_years, event) ~ cox_score, data = dat_val, reverse = TRUE)
c_h <- cc$concordance; c_se <- sqrt(cc$var)

# Harrell's C of the full SANT Cox linear predictor (exact log-HRs, no point
# rounding). External counterpart of the SANT Step 10 C-statistic.
dat_val$lp_model <- with(dat_val,
                         age_50_59 * b_sant[["age_50_59"]] + age_60_69 * b_sant[["age_60_69"]] +
                           age_70plus * b_sant[["age_70plus"]] + male * b_sant[["male"]] +
                           biliary_ap * b_sant[["biliary_ap"]] + alcohol_ap * b_sant[["alcohol_ap"]] +
                           unspecified_ap * b_sant[["unspecified_ap"]] + other_ap * b_sant[["other_ap"]])
cc_m  <- concordance(Surv(surv_years, event) ~ lp_model, data = dat_val, reverse = TRUE)
c_m   <- cc_m$concordance; c_m_se <- sqrt(cc_m$var)

discrimination <- tibble(
  Dataset = COHORT, Outcome = "Cancer within 1 year", Census_date = format(CENSUS_DATE),
  SANT_optimism_corrected_C_model = corrected_C,
  SANT_optimism_corrected_C_score = corrected_C_score,
  N = nrow(dat_val), Cancer_events = sum(dat_val$event),
  AUROC = round(auc_ci[2], 3), AUROC_lo = round(auc_ci[1], 3), AUROC_hi = round(auc_ci[3], 3),
  Harrell_C_model = round(c_m, 3),
  Harrell_C_model_lo = round(c_m - 1.96 * c_m_se, 3), Harrell_C_model_hi = round(c_m + 1.96 * c_m_se, 3),
  Harrell_C_score = round(c_h, 3),
  Harrell_C_score_lo = round(c_h - 1.96 * c_se, 3), Harrell_C_score_hi = round(c_h + 1.96 * c_se, 3)
)
cat(sprintf("Harrell's C, full SANT Cox model (compare to SANT Step 10): %.3f (95%% CI %.3f-%.3f)\n",
            c_m, c_m - 1.96 * c_m_se, c_m + 1.96 * c_m_se))
cat(sprintf("AUROC (integer score):  %.3f (95%% CI %.3f-%.3f, DeLong)\n", auc_ci[2], auc_ci[1], auc_ci[3]))
cat(sprintf("Harrell's C:            %.3f (95%% CI %.3f-%.3f)\n\n",
            c_h, c_h - 1.96 * c_se, c_h + 1.96 * c_se))
write.csv(discrimination, out(sprintf("Cox_discrimination_%s.csv", COHORT)), row.names = FALSE)
auroc_cox <- auc_ci[2]; auroc_cox_ci <- auc_ci   # names used in the final summary
cat(sprintf("\n=== COX 1-YEAR HARRELL'S C: SANT (optimism-corrected) vs %s (external) ===\n", COHORT))
cox_c_table <- tibble(
  Predictor = c("Cox linear predictor (exact log-HRs)", "PAnCReAS integer score"),
  `SANT optimism-corrected C` = sprintf("%.3f", c(corrected_C, corrected_C_score)),
  `External C (95% CI)`       = sprintf("%.3f (%.3f\u2013%.3f)", c(c_m, c_h),
                                        c(c_m, c_h) - 1.96 * c(c_m_se, c_se),
                                        c(c_m, c_h) + 1.96 * c(c_m_se, c_se)),
  Difference = sprintf("%.3f", c(corrected_C - c_m, corrected_C_score - c_h))
)
names(cox_c_table)[3] <- sprintf("%s external C (95%% CI)", COHORT)
print(as.data.frame(cox_c_table), row.names = FALSE)
write.csv(cox_c_table, out(sprintf("Cox_Harrell_C_SANT_vs_%s_1yr.csv", COHORT)),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Note: the %s C is an external estimate - no optimism correction is needed or applied;\n", COHORT))
cat("      the SANT value is the bootstrap optimism-corrected internal estimate.\n\n")

# ── 11B. KM cumulative incidence by band (3 / 6 / 9 / 12 months) ─────────────
cat("=== STEP 11B: KM cumulative incidence by band ===\n\n")
horizons <- c("3 months" = 3/12, "6 months" = 6/12, "9 months" = 9/12, "1 year" = 1)

km_by_band <- imap_dfr(horizons, function(h, lbl) {
  dat_val %>%
    group_by(risk_band) %>%
    group_modify(~ bind_cols(tibble(N_total = nrow(.x)), km_at(.x, h))) %>%
    ungroup() %>%
    mutate(Timepoint = factor(lbl, levels = names(horizons)))
}) %>%
  mutate(across(c(km_pct, km_lo95, km_hi95), ~ round(.x, 2)))

print(km_by_band %>% arrange(risk_band, Timepoint))
write_release(km_by_band,
              list(n_events = c("km_pct", "km_lo95", "km_hi95"),
                   n_at_risk = character(0)),
              sprintf("Cox_KM_by_band_%s.csv", COHORT))

band_summary_1yr <- km_by_band %>%
  filter(Timepoint == "1 year") %>%
  mutate(NNS = ifelse(km_pct > 0, round(100 / km_pct), NA)) %>%
  select(risk_band, N_total, n_events, km_pct, km_lo95, km_hi95, NNS)
cat("\n1-year cancer rate and NNS by band:\n")
print(band_summary_1yr)
write_release(band_summary_1yr,
              list(n_events = c("km_pct", "km_lo95", "km_hi95", "NNS")),
              sprintf("Cox_band_summary_1yr_%s.csv", COHORT))

# ── 11C. Calibration (TRIPOD: SANT-predicted vs observed) ────────────────────
cat("\n=== STEP 11C: Calibration ===\n\n")

calib_by_band <- map_dfr(c(1), function(t) {
  pcol <- "pred_1yr"
  bind_rows(
    dat_val %>% group_by(risk_band) %>%
      group_modify(~ bind_cols(tibble(N = nrow(.x), Predicted_pct = 100 * mean(.x[[pcol]])),
                               km_at(.x, t))) %>%
      ungroup() %>% mutate(risk_band = as.character(risk_band)),
    bind_cols(tibble(risk_band = "Overall", N = nrow(dat_val),
                     Predicted_pct = 100 * mean(dat_val[[pcol]])),
              km_at(dat_val, t))
  ) %>% mutate(Horizon_yr = t)
}) %>%
  transmute(Horizon_yr, risk_band, N, n_events,
            Predicted_pct = round(Predicted_pct, 2),
            Observed_KM_pct = round(km_pct, 2),
            Observed_lo95 = round(km_lo95, 2), Observed_hi95 = round(km_hi95, 2),
            O_E_ratio = round(km_pct / Predicted_pct, 2))

print(calib_by_band)
write_release(calib_by_band,
              list(n_events = c("Observed_KM_pct", "Observed_lo95", "Observed_hi95", "O_E_ratio")),
              sprintf("Cox_calibration_by_band_%s.csv", COHORT))

dat_val$lp_sant <- dat_val$cox_score * ref_unit
slope_fit <- coxph(Surv(surv_years, event) ~ lp_sant, data = dat_val)
cal_slope <- coef(slope_fit)[[1]]
cal_ci    <- confint(slope_fit)[1, ]
cat(sprintf("\nCalibration slope: %.3f (95%% CI %.3f-%.3f)  [1 = ideal; <1 = predictions too extreme]\n",
            cal_slope, cal_ci[1], cal_ci[2]))
oe_1 <- calib_by_band$O_E_ratio[calib_by_band$Horizon_yr == 1 & calib_by_band$risk_band == "Overall"]
cat(sprintf("Overall O/E ratio at 1 year = %.2f  [1 = ideal]\n\n", oe_1))

calib_summary <- tibble(
  Metric = c("Calibration slope", "Calibration slope 95% CI lower", "Calibration slope 95% CI upper",
             "O/E ratio at 1 year"),
  Value  = round(c(cal_slope, cal_ci[1], cal_ci[2], oe_1), 3)
)
write.csv(calib_summary, out(sprintf("Cox_calibration_summary_%s.csv", COHORT)), row.names = FALSE)

# ── 11D. Sensitivity / specificity by cut-off ────────────────────────────────
cat("=== STEP 11D: Sensitivity / specificity by cutoff ===\n\n")

thresh_df <- map_dfr(sort(unique(dat_val$cox_score)), function(th) {
  pos <- dat_val$cox_score >= th; ev <- dat_val$event == 1
  tp <- sum(pos & ev); fp <- sum(pos & !ev); fn <- sum(!pos & ev); tn <- sum(!pos & !ev)
  tibble(Cutoff = th,
         Sensitivity = round(tp / (tp + fn), 3), Specificity = round(tn / (tn + fp), 3),
         PPV = round(ifelse(tp + fp == 0, NA, tp / (tp + fp)), 4),
         NPV = round(ifelse(tn + fn == 0, NA, tn / (tn + fn)), 4),
         TP = tp, FP = fp, FN = fn, TN = tn,
         N_flagged = tp + fp,
         Pct_flagged = round(100 * (tp + fp) / nrow(dat_val), 1))
})
cat("At the band cut-offs:\n")
print(thresh_df %>% filter(Cutoff %in% c(BAND_LOW_MAX + 1, BAND_HIGH_MIN)))
# Rows are cumulative: adjacent rows can be differenced to recover per-score
# counts, so check the release copy by hand before sending.
write_release(thresh_df,
              list(TP = c("Sensitivity", "PPV", "NPV", "FN"),
                   FN = c("Sensitivity", "PPV", "NPV", "TP"),
                   FP = c("Specificity", "PPV", "TN", "N_flagged", "Pct_flagged"),
                   TN = c("Specificity", "NPV", "FP"),
                   N_flagged = c("Pct_flagged")),
              sprintf("Cox_cutoff_performance_%s.csv", COHORT))
cat("\n")

# ── 11E. Imaging strategies: High only vs High + Moderate (Wilson 95% CI) ───
cat("=== STEP 11E: Imaging strategy performance ===\n\n")
fmt_ci <- function(k, n, dp = 1) {
  if (n == 0) return(NA_character_)
  ci <- wilson_ci(k, n)
  sprintf(paste0("%.", dp, "f%% (%.", dp, "f–%.", dp, "f)"),
          100 * k / n, 100 * ci[["lo"]], 100 * ci[["hi"]])
}
imaging_strategies <- map_dfr(
  list(c("High band only", BAND_HIGH_MIN),
       c("High + Moderate bands", BAND_LOW_MAX + 1)),
  function(s) {
    th <- as.numeric(s[2]); pos <- dat_val$cox_score >= th; ev <- dat_val$event == 1
    tp <- sum(pos & ev); fp <- sum(pos & !ev); fn <- sum(!pos & ev); tn <- sum(!pos & !ev)
    tibble(Strategy = s[1], Score_cutoff = paste0("≥", th),
           Imaged_n = tp + fp, Imaged_pct = round(100 * (tp + fp) / nrow(dat_val), 1),
           Cancers_detected = tp, Cancers_missed = fn,
           Sensitivity = fmt_ci(tp, tp + fn), Specificity = fmt_ci(tn, tn + fp),
           PPV = fmt_ci(tp, tp + fp), NPV = fmt_ci(tn, tn + fn, dp = 2),
           Images_per_cancer = ifelse(tp > 0, round((tp + fp) / tp), NA),
           FP = fp, TN = tn)
  })
print(as.data.frame(imaging_strategies %>% select(-FP, -TN)), row.names = FALSE)
write_release(imaging_strategies,
              list(Cancers_detected = c("Sensitivity", "PPV", "NPV", "Images_per_cancer", "Cancers_missed"),
                   Cancers_missed   = c("Sensitivity", "PPV", "NPV", "Cancers_detected", "Images_per_cancer"),
                   FP               = c("Specificity", "PPV", "Imaged_n", "Imaged_pct", "Images_per_cancer"),
                   TN               = c("Specificity", "NPV"),
                   Imaged_n         = c("Imaged_pct", "PPV", "Images_per_cancer")),
              sprintf("Cox_imaging_strategies_%s.csv", COHORT),
              drop_cols = c("FP", "TN"))
cat("\n")

# ── 11F. Band performance table (manuscript format) ─────────────────────────
band_performance <- function(d, cohort = "Cohort",
                             low_max = BAND_LOW_MAX, high_min = BAND_HIGH_MIN,
                             horizon_yrs = 1, out_dir = NULL, score_min = SCORE_MIN) {
  
  stopifnot(all(c("surv_years", "event", "cox_score") %in% names(d)))
  if (min(d$cox_score, na.rm = TRUE) < score_min)
    warning(sprintf("cox_score has values below the SANT achievable minimum (%d) - check the scoring.",
                    score_min))
  
  lbl <- c(sprintf("Low (≤%d)", low_max),
           sprintf("Moderate (%d–%d)", low_max + 1, high_min - 1),
           sprintf("High (≥%d)", high_min))
  d <- d %>%
    mutate(band_num = case_when(cox_score <= low_max  ~ 1L,
                                cox_score >= high_min ~ 3L,
                                TRUE                  ~ 2L),
           band = factor(lbl[band_num], levels = lbl))
  
  # ---- per-band observed rate + 1-year KM ----------------------------------
  band_tab <- lapply(lbl, function(b) {
    x <- d[d$band == b, ]
    n <- nrow(x); k <- sum(x$event)
    ci <- if (n > 0) binom.test(k, n)$conf.int else c(NA, NA)
    if (n > 0 && k > 0) {
      s <- summary(survfit(Surv(surv_years, event) ~ 1, data = x),
                   times = horizon_yrs, extend = TRUE)
      km <- c(1 - s$surv, 1 - s$upper, 1 - s$lower) * 100
      n_risk <- s$n.risk
    } else {
      km <- c(0, NA, NA); n_risk <- if (n > 0) sum(x$surv_years >= horizon_yrs) else NA
    }
    tibble(Cohort = cohort, Band = b, N = n,
           Pct_of_cohort = 100 * n / nrow(d),
           Cancers = k,
           Observed_pct = 100 * k / max(n, 1),
           Observed_lo95 = 100 * ci[1], Observed_hi95 = 100 * ci[2],
           KM_1yr_pct = km[1], KM_1yr_lo95 = km[2], KM_1yr_hi95 = km[3],
           N_at_risk_1yr = n_risk)
  }) %>% bind_rows()
  
  # ---- C-statistics (Harrell's; larger score = shorter time) ---------------
  c_row <- function(x, what) {
    cc <- concordance(Surv(d$surv_years, d$event) ~ x, reverse = TRUE)
    se <- sqrt(cc$var)
    tibble(Cohort = cohort, Predictor = what,
           C = cc$concordance,
           C_lo95 = cc$concordance - 1.96 * se,
           C_hi95 = cc$concordance + 1.96 * se)
  }
  c_tab <- bind_rows(c_row(d$band_num,  "PAnCReAS band (3 levels)"),
                     c_row(d$cox_score, "PAnCReAS integer score"))
  
  # ---- formatted manuscript table (<10 rule applied) ------------------------
  f1 <- function(x) sprintf("%.2f", x)
  manuscript <- band_tab %>%
    mutate(supp = sdc_small(Cancers)) %>%
    transmute(Cohort, Band,
              `N (%)` = ifelse(sdc_small(N), "<10",
                               sprintf("%s (%.1f%%)", format(N, big.mark = ",", trim = TRUE), Pct_of_cohort)),
              Cancers = ifelse(supp, "<10", as.character(Cancers)),
              `Observed cancer rate, % (95% CI)` = ifelse(supp, "S",
                                                          sprintf("%s (%s–%s)", f1(Observed_pct), f1(Observed_lo95), f1(Observed_hi95))),
              `1-year KM incidence, % (95% CI)` = ifelse(supp, "S",
                                                         ifelse(is.na(KM_1yr_lo95), sprintf("%s", f1(KM_1yr_pct)),
                                                                sprintf("%s (%s–%s)", f1(KM_1yr_pct), f1(KM_1yr_lo95), f1(KM_1yr_hi95)))))
  c_text <- sprintf("%s: Harrell's C %.3f (95%% CI %.3f–%.3f)",
                    c_tab$Predictor, c_tab$C, c_tab$C_lo95, c_tab$C_hi95)
  
  cat(sprintf("\n=== PAnCReAS band performance - %s (N = %s, %d cancers) ===\n",
              cohort, format(nrow(d), big.mark = ","), sum(d$event)))
  print(as.data.frame(manuscript), row.names = FALSE)
  cat(paste0("  ", c_text, collapse = "\n"), "\n\n")
  
  if (!is.null(out_dir)) {
    tag <- gsub("[^A-Za-z0-9]", "", cohort)
    write_release(band_tab,
                  list(N = c("Pct_of_cohort"),
                       Cancers = c("Observed_pct", "Observed_lo95", "Observed_hi95",
                                   "KM_1yr_pct", "KM_1yr_lo95", "KM_1yr_hi95")),
                  sprintf("Band_performance_%s.csv", tag))
    write.csv(c_tab,      file.path(out_dir, sprintf("C_statistic_%s.csv", tag)), row.names = FALSE)
    write.csv(manuscript, file.path(out_dir, sprintf("Band_table_manuscript_%s.csv", tag)),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  invisible(list(band_table = band_tab, c_statistic = c_tab,
                 manuscript = manuscript, c_text = c_text))
}

bp_val <- band_performance(dat_val, cohort = COHORT, out_dir = out_dir)

# =============================================================================
# 12. FIGURES
# =============================================================================

cat("\n=== STEP 12: Figures ===\n\n")

theme_paper <- theme_bw(base_size = 13) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"),
        plot.caption  = element_text(colour = "grey40", size = 9),
        legend.position = "bottom")

# ── Analysis A: ROC ─────────────────────────────────────────────────────────
roc_df_A <- data.frame(fpr = 1 - roc_A$specificities, tpr = roc_A$sensitivities)
p_roc_A  <- ggplot(roc_df_A, aes(x=fpr, y=tpr)) +
  geom_line(colour="#2166AC", linewidth=1.3) +
  geom_abline(linetype="dashed", colour="grey50") +
  annotate("text", x=0.45, y=0.15,
           label = sprintf(
             "AUROC = %.3f (%.3f\u2013%.3f)\n\nSANT Platt applied directly\n(no %s outcome data used)",
             auroc_A, auroc_A_ci[1], auroc_A_ci[3], COHORT),
           hjust=0, size=3.2, colour="grey30") +
  labs(
    title    = sprintf("%s \u2014 Blind External Validation ROC (%s, 1-year)", ML_NAME, COHORT),
    subtitle = sprintf("N = %s  |  %s cancer events  |  SANT Platt applied directly",
                       format(nrow(dat_val), big.mark = ","), sdc_lab(sum(y_val))),
    x = "1 - Specificity", y = "Sensitivity",
    caption  = sprintf("No %s outcome labels used. Census date: %s.", COHORT, format(CENSUS_DATE))
  ) + theme_paper
ggsave(file.path(out_dir, "AnalysisA_roc.png"), p_roc_A, width=8, height=7, dpi=300)
cat("Saved: AnalysisA_roc.png\n")

# ── Analysis A: Calibration ──────────────────────────────────────────────────
calib_A <- tibble(prob=prob_A, label=y_val) %>%
  mutate(bin = quantile_bins(prob, 10)) %>%
  group_by(bin) %>%
  summarise(mean_pred=mean(prob), mean_obs=mean(label), n=n(),
            n_events = sum(label), .groups="drop") %>%
  filter(!is.na(bin))
warn_small_bins(calib_A$n_events, "AnalysisA_calibration.png")
p_cal_A <- ggplot(calib_A, aes(x=mean_pred, y=mean_obs)) +
  geom_abline(linetype="dashed", colour="grey50", linewidth=1) +
  geom_line(colour="#2166AC", linewidth=1.2) +
  geom_point(aes(size=n), colour="#2166AC", alpha=0.85) +
  scale_size_continuous(name="N in bin", range=c(2,8)) +
  annotate("text", x=min(calib_A$mean_pred)*1.1,
           y=max(calib_A$mean_obs, na.rm=TRUE)*0.9,
           label=sprintf("Brier = %.4f\nSkill = %.4f", brier_A, brier_skA),
           hjust=0, size=3.5, colour="grey30") +
  labs(
    title   = sprintf("%s \u2014 Calibration Curve (SANT Platt, no %s labels, 1-year)", ML_NAME, COHORT),
    subtitle = "SANT-derived Platt calibrator applied without modification",
    x = "Mean Predicted Probability (quantile bins)", y = "Observed Event Rate",
    caption = "Dashed = perfect calibration."
  ) + theme_paper
ggsave(file.path(out_dir, "AnalysisA_calibration.png"), p_cal_A, width=8, height=7, dpi=300)
cat("Saved: AnalysisA_calibration.png\n")

# ── PAnCReAS 1-year figures ─────────────────────────────────────────────────
cat("\n--- PAnCReAS figures ---\n")
theme_cox <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey40"),
        plot.caption = element_text(colour = "grey40", size = 9),
        legend.position = "bottom", panel.grid.minor = element_blank())
sub_txt <- sprintf("%s  |  N = %s  |  %s cancers within 1 year",
                   COHORT, format(nrow(dat_val), big.mark = ","), sdc_lab(sum(dat_val$event)))

# Fig A: 1-year rate by band with NNS (bands with 1-9 cancers: bar and rate hidden)
band_plot <- band_summary_1yr %>%
  mutate(supp  = sdc_small(n_events),
         y     = ifelse(supp, NA_real_, km_pct),
         lo    = ifelse(supp, NA_real_, km_lo95),
         hi    = ifelse(supp, NA_real_, km_hi95),
         lab_y = ifelse(supp, 0, km_pct),
         lab   = ifelse(supp, sprintf("<10 cancers\nN=%s", format(N_total, big.mark = ",")),
                        sprintf("%.2f%%\n(NNS=%s)\nN=%s", km_pct, ifelse(is.na(NNS), "-", NNS),
                                format(N_total, big.mark = ","))))
p_band <- ggplot(band_plot, aes(risk_band, y, fill = risk_band)) +
  geom_col(width = 0.6, colour = "grey20", linewidth = 0.4, na.rm = TRUE) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, colour = "grey20", na.rm = TRUE) +
  geom_text(aes(y = lab_y, label = lab), vjust = -0.3, size = 3.4, fontface = "bold") +
  scale_fill_manual(values = BAND_COLOURS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = sprintf("PAnCReAS - 1-year cancer incidence by band (%s)", COHORT), subtitle = sub_txt,
       x = "PAnCReAS band", y = "1-year KM cumulative incidence (%)",
       caption = "Error bars = 95% CI. NNS = number needed to scan to detect one cancer. Bands with 1-9 cancers suppressed.") +
  theme_cox
ggsave(out(sprintf("Cox_band_rate_1yr_%s.png", COHORT)), p_band, width = 9, height = 7, dpi = 400)

# Fig B: KM incidence by band over time (matches SANT Figure 1)
warn_small_bins(km_by_band$n_events, sprintf("Cox_KM_by_band_over_time_%s.png", COHORT))
p_km <- ggplot(km_by_band, aes(Timepoint, km_pct, colour = risk_band, group = risk_band)) +
  geom_line(linewidth = 1.1) + geom_point(size = 2.8) +
  geom_text(aes(label = ifelse(sdc_small(n_events), "<10 cases", sprintf("%.1f%%", km_pct))),
            vjust = -1, size = 3.2, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = BAND_COLOURS, name = "PAnCReAS band") +
  scale_y_continuous(limits = c(0, max(km_by_band$km_pct, na.rm = TRUE) * 1.15)) +
  labs(title = sprintf("PAnCReAS - KM cumulative incidence by band (%s)", COHORT), subtitle = sub_txt,
       x = "Time since AP", y = "KM cumulative incidence (%)") +
  theme_cox
ggsave(out(sprintf("Cox_KM_by_band_over_time_%s.png", COHORT)), p_km, width = 8, height = 5.5, dpi = 400)

# Fig C: ROC
roc_df <- tibble(fpr = 1 - roc_cox$specificities, tpr = roc_cox$sensitivities)
p_roc <- ggplot(roc_df, aes(fpr, tpr)) +
  geom_line(colour = "#2166AC", linewidth = 1.3) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.45, y = 0.15, hjust = 0, size = 3.5, colour = "grey30",
           label = sprintf("AUROC = %.3f (%.3f\u2013%.3f)\nHarrell's C = %.3f",
                           auc_ci[2], auc_ci[1], auc_ci[3], c_h)) +
  labs(title = sprintf("PAnCReAS - ROC curve (%s)", COHORT), subtitle = sub_txt,
       x = "1 - Specificity", y = "Sensitivity") +
  theme_cox
ggsave(out(sprintf("Cox_roc_%s.png", COHORT)), p_roc, width = 8, height = 7, dpi = 400)

# Fig D: Calibration - predicted vs observed by band at 1 year
cal_plot_dat <- calib_by_band %>%
  filter(risk_band != "Overall") %>%
  mutate(risk_band = factor(risk_band, levels = unname(BAND_LABELS)),
         Horizon = factor(paste0(Horizon_yr, "-year")))
warn_small_bins(cal_plot_dat$n_events, sprintf("Cox_calibration_%s.png", COHORT))
lim <- max(c(cal_plot_dat$Predicted_pct, cal_plot_dat$Observed_hi95), na.rm = TRUE) * 1.1
p_cal <- ggplot(cal_plot_dat, aes(Predicted_pct, Observed_KM_pct, colour = risk_band)) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = Observed_lo95, ymax = Observed_hi95), width = 0) +
  geom_point(size = 3.2) +
  facet_wrap(~ Horizon) +
  scale_colour_manual(values = BAND_COLOURS, name = "PAnCReAS band") +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
  labs(title = sprintf("PAnCReAS - calibration in %s", COHORT), subtitle = sub_txt,
       x = "SANT-predicted risk (%)", y = sprintf("Observed KM risk in %s (%%)", COHORT),
       caption = sprintf("Calibration slope %.2f (%.2f\u2013%.2f). O/E at 1 yr %.2f. Dashed = perfect calibration.",
                         cal_slope, cal_ci[1], cal_ci[2], oe_1)) +
  theme_cox
ggsave(out(sprintf("Cox_calibration_%s.png", COHORT)), p_cal, width = 7, height = 6, dpi = 400)

cat("Saved 4 PAnCReAS figures.\n\n")

# ── PAnCReAS decision curve analysis (time-to-event, 1-year risk) ───────────
# Net benefit uses KM risk at 1 year, so deaths before 1 year are handled
# as censored rather than counted as cancer-free.
{
  cat("--- PAnCReAS decision curve analysis (time-to-event) ---\n\n")
  dca_dat <- dat_val %>% select(surv_years, event, pred_1yr, clinical_rule)
  
  dca_result <- dcurves::dca(
    Surv(surv_years, event) ~ pred_1yr + clinical_rule,
    data       = dca_dat,
    time       = HORIZON_YRS,
    thresholds = seq(0.005, 0.10, by = 0.005),
    label      = list(pred_1yr      = "PAnCReAS (integer score, calibrated)",
                      clinical_rule = "Clinical rule")
  )
  
  strategy_labels <- c(pred_1yr = "PAnCReAS (integer score, calibrated)",
                       clinical_rule = "Clinical rule", all = "Scan-All", none = "Scan-None")
  dca_table <- as_tibble(dca_result$dca) %>%
    filter(variable %in% names(strategy_labels),
           sapply(threshold, function(t) any(abs(t - c(0.01, 0.02, 0.03, 0.04, 0.05)) < 1e-8))) %>%
    mutate(Strategy = factor(unname(strategy_labels[as.character(variable)]),
                             levels = unname(strategy_labels)),
           thresh_lbl = paste0("NB_at_", round(threshold * 100), "pct")) %>%
    select(Strategy, thresh_lbl, net_benefit) %>%
    pivot_wider(names_from = thresh_lbl, values_from = net_benefit) %>%
    arrange(Strategy)
  print(dca_table)
  write.csv(dca_table, out(sprintf("Cox_dca_summary_table_%s.csv", COHORT)), row.names = FALSE)
  
  p_dca <- plot(dca_result, smooth = TRUE) +
    labs(title = sprintf("PAnCReAS - decision curve analysis (%s)", COHORT), subtitle = sub_txt,
         x = "Threshold probability (1-year risk)", y = "Net benefit") + theme_cox
  ggsave(out(sprintf("Cox_dca_%s.png", COHORT)), p_dca, width = 10, height = 7, dpi = 300)
  cat("\n")
}

# =============================================================================
# 13. DEMOGRAPHICS (Table 1)
# =============================================================================

cat(sprintf("\n=== STEP 13: Table 1 demographics (%s) ===\n\n", COHORT))

t1_n_pct <- function(x, denom) {
  k <- sum(x, na.rm = TRUE)
  if (k > 0 && k < 10) return("<10")          # data custodian small-cell rule
  sprintf("%s (%.1f%%)", format(k, big.mark = ","), 100 * k / denom)
}
t1_med_iqr <- function(x, dp = 1) {
  if (sdc_small(sum(!is.na(x)))) return("S")  # median of fewer than 10 people
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", dp, "f (%.", dp, "f–%.", dp, "f)"), q[2], q[1], q[3])
}

table1_col <- function(d) {
  n   <- nrow(d)
  age <- d$nv_k85_firstage
  c(
    "N"                                   = sdc_lab(n),
    "Cancer events within 1 year, n (%)"  = t1_n_pct(d$event == 1, n),
    "Cancer events, any time, n (%)"      = t1_n_pct(d$event_full == 1, n),
    "Median age at first AP, years (IQR)" = t1_med_iqr(age),
    "Age <50 years, n (%)"                = t1_n_pct(age < 50, n),
    "Age 50–64 years, n (%)"         = t1_n_pct(age >= 50 & age < 65, n),
    "Age ≥65 years, n (%)"           = t1_n_pct(age >= 65, n),
    "Male, n (%)"                         = t1_n_pct(d$male == 1, n),
    "Female, n (%)"                       = t1_n_pct(d$male == 0, n),
    "Biliary AP, n (%)"                   = t1_n_pct(d$biliary_ap == 1, n),
    "Alcohol-induced AP, n (%)"           = t1_n_pct(d$alcohol_ap == 1, n),
    "Unspecified AP, n (%)"               = t1_n_pct(d$unspecified_ap == 1, n),
    "Other AP, n (%)"                     = t1_n_pct(d$other_ap == 1, n)
  )
}

ca_1yr <- dat_val %>% filter(event == 1)
ca_all <- dat_val %>% filter(event_full == 1)

c_ap   <- table1_col(dat_val)
table1 <- tibble(
  Characteristic   = names(c_ap),
  AP               = unname(c_ap),
  Cancer_1yr       = unname(table1_col(ca_1yr)),
  Cancer_any_time  = unname(table1_col(ca_all))
)
names(table1)[2:4] <- paste0(COHORT, c("_AP", "_cancer_1yr", "_cancer_any_time"))

print(as.data.frame(table1), row.names = FALSE)
write.csv(table1, out(sprintf("Table1_demographics_%s.csv", COHORT)),
          row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("\nSaved: Table1_demographics_%s.csv\n", COHORT))
cat("       Percentages use each column's own N. Counts 1-9 shown as '<10';\n")
cat("       check secondary suppression ('S') by hand before release.\n")
cat("       Aetiology categories are not mutually exclusive if patients have >1 code.\n\n")

# =============================================================================
# 14. FIGURE S3 PANEL — ML MODEL CALIBRATION
# =============================================================================

cat(sprintf("\n=== STEP 14: Figure S3 (C) calibration panel + %s calibration statistics ===\n\n", ML_NAME))

FIGS3_COLOURS <- c(grey = "#8C8C8C", red = "#D64541", blue = "#3D6DB5")

figS3_calibration <- function(prob, y, colour, file) {
  cd <- tibble(prob = prob, label = y) %>%
    mutate(bin = quantile_bins(prob, 10)) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(prob), mean_obs = mean(label), n = n(),
              n_events = sum(label), .groups = "drop") %>%
    filter(!is.na(bin))
  warn_small_bins(cd$n_events, file)
  p <- ggplot(cd, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(linetype = "dashed", colour = "grey50", linewidth = 1) +
    geom_line(colour = colour, linewidth = 1) +
    geom_point(aes(size = n), colour = colour) +
    scale_size_continuous(name = "N in bin", range = c(2, 8)) +
    labs(x = "Mean Predicted Probability (quantile bins)", y = "Observed Event Rate") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(out(file), p, width = 7, height = 6, dpi = 300)
  cat("Saved:", file, "\n")
  invisible(cd)
}

ml_calibration_stats <- function(prob, y, label) {
  lp    <- qlogis(pmin(pmax(prob, 1e-6), 1 - 1e-6))
  f_s   <- glm(y ~ lp, family = binomial)
  f_i   <- glm(y ~ 1 + offset(lp), family = binomial)
  ci_s  <- confint.default(f_s)["lp", ]
  ci_i  <- confint.default(f_i)[1, ]
  tibble(Model = ML_NAME, Probabilities = label, N = length(y), Events = sum(y),
         Calibration_slope = round(coef(f_s)[["lp"]], 3),
         Slope_lo95 = round(ci_s[1], 3), Slope_hi95 = round(ci_s[2], 3),
         CITL = round(coef(f_i)[1], 3),
         CITL_lo95 = round(ci_i[1], 3), CITL_hi95 = round(ci_i[2], 3),
         O_E_ratio = round(sum(y) / sum(prob), 3),
         Brier = round(mean((prob - y)^2), 4))
}

figS3_calibration(prob_A, y_val, FIGS3_COLOURS[["blue"]],
                  sprintf("FigS3_C-i_%s_SANT_calibration.png", COHORT))

ml_cal_stats <- ml_calibration_stats(prob_A, y_val, "SANT Platt (blind), 1-year outcome")
print(as.data.frame(ml_cal_stats), row.names = FALSE)
write_release(ml_cal_stats, list(Events = c("O_E_ratio")),
              sprintf("ML_calibration_statistics_%s.csv", COHORT))
cat(sprintf("Saved: ML_calibration_statistics_%s.csv\n", COHORT))
cat("Slope: 1 = ideal (<1 predictions too extreme). CITL: 0 = ideal (>0 under-prediction).\n")

# =============================================================================
# 15. FINAL SUMMARY (runs last, after every output has been written)
# =============================================================================

cat("\n")
cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("  FINAL SUMMARY — %s EXTERNAL VALIDATION\n", COHORT))
cat(sprintf("  N = %d  |  %d cancer events (%.2f%%)\n",
            nrow(dat_val), sum(y_val), mean(y_val)*100))
cat(sprintf("  Census date: %s (hardcoded) | first AP on/before %s\n",
            format(CENSUS_DATE), format(FU_CUTOFF)))
cat("════════════════════════════════════════════════════════════\n\n")

cat("Flowchart:\n")
print(as.data.frame(flowchart), row.names = FALSE)

cat(sprintf("\n%s — BLIND EXTERNAL VALIDATION:\n", toupper(ML_NAME)))
cat(sprintf("  AUROC (raw %s):     %.4f (95%% CI: %.3f-%.3f)\n",
            ML_NAME, auroc_raw, auroc_raw_ci[1], auroc_raw_ci[3]))
cat(sprintf("  AUROC (SANT Platt):      %.4f (95%% CI: %.3f-%.3f)\n",
            auroc_A, auroc_A_ci[1], auroc_A_ci[3]))
cat(sprintf("  Sensitivity at %.4f:   %.4f (%d / %d cancers caught)\n",
            THRESHOLD_A, perf_A$sens, perf_A$tp, perf_A$tp + perf_A$fn))
cat(sprintf("  Specificity at %.4f:   %.4f\n\n", THRESHOLD_A, perf_A$spec))

cat("PAnCReAS COX SCORE (1-YEAR) — RISK STRATIFICATION:\n")
cat(sprintf("  Harrell's C (full SANT model): %.3f (%.3f-%.3f)\n",
            c_m, c_m - 1.96 * c_m_se, c_m + 1.96 * c_m_se))
cat(sprintf("  Harrell's C (integer score):   %.3f (%.3f-%.3f)\n", c_h, c_h - 1.96 * c_se, c_h + 1.96 * c_se))
cat(sprintf("  AUROC (integer score):         %.3f (%.3f-%.3f)\n", auc_ci[2], auc_ci[1], auc_ci[3]))
cat(sprintf("  Calibration slope:             %.3f (%.3f-%.3f)\n", cal_slope, cal_ci[1], cal_ci[2]))
cat(sprintf("  O/E ratio at 1 year:           %.2f\n", oe_1))
cat(sprintf("\n  Harrell's C, SANT optimism-corrected vs %s external:\n", COHORT))
print(as.data.frame(cox_c_table), row.names = FALSE)
cat(paste0("  ", bp_val$c_text, collapse = "\n"), "\n")
cat("\n  Band performance (manuscript format, <10 rule applied):\n")
print(as.data.frame(bp_val$manuscript), row.names = FALSE)
cat("\n  Imaging strategies (unsuppressed - console only):\n")
print(as.data.frame(imaging_strategies %>% select(-FP, -TN)), row.names = FALSE)
cat("\n  1-year incidence by band (unsuppressed - console only):\n")
print(band_summary_1yr)
cat("\n")

cat("════════════════════════════════════════════════════════════\n")
cat(sprintf("  KEY OUTPUT FILES — %s\n", out_dir))
cat(sprintf("  Cohort:  Flowchart_%s_1yr.csv, Table1_demographics_%s.csv\n", COHORT, COHORT))
cat(sprintf("  %s: Table_S6_%s_1yr.csv, AnalysisA_blind_validation.csv,\n", ML_NAME, COHORT))
cat(sprintf("           AnalysisA_roc.png, FigS3_C-i_%s_SANT_calibration.png,\n", COHORT))
cat(sprintf("           ML_calibration_statistics_%s.csv\n", COHORT))
cat(sprintf("  Cox:     Cox_Harrell_C_SANT_vs_%s_1yr.csv, Cox_discrimination_%s.csv,\n", COHORT, COHORT))
cat(sprintf("           Band_table_manuscript_%s.csv, Cox_KM_by_band_%s.csv,\n", COHORT, COHORT))
cat(sprintf("           Cox_calibration_*_%s.csv, Cox_dca_summary_table_%s.csv + figures\n", COHORT, COHORT))
cat("════════════════════════════════════════════════════════════\n")
cat("All outputs are aggregate - no patient-level data are saved.\n")
cat("Counts of 1-9 are shown as <10 and linked values as S in every released table.\n")
cat(sprintf("=== %s VALIDATION COMPLETE ===\n\n", COHORT))
cat(sprintf("Please send the %s folder back to the SANT research team.\n", out_dir))
cat(sprintf("Do NOT send %s (unsuppressed, internal only).\n", internal_dir))