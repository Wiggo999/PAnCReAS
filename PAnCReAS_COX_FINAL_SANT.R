# =============================================================================
# PAnCReAS Cox Model - FITTED ON 1-YEAR FOLLOW-UP
# =============================================================================

library(tidyverse)
library(survival)
library(glmnet)
library(pROC)
library(dcurves)
library(ggplot2)
library(rms)

# =============================================================================
# STEP 0 — OUTPUT FOLDER + HORIZON
# =============================================================================

OUT_DIR <- "~/Desktop/PAnCReAS_16-10"
OUT_DIR <- path.expand(OUT_DIR)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output folder ready:", OUT_DIR, "\n\n")

HORIZON_YRS <- 1   # single source of truth for the follow-up limit (years)


# =============================================================================
# STEP 0B — PAnCReAS RISK-BAND CUT-OFFS, 1-YEAR MODEL (single source of truth)
# Low <= 1 | Moderate 2-8 | High 9-12
# Set from the point-by-point 1-year risk table (risk_by_point_1yr.csv).
# =============================================================================

BAND_LOW_MAX  <- 1    # Low      : score <= 1
BAND_HIGH_MIN <- 9    # High     : score >= 9   (Moderate = 2 to 8)
SCORE_MAX     <- 12   # highest achievable score (used in the High label only)

BAND_LABELS <- c(
  Low      = sprintf("Low (\u2264%d)", BAND_LOW_MAX),
  Moderate = sprintf("Moderate (%d\u2013%d)", BAND_LOW_MAX + 1, BAND_HIGH_MIN - 1),
  High     = sprintf("High (%d\u2013%d)", BAND_HIGH_MIN, SCORE_MAX)
)
BAND_COLOURS <- c(Low = "#0D0887", Moderate = "#CC4778", High = "#F0A227")

assign_band <- function(score) {
  factor(case_when(
    score <= BAND_LOW_MAX  ~ "Low",
    score >= BAND_HIGH_MIN ~ "High",
    TRUE                   ~ "Moderate"
  ), levels = c("Low", "Moderate", "High"))
}



# =============================================================================
# STEP 0C — COHORT DEFINITION WITH EXCLUSIONS (single source of truth)
# Used for every cohort built in this script, so all outputs share the same
# patients. Exclusions, applied in this order:
# Starts from the WHOLE cohort (no nv_k85 filter); patients without an AP
# episode drop out at the missing first-AP date step.
#   1. K86 == 1 (nv_k86) - removed from the cohort
#      entirely (neither event nor control)
#   2. missing/invalid first-AP date
#   3. pancreatic cancer diagnosed BEFORE first AP (diffweeks < 0) - removed,
#      not just kept as a non-event
#   4. pancreatic cancer flagged but no diagnosis date (diffweeks or pc_date
#      missing) - cannot be placed in time, so neither event nor control.
#   5. first AP after CENSUS_DATE - 365 days: < 1 year of possible follow-up,
#      so "no cancer within 1 year" could not be verified
# =============================================================================

CENSUS_DATE <- as.Date("2022-01-01")
FU_CUTOFF   <- CENSUS_DATE - 365   # last first-AP date with a full year of follow-up

build_ap_cohort <- function(raw, census = CENSUS_DATE, horizon = HORIZON_YRS) {
  cutoff <- census - 365
  d0 <- raw %>%
    mutate(
      ap_date    = as.Date(ap_date_clean,    format = "%Y-%m-%d"),
      pc_date    = as.Date(pc_date_clean,    format = "%Y-%m-%d"),
      death_date = as.Date(death_date_clean, format = "%Y-%m-%d"),
      cp_flag    = coalesce(nv_k86 == 1, FALSE),
      ca_flag    = coalesce(k85_panctcancer == 1, FALSE),
      prior_ca   = ca_flag & coalesce(k85_panctcancer_diffweeks < 0, FALSE),
      ca_nodate  = ca_flag & !prior_ca &
        (is.na(k85_panctcancer_diffweeks) | is.na(pc_date))
    )
  d0b <- d0  %>% filter(!cp_flag)
  d1  <- d0b %>% filter(!is.na(ap_date))
  d2  <- d1  %>% filter(!prior_ca)
  d3  <- d2  %>% filter(!ca_nodate)
  d4  <- d3  %>% filter(ap_date <= cutoff)
  
  flow <- tibble(
    Step = c("Total cohort",
             "Excluded: K86 == 1",
             "Excluded: missing or invalid first-AP date",
             "Excluded: pancreatic cancer diagnosed before first AP",
             "Excluded: pancreatic cancer recorded with no diagnosis date",
             sprintf("Excluded: first AP after %s (< 1 year possible follow-up)",
                     format(cutoff)),
             "Cohort after exclusions"),
    N = c(nrow(d0), nrow(d0) - nrow(d0b), nrow(d0b) - nrow(d1),
          nrow(d1) - nrow(d2), nrow(d2) - nrow(d3), nrow(d3) - nrow(d4), nrow(d4))
  )
  
  out <- d4 %>%
    mutate(
      event_full = as.integer(ca_flag & !is.na(pc_date) &
                                k85_panctcancer_diffweeks >= 0),
      end_date = case_when(
        event_full == 1                            ~ pc_date,
        sant_deaths_flag == 1 & !is.na(death_date) ~ death_date,
        TRUE                                       ~ census),
      # pmax floors at 1 day for same-day/same-admission cancer diagnoses
      surv_years_full = pmax(as.numeric(end_date - ap_date) / 365.25, 1 / 365.25),
      event      = as.integer(event_full == 1 & surv_years_full <= horizon),
      surv_years = pmin(surv_years_full, horizon)
    )
  attr(out, "flow") <- flow
  out
}


# =============================================================================
# STEP 0D — SMALL-CELL SUPPRESSION (custodian rule: counts 1-9 -> "<10")
# Release copies go to OUT_DIR. Unsuppressed copies go to INTERNAL_DIR for
# your own checking only - do not request release of that folder.
# Linked cells that would let a suppressed count be back-calculated (rates,
# %, sens/spec, NNI, the complementary count) are set to "S". Differencing
# across rows (e.g. cumulative cut-off tables) still needs a check by hand.
# =============================================================================

SDC_MIN      <- 10
INTERNAL_DIR <- file.path(OUT_DIR, "internal_unsuppressed_DO_NOT_RELEASE")
dir.create(INTERNAL_DIR, showWarnings = FALSE, recursive = TRUE)

sdc_small <- function(x) {
  x <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
  !is.na(x) & x > 0 & x < SDC_MIN
}
sdc_lab <- function(x) ifelse(sdc_small(x), "<10", format(x, big.mark = ",", trim = TRUE))

# rules: named list, count column -> linked columns to blank with "S"
sdc <- function(df, rules, label = "") {
  out <- df %>% mutate(across(everything(), as.character))
  n_supp <- 0
  for (cc in names(rules)) {
    f <- sdc_small(df[[cc]])
    if (!any(f)) next
    out[[cc]][f] <- "<10"
    for (lc in intersect(rules[[cc]], names(out))) out[[lc]][f] <- "S"
    n_supp <- n_supp + sum(f)
  }
  if (n_supp > 0)
    cat(sprintf("  [SDC] %s: %d cell(s) suppressed - check secondary suppression by hand.\n",
                label, n_supp))
  out
}

write_release <- function(df, rules, file, drop_cols = character(0)) {
  write.csv(df, file.path(INTERNAL_DIR, file), row.names = FALSE, fileEncoding = "UTF-8")
  rel <- sdc(df %>% select(-any_of(drop_cols)), rules, file)
  write.csv(rel, file.path(OUT_DIR, file), row.names = FALSE, fileEncoding = "UTF-8")
  invisible(rel)
}


# =============================================================================
# STEP 0A — GUARD CLAUSE: VERIFY REQUIRED COLUMNS EXIST BEFORE STEP 1
# =============================================================================

required_cols <- c(
  "nv_k85", "ap_date_clean", "pc_date_clean", "death_date_clean",
  "k85_panctcancer", "k85_panctcancer_diffweeks", "nv_k86",
  "sant_deaths_flag", "nv_k85_firstage", "sant_nv_sex",
  "nv_k851", "nv_k852", "nv_k859", "nv_k85own", "nv_k85other", "nv_k858", "nv_k853"
)

missing_cols <- setdiff(required_cols, names(cmbspine_pancreatitis))

if (length(missing_cols) > 0) {
  near_matches <- lapply(missing_cols, function(col) {
    hits <- grep(gsub("_diffweeks$|_clean$", "", col),
                 names(cmbspine_pancreatitis), value = TRUE, ignore.case = TRUE)
    if (length(hits) == 0) "  (no similar column names found)" else
      paste0("  possible matches: ", paste(hits, collapse = ", "))
  })
  msg <- paste0(
    "cmbspine_pancreatitis is missing ", length(missing_cols),
    " required column(s):\n",
    paste0("  - ", missing_cols, "\n", near_matches, collapse = "\n"),
    "\n\nThis usually means cmbspine_pancreatitis in the current session is ",
    "stale, was reloaded without the cancer-linkage merge step, or has a ",
    "renamed/duplicated column. Re-check names(cmbspine_pancreatitis) and ",
    "grep(\"panctcancer\", names(cmbspine_pancreatitis), value = TRUE) ",
    "before re-running."
  )
  stop(msg)
}

cat("Guard check passed: all", length(required_cols),
    "required columns found in cmbspine_pancreatitis.\n\n")


# =============================================================================
# STEP 1 — BUILD SURVIVAL DATASET (TRUNCATED AT 1 YEAR)
# =============================================================================
# Entry criterion:  whole cohort with a first AP on/before FU_CUTOFF
#                   (CENSUS_DATE - 365 days): everyone can complete 1 year
# Exclusions:       Step 0C (chronic pancreatitis nv_k86 == 1, prior
#                   cancer, cancer with no date, < 1 year possible
#                   follow-up), then complete-case on predictors
# Event definition: k85_panctcancer == 1 AND diffweeks >= 0
#                   AND diagnosed within HORIZON_YRS of first AP
# Censoring:        death or HORIZON_YRS, whichever first (the census date
#                   can no longer fall inside the first year)
# =============================================================================

census_date <- CENSUS_DATE
cat("Census date (end of follow-up):", format(census_date), "\n")
cat("Follow-up limit (horizon):     ", HORIZON_YRS, "year(s)\n")
cat("Last eligible first-AP date:   ", format(FU_CUTOFF), "\n\n")

cohort_base <- build_ap_cohort(cmbspine_pancreatitis)
flow_base   <- attr(cohort_base, "flow")

dat_cox <- cohort_base %>%
  mutate(
    # -- predictors ----------------------------------------------------------
    age_cat = case_when(
      nv_k85_firstage < 50                          ~ 0,
      nv_k85_firstage >= 50 & nv_k85_firstage < 60 ~ 1,
      nv_k85_firstage >= 60 & nv_k85_firstage < 70 ~ 2,
      nv_k85_firstage >= 70                         ~ 3
    ),
    age_50_59  = as.integer(age_cat == 1),
    age_60_69  = as.integer(age_cat == 2),
    age_70plus = as.integer(age_cat == 3),
    male           = as.integer(sant_nv_sex == min(sant_nv_sex, na.rm = TRUE)),
    biliary_ap     = as.integer(nv_k851  == 1),
    alcohol_ap     = as.integer(nv_k852  == 1),
    unspecified_ap = as.integer(nv_k859  == 1),
    other_ap       = as.integer(
      coalesce(nv_k85own,   0L) == 1 |
        coalesce(nv_k85other, 0L) == 1 |
        coalesce(nv_k858,     0L) == 1 |
        coalesce(nv_k853,     0L) == 1
    )
  ) %>%
  select(event, surv_years, event_full, surv_years_full,
         age_cat, age_50_59, age_60_69, age_70plus, male,
         biliary_ap, alcohol_ap, unspecified_ap, other_ap) %>%
  drop_na()

# -- Flowchart (SANT) -------------------------------------------------------
flow_sant <- bind_rows(
  flow_base,
  tibble(Step = c("Excluded: missing predictor data (complete-case)",
                  "Analytic cohort",
                  "  Pancreatic cancer within 1 year (model outcome)",
                  "  Pancreatic cancer at any time after AP (reference)"),
         N = c(nrow(cohort_base) - nrow(dat_cox), nrow(dat_cox),
               sum(dat_cox$event), sum(dat_cox$event_full)))
)
cat("=== FLOWCHART (SANT) ===\n")
print(as.data.frame(flow_sant), row.names = FALSE)
write_release(flow_sant, list(N = character(0)), "flowchart_SANT_1yr.csv")
cat("Saved: flowchart_SANT_1yr.csv\n")
cat("NOTE: N and event counts now differ from 19,572 / 235 / 185 because of the\n")
cat("      Step 0C exclusions - update the manuscript numbers from this table.\n\n")

cat("=== SURVIVAL DATASET SUMMARY (1-YEAR FOLLOW-UP) ===\n")
cat("Total N:                       ", nrow(dat_cox), "\n")
cat("Cancer events, full follow-up: ", sum(dat_cox$event_full), "\n")
cat("Cancer events within 1 year:   ", sum(dat_cox$event), "\n")
cat(sprintf("  -> %.1f%% of all post-AP cancers occur within 1 year\n",
            100 * sum(dat_cox$event) / sum(dat_cox$event_full)))
cat("Cancers after 1 yr (censored): ", sum(dat_cox$event_full) - sum(dat_cox$event), "\n")
cat("No cancer within 1 year:       ", sum(dat_cox$event == 0), "\n")
cat("Died within 1 year, no cancer: ",
    sum(dat_cox$event == 0 & dat_cox$surv_years < HORIZON_YRS), "\n")
cat("1-year event rate:             ", round(mean(dat_cox$event) * 100, 2), "%\n\n")

cat("Same-day cancer diagnoses (surv_years = 1 day):",
    sum(dat_cox$event == 1 & dat_cox$surv_years <= 0.01), "\n")
cat("Later cancer diagnoses within 1 year:          ",
    sum(dat_cox$event == 1 & dat_cox$surv_years > 0.01), "\n\n")

# event_full is kept ONLY for the "% of all cancers" columns in Step 8A and
# Step 11B; the model, C-statistics, KM and DCA all use the 1-year `event`.
dat_cox <- dat_cox %>% select(-surv_years_full)


# =============================================================================
# STEP 2 — KAPLAN-MEIER OVERVIEW (0-12 MONTHS)
# =============================================================================

km_fit <- survfit(Surv(surv_years, event) ~ 1, data = dat_cox)

cat("=== KAPLAN-MEIER: CUMULATIVE CANCER RISK (1-YEAR) ===\n")
summary_times <- c(1/12, 3/12, 6/12, 9/12, 1)
km_summary    <- summary(km_fit, times = summary_times, extend = TRUE)
km_table <- data.frame(
  Months     = round(km_summary$time * 12),
  N_risk     = km_summary$n.risk,
  N_events   = km_summary$n.event,
  Cum_risk   = round((1 - km_summary$surv) * 100, 2),
  Lower_95CI = round((1 - km_summary$upper) * 100, 2),
  Upper_95CI = round((1 - km_summary$lower) * 100, 2)
)
print(km_table)

km_ymax <- max(0.005, (1 - min(km_fit$lower, na.rm = TRUE)) * 1.15)

png(file.path(OUT_DIR, "cox_km_overall.png"), width = 800, height = 550, res = 120)
plot(km_fit,
     fun      = "event",
     xlab     = "Years from first AP episode",
     ylab     = "Cumulative incidence of pancreatic cancer",
     main     = "Cumulative cancer incidence in the first year after AP\n(cmbspine_pancreatitis)",
     col      = "#2166AC", lwd = 2,
     conf.int = TRUE,
     xlim     = c(0, HORIZON_YRS),
     ylim     = c(0, km_ymax))
dev.off()


# =============================================================================
# STEP 3 — COX PROPORTIONAL HAZARDS REGRESSION (1-YEAR DATA)
# =============================================================================

cox_model <- coxph(
  Surv(surv_years, event) ~
    age_50_59 + age_60_69 + age_70plus + male +
    biliary_ap + alcohol_ap + unspecified_ap + other_ap,
  data = dat_cox,
  x    = TRUE
)

cat("\n=== COX MODEL SUMMARY (1-YEAR FOLLOW-UP) ===\n")
print(summary(cox_model))


# =============================================================================
# STEP 4 — TEST PROPORTIONAL HAZARDS ASSUMPTION
# =============================================================================

ph_test <- cox.zph(cox_model)
cat("\n=== PROPORTIONAL HAZARDS TEST ===\n")
print(ph_test)
cat("Variables with p < 0.05 may violate PH assumption.\n")
cat("Global p-value indicates whether model-wide PH holds.\n\n")

ph_table <- data.frame(
  Variable = rownames(ph_test$table),
  Chisq    = round(ph_test$table[, "chisq"], 2),
  df       = ph_test$table[, "df"],
  p        = format.pval(ph_test$table[, "p"], digits = 2, eps = 0.001),
  row.names = NULL
)
cat("=== SCHOENFELD RESIDUAL TEST (cox.zph) ===\n")
print(ph_table, row.names = FALSE)
write.csv(ph_table, file.path(OUT_DIR, "cox_schoenfeld_PH_test.csv"), row.names = FALSE)
cat("Saved: cox_schoenfeld_PH_test.csv\n\n")

png(file.path(OUT_DIR, "cox_schoenfeld.png"), width = 1000, height = 700, res = 120)
par(mfrow = c(3, 3))
plot(ph_test)
dev.off()


# =============================================================================
# STEP 5 — DERIVE INTEGER SCORE FROM COX COEFFICIENTS
# =============================================================================

cox_coefs <- coef(cox_model)

cat("=== RAW COX COEFFICIENTS ===\n")
print(round(data.frame(
  Log_HR = cox_coefs,
  HR     = exp(cox_coefs)
), 3))

ref_unit  <- min(abs(cox_coefs[abs(cox_coefs) > 0.10]))
ref_label <- names(cox_coefs)[abs(cox_coefs) == ref_unit][1]
cat(sprintf("\nReference unit: %.4f  [smallest |log-HR| > 0.10: %s]\n\n",
            ref_unit, ref_label))

# Previous (full-follow-up) score used a reference unit of 0.115.
if (abs(ref_unit - 0.115) / 0.115 > 0.20) {
  cat(sprintf("WARNING: reference unit changed from 0.115 to %.4f (>20%%).\n", ref_unit))
  cat("         Points are on a different scale - re-check the\n")
  cat("         Low/Moderate/High cut-offs before using them.\n\n")
} else {
  cat("Reference unit is within 20% of the previous 0.115 - point scale comparable.\n\n")
}

score_derivation <- data.frame(
  Variable = names(cox_coefs),
  Log_HR   = cox_coefs,
  HR       = exp(cox_coefs)
) %>%
  mutate(
    ref_unit   = ref_unit,
    Raw_points = Log_HR / ref_unit,
    Points     = as.integer(round(Raw_points))
  )

cat("=== SCORE DERIVATION TABLE ===\n")
print(score_derivation %>%
        select(Variable, Log_HR, HR, Raw_points, Points))
write.csv(score_derivation, file.path(OUT_DIR, "cox_score_derivation.csv"),
          row.names = FALSE)


# =============================================================================
# STEP 5A — TABLE S3: SCORE DERIVATION WITH HR (95% CI), P-VALUES AND
#           SCHOENFELD (PROPORTIONAL HAZARDS) TEST P-VALUES
# =============================================================================

s_cox <- summary(cox_model)

table_s3 <- data.frame(
  Variable     = rownames(s_cox$coefficients),
  Log_HR       = sprintf("%+.3f", s_cox$coefficients[, "coef"]),
  HR_CI        = sprintf("%.2f (%.2f\u2013%.2f)",
                         s_cox$conf.int[, 1], s_cox$conf.int[, 3], s_cox$conf.int[, 4]),
  p            = format.pval(s_cox$coefficients[, 5], digits = 2, eps = 0.001),
  row.names    = NULL
) %>%
  left_join(score_derivation %>%
              transmute(Variable,
                        Raw_points   = sprintf("%+.2f", Raw_points),
                        Final_points = sprintf("%+d", Points)),
            by = "Variable") %>%
  left_join(ph_table %>% transmute(Variable, PH_test_p = p),
            by = "Variable")

ph_global_p <- ph_table$p[ph_table$Variable == "GLOBAL"]

cat("=== TABLE S3: COX INTEGER RISK SCORE DERIVATION (1-YEAR MODEL) ===\n")
print(table_s3, row.names = FALSE)
cat(sprintf("\nReference unit (%s log-HR): %.4f\n", ref_label, ref_unit))
cat(sprintf("Global Schoenfeld test (proportional hazards): p = %s\n", ph_global_p))
cat("Baseline 1-year survival S0(1) is printed in Step 10A.\n\n")

var_labels <- c(
  age_50_59      = "Age 50\u201359 (ref < 50)",
  age_60_69      = "Age 60\u201369",
  age_70plus     = "Age \u226570",
  male           = "Male (ref female)",
  biliary_ap     = "Biliary AP",
  alcohol_ap     = "Alcohol AP",
  unspecified_ap = "Unspecified AP",
  other_ap       = "Other AP"
)

table_s3_out <- table_s3 %>%
  transmute(
    Variable        = unname(var_labels[Variable]),
    `Log-HR`        = Log_HR,
    `HR (95% CI)`   = HR_CI,
    `Raw points`    = Raw_points,
    `Final points`  = Final_points,
    `PH test p`     = PH_test_p,
    `Cox p`         = p
  )

cat("--- Table S3 (manuscript layout) ---\n")
print(table_s3_out, row.names = FALSE)
cat(sprintf("Footnote: reference unit = %s log-HR = %.4f; 1-year model; global PH test p = %s.\n\n",
            ref_label, ref_unit, ph_global_p))

write.csv(table_s3_out, file.path(OUT_DIR, "Table_S3_score_derivation.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table_S3_score_derivation.csv\n\n")


# =============================================================================
# STEP 6 — SCORING CARD DEFINITION AND PRINT
# =============================================================================

point_map <- setNames(score_derivation$Points, score_derivation$Variable)

# Achievable score range from THIS point map. Age bands are mutually
# exclusive; aetiology codes are not, so all negative codes can co-occur.
# The old -11 floor came from the full-follow-up model and is dropped: the
# score cannot fall below SCORE_MIN_ACHIEVABLE, so no floor is needed.
age_pts <- point_map[c("age_50_59", "age_60_69", "age_70plus")]
aet_pts <- point_map[c("biliary_ap", "alcohol_ap", "unspecified_ap", "other_ap")]
SCORE_MIN_ACHIEVABLE <- unname(min(0, age_pts) + min(0, point_map["male"]) + sum(pmin(0, aet_pts)))
SCORE_MAX_ACHIEVABLE <- unname(max(0, age_pts) + max(0, point_map["male"]) + sum(pmax(0, aet_pts)))
cat(sprintf("Achievable score range (1-year point map): %d to %d - no floor applied.\n\n",
            SCORE_MIN_ACHIEVABLE, SCORE_MAX_ACHIEVABLE))

print_score_card <- function(title) {
  cat("\n")
  cat("================================================================\n")
  cat(sprintf("   %s\n", title))
  cat(sprintf("   cmbspine_pancreatitis  (N = %d | %d cancers within 1 year)\n",
              nrow(dat_cox), sum(dat_cox$event)))
  cat("================================================================\n\n")
  cat("VARIABLE          CONDITION                          POINTS\n")
  cat("----------------------------------------------------------------\n")
  cat("Age at first AP   < 50 years                           0\n")
  cat(sprintf("                  50-59 years                         %+d\n", point_map["age_50_59"]))
  cat(sprintf("                  60-69 years                         %+d\n", point_map["age_60_69"]))
  cat(sprintf("                  >= 70 years                         %+d\n", point_map["age_70plus"]))
  cat("\n")
  cat("Sex               Female                               0\n")
  cat(sprintf("                  Male                                %+d\n", point_map["male"]))
  cat("\n")
  cat("AP aetiology      None of the codes below (reference)     0\n")
  cat(sprintf("                  Other                               %+d\n", point_map["other_ap"]))
  cat(sprintf("                  Alcohol-induced                     %+d\n", point_map["alcohol_ap"]))
  cat(sprintf("                  Unspecified                         %+d\n", point_map["unspecified_ap"]))
  cat(sprintf("                  Biliary                             %+d\n", point_map["biliary_ap"]))
  cat("----------------------------------------------------------------\n")
}

print_score_card("AP-CANCER COX RISK SCORE (1-YEAR MODEL) — SCORING CARD")
cat("\n")


# =============================================================================
# STEP 7 — APPLY SCORE TO DATASET (no floor; see Step 6)
# =============================================================================

dat_cox <- dat_cox %>%
  mutate(
    cox_score_raw =
      (age_50_59      * point_map["age_50_59"])      +
      (age_60_69      * point_map["age_60_69"])      +
      (age_70plus     * point_map["age_70plus"])     +
      (male           * point_map["male"])            +
      (biliary_ap     * point_map["biliary_ap"])      +
      (alcohol_ap     * point_map["alcohol_ap"])      +
      (unspecified_ap * point_map["unspecified_ap"])  +
      (other_ap       * point_map["other_ap"]),
    cox_score = cox_score_raw     # no floor: see Step 6
  )
stopifnot(all(dat_cox$cox_score >= SCORE_MIN_ACHIEVABLE),
          all(dat_cox$cox_score <= SCORE_MAX_ACHIEVABLE))

cat("Score distribution (counts 1-9 shown as <10):\n")
print(as.data.frame(sdc(dat_cox %>% count(cox_score, name = "n"),
                        list(n = character(0)), "score distribution")), row.names = FALSE)
cat("\nScore range:", min(dat_cox$cox_score), "to", max(dat_cox$cox_score), "\n\n")
if (max(dat_cox$cox_score) != SCORE_MAX)
  warning(sprintf("Max score is %d, not %d - update SCORE_MAX in Step 0B so the High label is correct.",
                  max(dat_cox$cox_score), SCORE_MAX))


# =============================================================================
# STEP 8 — 1-YEAR CANCER RATE BY SCORE
# =============================================================================

cancer_by_score <- dat_cox %>%
  group_by(cox_score) %>%
  summarise(
    n           = n(),
    n_cancer    = sum(event),
    cancer_rate = round(mean(event) * 100, 2),
    .groups     = "drop"
  ) %>%
  mutate(
    NNS = round(100 / pmax(cancer_rate, 0.001)),
    Recommendation = case_when(
      cancer_rate >= 2.0 ~ "Consider CT",
      cancer_rate >= 1.0 ~ "Enhanced follow-up",
      TRUE               ~ "Standard care"
    )
  )

cat("=== 1-YEAR CANCER RATE BY SCORE ===\n")
print(as.data.frame(cancer_by_score))
cat("\n")
write_release(cancer_by_score,
              list(n        = c("n_cancer", "cancer_rate", "NNS"),
                   n_cancer = c("cancer_rate", "NNS")),
              "cox_score_performance.csv")


# =============================================================================
# STEP 8A — POINT-BY-POINT 1-YEAR RISK TABLE (used to set the bands)
# =============================================================================

bh_tmp <- basehaz(cox_model, centered = FALSE)
S0_1y  <- exp(-max(bh_tmp$hazard[bh_tmp$time <= HORIZON_YRS], default = 0))
N_all  <- nrow(dat_cox); E_all <- sum(dat_cox$event); E_full <- sum(dat_cox$event_full)

risk_by_point <- dat_cox %>%
  group_by(score = cox_score) %>%
  group_modify(~ {
    s <- summary(survfit(Surv(surv_years, event) ~ 1, data = .x),
                 times = HORIZON_YRS, extend = TRUE)
    tibble(n                  = nrow(.x),
           cancers_1yr        = sum(.x$event),
           pct_of_1yr_cancers = 100 * sum(.x$event) / E_all,
           obs_pct            = 100 * mean(.x$event),
           cancers_all        = sum(.x$event_full),
           pct_of_all_cancers = 100 * sum(.x$event_full) / E_full,
           rate_all_pct       = 100 * mean(.x$event_full),
           km_pct  = 100 * (1 - s$surv),
           km_lo95 = 100 * (1 - s$upper),
           km_hi95 = 100 * (1 - s$lower))
  }) %>%
  ungroup() %>%
  arrange(score) %>%
  mutate(
    band            = assign_band(score),
    pred_pct        = 100 * (1 - S0_1y ^ exp(score * ref_unit)),
    cum_n_le        = cumsum(n),
    cum_rate_le_pct = 100 * cumsum(cancers_1yr) / cumsum(n),
    pct_flagged_ge  = 100 * rev(cumsum(rev(n))) / N_all,
    cancers_ge      = rev(cumsum(rev(cancers_1yr))),
    sens_ge_pct     = 100 * cancers_ge / E_all,
    rate_ge_pct     = 100 * cancers_ge / rev(cumsum(rev(n)))
  ) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

cat("=== STEP 8A: POINT-BY-POINT 1-YEAR RISK ===\n")
print(as.data.frame(risk_by_point))
# NOT SUPPRESSED (by request): full counts and cumulative columns, used to set
# the bands. Contains cells of 1-9 - keep it internal or clear it with the
# custodian before release.
write.csv(risk_by_point, file.path(OUT_DIR, "risk_by_point_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: risk_by_point_1yr.csv\n\n")


# =============================================================================
# STEP 9 — SENSITIVITY / SPECIFICITY BY CUTOFF (1-YEAR EVENTS)
# =============================================================================

total_events <- sum(dat_cox$event)
score_vals   <- sort(unique(dat_cox$cox_score))

thresh_df <- map_dfr(score_vals, function(thresh) {
  pred_pos <- as.integer(dat_cox$cox_score >= thresh)
  tp <- sum(pred_pos == 1 & dat_cox$event == 1)
  fp <- sum(pred_pos == 1 & dat_cox$event == 0)
  fn <- sum(pred_pos == 0 & dat_cox$event == 1)
  tn <- sum(pred_pos == 0 & dat_cox$event == 0)
  
  tibble(
    Cutoff         = thresh,
    Sensitivity    = round(tp / (tp + fn), 3),
    Specificity    = round(tn / (tn + fp), 3),
    PPV            = round(ifelse((tp+fp)==0, NA, tp/(tp+fp)), 3),
    NPV            = round(ifelse((tn+fn)==0, NA, tn/(tn+fn)), 3),
    TP = tp, FP = fp, FN = fn, TN = tn,
    N_flagged      = tp + fp,
    Pct_flagged    = round((tp + fp) / nrow(dat_cox) * 100, 1),
    Cancers_caught = tp,
    Cancers_missed = fn
  )
})

cat("\n=== SENSITIVITY / SPECIFICITY BY CUTOFF (1-YEAR) ===\n")
print(thresh_df %>%
        select(Cutoff, Sensitivity, Specificity,
               PPV, NPV, Cancers_caught, Cancers_missed, Pct_flagged) %>% as.data.frame())
# Cut-off rows are cumulative, so adjacent rows can be differenced to
# recover per-score counts - check the release copy by hand, or release
# Table S4 (band thresholds) instead of this table.
write_release(thresh_df,
              list(TP             = c("Sensitivity", "PPV", "NPV", "FN", "Cancers_missed"),
                   Cancers_caught = c("Sensitivity", "PPV", "NPV", "FN", "Cancers_missed"),
                   FN             = c("Sensitivity", "PPV", "NPV", "TP", "Cancers_caught"),
                   Cancers_missed = c("Sensitivity", "PPV", "NPV", "TP", "Cancers_caught"),
                   FP             = c("Specificity", "PPV", "TN", "N_flagged", "Pct_flagged"),
                   TN             = c("Specificity", "NPV", "FP"),
                   N_flagged      = c("Pct_flagged")),
              "cox_cutoff_performance.csv")


# =============================================================================
# STEP 9A — TABLE S4: DIAGNOSTIC PERFORMANCE OF MODERATE+HIGH AND HIGH-ONLY
# thresholds for cancer within 1 year, exact (Clopper-Pearson) 95% CIs.
# Patients censored before 1 year without cancer are classified negative.
# =============================================================================

diag_row <- function(label, flag) {
  y  <- dat_cox$event
  tp <- sum(flag & y == 1); fn <- sum(!flag & y == 1)
  fp <- sum(flag & y == 0); tn <- sum(!flag & y == 0)
  fmt <- function(x, n) {
    if (sdc_small(x) || sdc_small(n - x)) return("S")   # small-cell rule
    ci <- binom.test(x, n)$conf.int * 100
    sprintf("%.1f%% (%.1f\u2013%.1f)", 100 * x / n, ci[1], ci[2])
  }
  tibble(
    Threshold         = label,
    N_flagged         = if (sdc_small(sum(flag))) "<10" else
      if (sdc_small(fp)) "S" else format(sum(flag), big.mark = ","),
    Pct_flagged       = if (sdc_small(sum(flag)) || sdc_small(fp)) "S" else
      sprintf("%.1f%%", 100 * mean(flag)),
    Cancers_caught    = sprintf("%s / %d",
                                if (sdc_small(tp)) "<10" else if (sdc_small(fn)) "S" else tp,
                                tp + fn),
    `Sensitivity (95% CI)` = fmt(tp, tp + fn),
    `Specificity (95% CI)` = fmt(tn, tn + fp),
    `PPV (95% CI)`         = if (sdc_small(fn)) "S" else fmt(tp, tp + fp),
    `NPV (95% CI)`         = if (sdc_small(tp)) "S" else fmt(tn, tn + fn)
  )
}

table_s4 <- bind_rows(
  diag_row(sprintf("Moderate + High (\u2265%d)", BAND_LOW_MAX + 1),
           dat_cox$cox_score >= BAND_LOW_MAX + 1),
  diag_row(sprintf("High only (\u2265%d)", BAND_HIGH_MIN),
           dat_cox$cox_score >= BAND_HIGH_MIN)
)

cat("\n=== TABLE S4: DIAGNOSTIC PERFORMANCE (CANCER WITHIN 1 YEAR) ===\n")
print(as.data.frame(table_s4), row.names = FALSE)
write.csv(table_s4, file.path(OUT_DIR, "Table_S4_diagnostic_performance.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table_S4_diagnostic_performance.csv\n\n")


# =============================================================================
# STEP 10 — CONCORDANCE (C-STATISTIC) — 1-YEAR
# =============================================================================

c_stat <- concordance(cox_model)
cat(sprintf("\n1-year Harrell's C (Cox model, apparent): %.3f (SE = %.3f)\n",
            c_stat$concordance, sqrt(c_stat$var)))

roc_score <- roc(dat_cox$event, dat_cox$cox_score, quiet = TRUE)
cat(sprintf("Integer score AUROC (1-year events):      %.3f\n\n",
            as.numeric(auc(roc_score))))

png(file.path(OUT_DIR, "cox_roc.png"), width = 700, height = 600, res = 120)
plot(roc_score, col = "#2166AC", lwd = 2,
     main = sprintf("ROC - Cox Integer Score, 1-year events (AUC = %.3f)\ncmbspine_pancreatitis",
                    as.numeric(auc(roc_score))))
dev.off()


# =============================================================================
# STEP 10A — BASELINE SURVIVAL S0(t), t <= 1 YEAR (TRIPOD item 15a)
# P(cancer by t) = 1 - S0(t) ^ exp(LP)
# =============================================================================

cat("=== STEP 10A: BASELINE SURVIVAL S0(t) ===\n\n")

base_surv <- basehaz(cox_model, centered = FALSE)

time_points <- c(3/12, 6/12, 9/12, 1)
s0_table <- map_dfr(time_points, function(t) {
  h0 <- max(base_surv$hazard[base_surv$time <= t], default = 0)
  tibble(
    Months                 = round(t * 12),
    Years                  = t,
    Cum_Baseline_Hazard_H0 = round(h0, 6),
    Baseline_Survival_S0   = round(exp(-h0), 6)
  )
})

cat("Baseline survival at key time points (centered = FALSE):\n")
print(s0_table)
cat("\nInterpretation:\n")
cat("  LP = score * reference_unit\n")
cat("  P(cancer by t) = 1 - S0(t) ^ exp(LP)\n\n")

s0_1yr <- s0_table$Baseline_Survival_S0[s0_table$Years == 1]
lp_ex  <- 2 * ref_unit
p_ex   <- 1 - s0_1yr ^ exp(lp_ex)
cat(sprintf("  Example: score = 2  ->  LP = %.4f  ->  1-yr cancer risk = %.2f%%\n\n",
            lp_ex, p_ex * 100))

# Confirm S0 is UNCENTRED (= survival for a score-0 patient: female, <50,
# no aetiology code). survfit() for that patient must equal S0(1yr) exactly;
# R's default centred baseline would not.
ref_pt <- dat_cox[1, ] %>%
  mutate(across(c(age_50_59, age_60_69, age_70plus, male,
                  biliary_ap, alcohol_ap, unspecified_ap, other_ap), ~ 0L))
s_ref <- summary(survfit(cox_model, newdata = ref_pt),
                 times = HORIZON_YRS, extend = TRUE)$surv
stopifnot(abs(s_ref - S0_1y) < 1e-10)
cat(sprintf("Confirmed: S0(1yr) = %.6f is UNCENTRED (matches survfit for score 0).\n\n", S0_1y))

write.csv(s0_table, file.path(OUT_DIR, "cox_baseline_survival.csv"), row.names = FALSE)
saveRDS(s0_table,   file.path(OUT_DIR, "cox_baseline_survival.rds"))
cat("Baseline survival saved to cox_baseline_survival.csv / .rds\n\n")


# =============================================================================
# STEP 10B — BOOTSTRAP INTERNAL VALIDATION, COX MODEL (TRIPOD item 10b)
# Because the data are truncated at 1 year, this is the optimism-corrected
# 1-YEAR Harrell's C for the Cox model.
# =============================================================================

cat("=== STEP 10B: BOOTSTRAP INTERNAL VALIDATION (B = 500) ===\n\n")

set.seed(123)

dd <- datadist(dat_cox)
options(datadist = "dd")

cox_rms <- cph(
  Surv(surv_years, event) ~
    age_50_59 + age_60_69 + age_70plus + male +
    biliary_ap + alcohol_ap + unspecified_ap + other_ap,
  data = dat_cox,
  x    = TRUE,
  y    = TRUE,
  surv = TRUE
)

cat("rms::cph model fitted. Running bootstrap validation (B = 500)...\n")
boot_val <- validate(cox_rms, method = "boot", B = 500)

cat("\n=== BOOTSTRAP VALIDATION RESULTS ===\n")
print(boot_val)

# Dxy (Somers' D) to C-statistic: C = (Dxy / 2) + 0.5
apparent_Dxy   <- boot_val["Dxy", "index.orig"]
corrected_Dxy  <- boot_val["Dxy", "index.corrected"]
optimism_Dxy   <- boot_val["Dxy", "optimism"]

apparent_C     <- (apparent_Dxy  / 2) + 0.5
corrected_C    <- (corrected_Dxy / 2) + 0.5
optimism_C     <- optimism_Dxy   / 2

shrinkage      <- boot_val["Slope", "index.corrected"]

cat("\n=== OPTIMISM-CORRECTED PERFORMANCE (1-YEAR, COX MODEL) ===\n")
cat(sprintf("  Apparent 1-yr Harrell's C:           %.3f\n", apparent_C))
cat(sprintf("  Optimism:                            %.3f\n", optimism_C))
cat(sprintf("  Optimism-corrected 1-yr Harrell's C: %.3f\n", corrected_C))
cat(sprintf("  Calibration slope (shrinkage):       %.3f\n", shrinkage))
if (shrinkage < 0.95) {
  cat(sprintf("  NOTE: Shrinkage = %.3f < 0.95 - moderate overfitting.\n", shrinkage))
  cat("        Consider multiplying all log-HRs by the shrinkage factor\n")
  cat("        before deriving integer score points.\n")
} else {
  cat("  NOTE: Shrinkage >= 0.95 - minimal overfitting, coefficients stable.\n")
}
cat("\n")

boot_summary <- data.frame(
  Metric = c("Apparent 1-yr Harrell's C",
             "Optimism",
             "Optimism-corrected 1-yr Harrell's C",
             "Calibration slope (shrinkage)"),
  Value  = round(c(apparent_C, optimism_C, corrected_C, shrinkage), 3)
)
write.csv(boot_summary, file.path(OUT_DIR, "cox_bootstrap_validation.csv"),
          row.names = FALSE)
cat("Saved: cox_bootstrap_validation.csv\n\n")


# =============================================================================
# STEP 10C — OPTIMISM-CORRECTED 1-YEAR HARRELL'S C FOR THE INTEGER SCORE
# rms::validate (10B) covers the continuous Cox model only. This repeats the
# whole score-building process in each bootstrap sample (refit -> reference
# unit -> round), so the optimism includes the rounding.
# =============================================================================

cat("=== STEP 10C: 1-YEAR HARRELL'S C, INTEGER SCORE (B = 500) ===\n\n")

B_C1Y <- 500

cox_form <- Surv(surv_years, event) ~
  age_50_59 + age_60_69 + age_70plus + male +
  biliary_ap + alcohol_ap + unspecified_ap + other_ap

pred_vars <- c("age_50_59", "age_60_69", "age_70plus", "male",
               "biliary_ap", "alcohol_ap", "unspecified_ap", "other_ap")

harrell_c <- function(d, risk) {
  cc <- concordance(Surv(surv_years, event) ~ risk,
                    data = data.frame(surv_years = d$surv_years,
                                      event = d$event, risk = risk),
                    reverse = TRUE)
  c(C = unname(cc$concordance), SE = sqrt(unname(cc$var)))
}

build_scores <- function(fit, newdat) {
  b   <- coef(fit)[pred_vars]
  X   <- as.matrix(newdat[, pred_vars])
  big <- abs(b[abs(b) > 0.10])
  ru  <- if (length(big) > 0) min(big) else ref_unit
  pts <- round(b / ru)
  list(lp    = as.vector(X %*% b),
       score = as.vector(X %*% pts))
}

app_sc <- build_scores(cox_model, dat_cox)
stopifnot(isTRUE(all.equal(app_sc$score, as.numeric(dat_cox$cox_score),
                           check.attributes = FALSE)))

app_lp_c    <- harrell_c(dat_cox, app_sc$lp)
app_score_c <- harrell_c(dat_cox, app_sc$score)

set.seed(2027)
n_c1y   <- nrow(dat_cox)
opt_c1y <- matrix(NA_real_, B_C1Y, 2, dimnames = list(NULL, c("lp", "score")))

cat(sprintf("Running bootstrap (B = %d)...\n", B_C1Y))
for (i in seq_len(B_C1Y)) {
  bd    <- dat_cox[sample.int(n_c1y, n_c1y, replace = TRUE), ]
  fit_b <- tryCatch(coxph(cox_form, data = bd), error = function(e) NULL)
  if (is.null(fit_b) || anyNA(coef(fit_b))) next
  
  s_boot <- build_scores(fit_b, bd)
  s_orig <- build_scores(fit_b, dat_cox)
  
  opt_c1y[i, "lp"]    <- harrell_c(bd, s_boot$lp)["C"]    - harrell_c(dat_cox, s_orig$lp)["C"]
  opt_c1y[i, "score"] <- harrell_c(bd, s_boot$score)["C"] - harrell_c(dat_cox, s_orig$score)["C"]
  
  if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, B_C1Y))
}

optim_lp    <- mean(opt_c1y[, "lp"],    na.rm = TRUE)
optim_score <- mean(opt_c1y[, "score"], na.rm = TRUE)
corr_C1y_lp    <- app_lp_c["C"]    - optim_lp
corr_C1y_score <- app_score_c["C"] - optim_score

c1y_table <- data.frame(
  Predictor         = c("Cox linear predictor", "PAnCReAS integer score"),
  Events_1yr        = sum(dat_cox$event),
  Apparent_C        = round(c(app_lp_c["C"], app_score_c["C"]), 3),
  Apparent_SE       = round(c(app_lp_c["SE"], app_score_c["SE"]), 3),
  Optimism          = round(c(optim_lp, optim_score), 3),
  Corrected_C       = round(c(corr_C1y_lp, corr_C1y_score), 3),
  Corrected_lo95    = round(c(corr_C1y_lp    - 1.96 * app_lp_c["SE"],
                              corr_C1y_score - 1.96 * app_score_c["SE"]), 3),
  Corrected_hi95    = round(c(corr_C1y_lp    + 1.96 * app_lp_c["SE"],
                              corr_C1y_score + 1.96 * app_score_c["SE"]), 3),
  Bootstrap_reps_ok = sum(complete.cases(opt_c1y)),
  row.names = NULL
)

cat("\n=== 1-YEAR HARRELL'S C: APPARENT vs OPTIMISM-CORRECTED ===\n")
print(c1y_table, row.names = FALSE)
cat("\nThe Cox-LP row should closely match Step 10B (different bootstrap draws).\n")
cat("Tied integer scores count as 0.5, so the integer score's C sits slightly lower.\n")
cat("95% CI = apparent CI shifted by the optimism (approximate).\n\n")

write.csv(c1y_table, file.path(OUT_DIR, "cox_1yr_harrell_C_optimism_corrected.csv"),
          row.names = FALSE)
saveRDS(list(table = c1y_table, optimism_draws = opt_c1y),
        file.path(OUT_DIR, "cox_1yr_harrell_C_optimism_corrected.rds"))
cat("Saved: cox_1yr_harrell_C_optimism_corrected.csv / .rds\n\n")


# =============================================================================
# STEP 11 — KM CURVES BY SCORE GROUP (0-12 MONTHS)
# =============================================================================

dat_cox <- dat_cox %>%
  mutate(
    risk_band   = assign_band(cox_score),
    score_group = factor(unname(BAND_LABELS[as.character(risk_band)]),
                         levels = unname(BAND_LABELS))
  )

km_group <- survfit(Surv(surv_years, event) ~ score_group, data = dat_cox)

km_grp_ymax <- max(0.005, (1 - min(km_group$surv, na.rm = TRUE)) * 1.2)

png(file.path(OUT_DIR, "cox_km_by_score.png"), width = 900, height = 620, res = 120)
plot(km_group,
     fun      = "event",
     xlab     = "Years from first AP episode",
     ylab     = "Cumulative incidence of pancreatic cancer",
     main     = "1-year cumulative cancer incidence by risk score group\n(cmbspine_pancreatitis)",
     col      = unname(BAND_COLOURS),
     lwd      = 2,
     conf.int = FALSE,
     xlim     = c(0, HORIZON_YRS),
     ylim     = c(0, km_grp_ymax))
legend("topleft",
       legend = levels(dat_cox$score_group),
       col    = unname(BAND_COLOURS),
       lwd    = 2, bty = "n", cex = 0.9)
dev.off()

cat("1-year cumulative cancer risk by score group:\n")
print(summary(km_group, times = HORIZON_YRS, extend = TRUE)$table)


# =============================================================================
# STEP 11B — 1-YEAR RESTRICTED CUMULATIVE INCIDENCE BY RISK BAND
# =============================================================================

cat("\n=== STEP 11B: 1-YEAR RESTRICTED CUMULATIVE INCIDENCE BY BAND ===\n\n")

band_cancer_dist <- dat_cox %>%
  group_by(risk_band) %>%
  summarise(
    n                  = n(),
    pct_of_cohort      = 100 * n() / nrow(dat_cox),
    cancers_1yr        = sum(event),
    pct_of_1yr_cancers = 100 * sum(event) / sum(dat_cox$event),
    rate_1yr_pct       = 100 * mean(event),
    cancers_all        = sum(event_full),
    pct_of_all_cancers = 100 * sum(event_full) / sum(dat_cox$event_full),
    rate_all_pct       = 100 * mean(event_full),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.double), ~ round(.x, 2)))

cat("=== Risk band: share of 1-year cancers and of ALL cancers ===\n")
cat(sprintf("(1-year cancers = %d; all post-AP cancers, full follow-up = %d)\n",
            sum(dat_cox$event), sum(dat_cox$event_full)))
print(as.data.frame(band_cancer_dist))
write_release(band_cancer_dist,
              list(n           = c("pct_of_cohort"),
                   cancers_1yr = c("pct_of_1yr_cancers", "rate_1yr_pct"),
                   cancers_all = c("pct_of_all_cancers", "rate_all_pct")),
              "band_cancer_distribution.csv")
cat("Saved: band_cancer_distribution.csv\n\n")

horizon_yrs <- HORIZON_YRS

restricted_ci_by_band <- dat_cox %>%
  group_by(risk_band) %>%
  group_modify(~ {
    band_dat <- .x
    n_band   <- nrow(band_dat)
    
    naive_1yr_cases <- sum(band_dat$event == 1 & band_dat$surv_years <= horizon_yrs)
    naive_1yr_rate  <- naive_1yr_cases / n_band * 100
    
    fit <- tryCatch(
      survfit(Surv(surv_years, event) ~ 1, data = band_dat),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      km_rate <- NA_real_; km_lo <- NA_real_; km_hi <- NA_real_; n_risk_1yr <- NA_integer_
    } else {
      s <- summary(fit, times = horizon_yrs, extend = TRUE)
      km_rate    <- (1 - s$surv)  * 100
      km_lo      <- (1 - s$upper) * 100
      km_hi      <- (1 - s$lower) * 100
      n_risk_1yr <- s$n.risk
    }
    
    tibble(
      n_band           = n_band,
      n_cases_1yr      = naive_1yr_cases,
      naive_1yr_rate   = round(naive_1yr_rate, 2),
      km_1yr_rate      = round(km_rate, 2),
      km_1yr_lo95      = round(km_lo, 2),
      km_1yr_hi95      = round(km_hi, 2),
      n_at_risk_1yr    = n_risk_1yr,
      naive_vs_km_diff = round(naive_1yr_rate - km_rate, 3)
    )
  }) %>%
  ungroup()

cat("1-year restricted cumulative incidence by band (naive vs KM):\n")
print(restricted_ci_by_band)
cat("\nnaive_vs_km_diff near 0 = minimal early censoring (deaths) in that band.\n\n")

write_release(restricted_ci_by_band,
              list(n_cases_1yr = c("naive_1yr_rate", "km_1yr_rate", "km_1yr_lo95",
                                   "km_1yr_hi95", "naive_vs_km_diff"),
                   n_at_risk_1yr = character(0)),
              "cox_1yr_restricted_incidence_by_PAnCReAS_band.csv")
cat("Saved: cox_1yr_restricted_incidence_by_PAnCReAS_band.csv\n\n")

# -----------------------------------------------------------------------------
# STEP 11B-ii — TABLE 3: OBSERVED 1-YEAR CANCER RISK BY RISK BAND
# Score | Band | N | Cancers (1-yr) | AP cohort assigned (%) |
# 1-year risk, KM (95% CI) | Number needed to image (NNI = 100 / KM %)
# -----------------------------------------------------------------------------
score_range_lab <- c(
  Low      = sprintf("≤%d", BAND_LOW_MAX),
  Moderate = sprintf("%d–%d", BAND_LOW_MAX + 1, BAND_HIGH_MIN - 1),
  High     = sprintf("%d–%d", BAND_HIGH_MIN, SCORE_MAX)
)

table_3 <- dat_cox %>%
  group_by(risk_band) %>%
  summarise(N = n(), Cancers = sum(event), .groups = "drop") %>%
  left_join(restricted_ci_by_band %>%
              select(risk_band, km_1yr_rate, km_1yr_lo95, km_1yr_hi95),
            by = "risk_band") %>%
  transmute(
    Score                       = unname(score_range_lab[as.character(risk_band)]),
    Band                        = as.character(risk_band),
    N                           = format(N, big.mark = ","),
    Cancers,
    `AP cohort assigned (%)`    = sprintf("%.1f", 100 * as.numeric(gsub(",", "", N)) / nrow(dat_cox)),
    `1-year risk, KM (95% CI)`  = sprintf("%.2f%% (%.2f–%.2f%%)",
                                          km_1yr_rate, km_1yr_lo95, km_1yr_hi95),
    `Number needed to image (NNI)` = ifelse(km_1yr_rate > 0,
                                            as.character(round(100 / km_1yr_rate)), "—")
  )

ph_viol <- rownames(ph_test$table)[ph_test$table[, "p"] < 0.05 &
                                     rownames(ph_test$table) != "GLOBAL"]
table_3_footnote <- sprintf(paste0(
  "NNI, number needed to image; KM, Kaplan-Meier. Cox proportional hazards model ",
  "fitted on 1-year follow-up (SANT); apparent C-statistic = %.3f; optimism-corrected ",
  "C-statistic = %.3f (bootstrap B = 500); calibration slope (shrinkage) = %.3f. ",
  "Reference unit = %s log-HR (%.4f). PH violations (Schoenfeld p < 0.05): %s; global p = %s. ",
  "The 1-year rate is the KM estimate, accounting for censoring and deaths. ",
  "Bands: Low %s, Moderate %s, High %s. AP cohort assigned (%%) = percentage of the ",
  "total AP cohort (N = %s) in each band; %d cancers within 1 year."),
  apparent_C, corrected_C, shrinkage, ref_label, ref_unit,
  if (length(ph_viol)) paste(ph_viol, collapse = ", ") else "none", ph_global_p,
  score_range_lab["Low"], score_range_lab["Moderate"], score_range_lab["High"],
  format(nrow(dat_cox), big.mark = ","), sum(dat_cox$event))

cat("=== TABLE 3: COX INTEGER RISK SCORE - OBSERVED 1-YEAR CANCER RISK BY BAND ===\n")
print(as.data.frame(table_3), row.names = FALSE)
cat("\nFootnote:", table_3_footnote, "\n\n")

write_release(table_3,
              list(N       = c("AP cohort assigned (%)"),
                   Cancers = c("1-year risk, KM (95% CI)", "Number needed to image (NNI)")),
              "Table3_risk_bands_1yr.csv")
writeLines(table_3_footnote, file.path(OUT_DIR, "Table3_footnote.txt"))
cat("Saved: Table3_risk_bands_1yr.csv / Table3_footnote.txt\n\n")

horizons       <- c(3/12, 6/12, 9/12, 1)
horizon_labels <- c("3_months", "6_months", "9_months", "1_year")

multi_horizon <- map2_dfr(horizons, horizon_labels, function(h, lbl) {
  dat_cox %>%
    group_by(risk_band) %>%
    group_modify(~ {
      fit <- survfit(Surv(surv_years, event) ~ 1, data = .x)
      s <- summary(fit, times = h, extend = TRUE)
      tibble(
        horizon     = lbl,
        n_band      = nrow(.x),
        n_events    = sum(.x$event == 1 & .x$surv_years <= h),
        km_rate_raw = (1 - s$surv) * 100,            # unrounded, for figure labels
        km_rate_pct = round((1 - s$surv) * 100, 2),
        km_lo95     = round((1 - s$upper) * 100, 2),
        km_hi95     = round((1 - s$lower) * 100, 2)
      )
    }) %>%
    ungroup()
})

cat("=== Multi-horizon KM cumulative incidence by band (3/6/9/12 months) ===\n")
print(multi_horizon %>% arrange(risk_band, horizon))
write_release(multi_horizon,
              list(n_events = c("km_rate_pct", "km_lo95", "km_hi95")),
              "cox_multihorizon_incidence_by_PAnCReAS_band.csv",
              drop_cols = "km_rate_raw")
cat("Saved: cox_multihorizon_incidence_by_PAnCReAS_band.csv\n\n")


# =============================================================================
# STEP 11C — FIGURE 1: KM CUMULATIVE INCIDENCE BY BAND (3, 6, 9, 12 MONTHS)
# =============================================================================

# Round half up (R's sprintf/round use round-half-to-even, so 4.25 -> 4.2).
# Labels are built from the UNROUNDED KM estimate, not the 2-dp value, to
# avoid double rounding.
round_half_up <- function(x, digits = 1) {
  m <- 10^digits
  sign(x) * floor(abs(x) * m + 0.5 + 1e-9) / m
}

cat("Figure 1 labels (unrounded KM % -> label):\n")
print(multi_horizon %>%
        transmute(risk_band, horizon, km_rate_raw = round(km_rate_raw, 4),
                  label = sprintf("%.1f%%", round_half_up(km_rate_raw, 1))) %>%
        as.data.frame())

fig1_dat <- multi_horizon %>%
  mutate(
    horizon = factor(horizon,
                     levels = c("3_months", "6_months", "9_months", "1_year"),
                     labels = c("3 months", "6 months", "9 months", "1 year")),
    band_lab = factor(unname(BAND_LABELS[as.character(risk_band)]),
                      levels = unname(BAND_LABELS))
  )

if (any(sdc_small(fig1_dat$n_events)))
  warning("Figure 1 has band x horizon points with 1-9 cancers: labels are suppressed, ",
          "but the points can still be read off the axis - check with the custodian.")

p_fig1 <- ggplot(fig1_dat,
                 aes(x = horizon, y = km_rate_pct,
                     colour = band_lab, group = band_lab)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  geom_text(aes(label = ifelse(sdc_small(n_events), "<10 cases",
                               sprintf("%.1f%%", round_half_up(km_rate_raw, 1)))),
            vjust = -1, size = 3.2, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = setNames(unname(BAND_COLOURS), unname(BAND_LABELS)),
                      name = "PAnCReAS band") +
  scale_y_continuous(limits = c(0, max(fig1_dat$km_rate_pct) * 1.15),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Time since AP", y = "KM cumulative incidence (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(OUT_DIR, "Figure1_KM_incidence_by_PAnCReAS_band.png"),
       p_fig1, width = 8, height = 5.5, dpi = 400)
cat("Saved: Figure1_KM_incidence_by_PAnCReAS_band.png\n\n")


# =============================================================================
# STEP 12 — 1-YEAR CANCER RATE BY SCORE PLOT
# =============================================================================

# Scores with 1-9 cancers (or 1-9 patients): bar and rate hidden
cancer_by_score_plot <- cancer_by_score %>%
  mutate(supp      = sdc_small(n_cancer) | sdc_small(n),
         rate_plot = ifelse(supp, NA_real_, cancer_rate),
         lab_y     = ifelse(supp, 0, cancer_rate),
         lab       = ifelse(supp,
                            paste0("<10 cancers\nn=", sdc_lab(n)),
                            paste0(cancer_rate, "%\n(NNS=", NNS, ")\nn=",
                                   format(n, big.mark = ","))))

p_score <- ggplot(
  cancer_by_score_plot,
  aes(x    = factor(cox_score),
      y    = rate_plot,
      fill = Recommendation)
) +
  geom_col(alpha = 0.85, colour = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_text(
    aes(y = lab_y, label = lab),
    vjust = -0.25, size = 2.9, lineheight = 0.9, show.legend = FALSE
  ) +
  geom_hline(yintercept = 1.0, linetype = "dotted",
             colour = "#2166AC", linewidth = 0.8) +
  geom_hline(yintercept = 2.0, linetype = "dashed",
             colour = "#FF7F00", linewidth = 0.8) +
  scale_fill_manual(values = c(
    "Consider CT"        = "#FEF9E7",
    "Enhanced follow-up" = "#EBF5FB",
    "Standard care"      = "#F2F3F4"
  )) +
  scale_y_continuous(
    limits = c(0, max(cancer_by_score$cancer_rate, na.rm = TRUE) * 1.45)
  ) +
  labs(
    title    = "AP-Cancer Cox Risk Score (1-year model) - Observed 1-Year Cancer Rate by Score",
    subtitle = sprintf(
      "cmbspine_pancreatitis  |  N = %d  |  %d cancers within 1 year (%.2f%%)",
      nrow(dat_cox), sum(dat_cox$event), mean(dat_cox$event) * 100
    ),
    x       = "Cox Integer Risk Score",
    y       = "Observed 1-Year Cancer Rate (%)",
    fill    = "Recommendation",
    caption = paste0(
      "Dotted = 1% threshold (enhanced follow-up).  ",
      "Dashed = 2% threshold (consider CT).\n",
      "NNS = number of AP patients needed to scan to detect one cancer.  ",
      "Reference unit = smallest log-HR > 0.10.  ",
      "Bars with 1-9 cancers suppressed (<10 rule)."
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.caption    = element_text(colour = "grey40", size = 9),
    plot.title      = element_text(face = "bold"),
    axis.text.x     = element_text(size = 9)
  )

ggsave(file.path(OUT_DIR, "cox_score_cancer_rate.png"), p_score,
       width = 10, height = 7, dpi = 300)


# =============================================================================
# STEP 13 — SENSITIVITY / SPECIFICITY PLOT
# =============================================================================

p_thresh <- thresh_df %>%
  pivot_longer(cols      = c(Sensitivity, Specificity),
               names_to  = "Metric",
               values_to = "Value") %>%
  mutate(Value = case_when(
    Metric == "Sensitivity" & (sdc_small(TP) | sdc_small(FN)) ~ NA_real_,
    Metric == "Specificity" & (sdc_small(FP) | sdc_small(TN)) ~ NA_real_,
    TRUE ~ Value)) %>%
  ggplot(aes(x      = factor(Cutoff),
             y      = Value,
             colour = Metric,
             group  = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_text(
    data = thresh_df,
    aes(x     = factor(Cutoff),
        y     = -0.12,
        label = paste0(ifelse(sdc_small(Cancers_caught), "<10",
                              ifelse(sdc_small(Cancers_missed), "S", Cancers_caught)), "/",
                       total_events, "\n",
                       ifelse(sdc_small(N_flagged) | sdc_small(FP), "S",
                              paste0(Pct_flagged, "%")))),
    inherit.aes = FALSE,
    size = 2.8, colour = "grey30", lineheight = 0.85, show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Sensitivity" = "#E41A1C", "Specificity" = "#2166AC"),
    name   = "Metric"
  ) +
  scale_y_continuous(limits = c(-0.18, 1.05),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(
    title   = "AP-Cancer Cox Risk Score (1-year model) - Sensitivity and Specificity by Cutoff",
    x       = "Score Cutoff (flag for CT if score >= cutoff)",
    y       = "Metric Value",
    caption = paste0(
      "Numbers below axis: 1-year cancers caught / total  |  % of patients flagged.\n",
      "Select cutoff balancing sensitivity requirement vs CT resource constraints."
    )
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.caption    = element_text(colour = "grey40", size = 9),
    plot.title      = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "cox_score_sensitivity_specificity.png"), p_thresh,
       width = 10, height = 7, dpi = 300)


# =============================================================================
# STEP 14 — DECISION CURVE ANALYSIS (PREDICTED 1-YEAR RISK, TIME-TO-EVENT)
# Net benefit uses KM risk at 1 year, so deaths/censoring before 1 year are
# handled properly rather than counted as cancer-free.
# P(cancer by 1yr) = 1 - exp(-H0(1) * exp(LP))
# =============================================================================

bh    <- basehaz(cox_model, centered = FALSE)
H0_1y <- max(bh$hazard[bh$time <= HORIZON_YRS], default = 0)

b <- coef(cox_model)
lp_raw <- with(dat_cox,
               age_50_59      * b["age_50_59"] +
                 age_60_69      * b["age_60_69"] +
                 age_70plus     * b["age_70plus"] +
                 male           * b["male"] +
                 biliary_ap     * b["biliary_ap"] +
                 alcohol_ap     * b["alcohol_ap"] +
                 unspecified_ap * b["unspecified_ap"] +
                 other_ap       * b["other_ap"])

pred_1yr <- 1 - exp(-H0_1y * exp(lp_raw))

chk <- 1 - summary(survfit(cox_model, newdata = dat_cox[1:5, ]),
                   times = HORIZON_YRS, extend = TRUE)$surv
cat("Check: direct 1-yr risk vs survfit() (first 5 patients) - columns should match\n")
print(round(cbind(survfit = as.vector(chk), direct = as.vector(pred_1yr[1:5])), 5))
stopifnot(max(abs(as.vector(chk) - as.vector(pred_1yr[1:5]))) < 1e-8)
cat("\n")

dat_cox <- dat_cox %>%
  mutate(
    cox_prob_1yr = as.vector(pred_1yr),
    # Integer score mapped through the model's own S0(1yr) (Step 10A)
    cox_score_calibrated = 1 - s0_1yr ^ exp(cox_score * ref_unit),
    # Clinical rule: non-gallstone/alcohol aetiology OR age >= 50y
    clinical_rule = as.integer(
      (biliary_ap == 0 & alcohol_ap == 0) | age_cat >= 1
    )
  )

cat("Clinical rule positive (flagged for scan):",
    sum(dat_cox$clinical_rule == 1), "/", nrow(dat_cox),
    sprintf("(%.1f%%)\n\n", mean(dat_cox$clinical_rule) * 100))

cat("Calibrated integer score -> predicted 1-year risk, by score value:\n")
print(dat_cox %>%
        distinct(cox_score, cox_score_calibrated) %>%
        arrange(cox_score) %>%
        mutate(cox_score_calibrated = round(cox_score_calibrated * 100, 3)) %>%
        rename(predicted_1yr_risk_pct = cox_score_calibrated) %>% as.data.frame())
cat("\n")

dca_result <- dca(
  Surv(surv_years, event) ~ cox_prob_1yr + cox_score_calibrated + clinical_rule,
  data       = dat_cox %>%
    select(surv_years, event, cox_prob_1yr, cox_score_calibrated, clinical_rule) %>%
    drop_na(),
  time       = HORIZON_YRS,
  thresholds = seq(0.005, 0.10, by = 0.005),
  label      = list(
    cox_prob_1yr          = "PAnCReAS (Cox 1-yr probability)",
    cox_score_calibrated  = "PAnCReAS (integer score, calibrated)",
    clinical_rule         = "Clinical rule"
  )
)

cat("\n=== DECISION CURVE ANALYSIS SUMMARY (1-YEAR) ===\n")
print(summary(dca_result))

p_dca <- plot(dca_result, smooth = TRUE) +
  labs(
    title    = "Decision Curve Analysis - AP-Cancer Cox Risk Score (1-year)",
    subtitle = sprintf(
      "cmbspine_pancreatitis  |  N = %d  |  %d cancers within 1 year (%.2f%%)",
      nrow(dat_cox), sum(dat_cox$event), mean(dat_cox$event) * 100
    ),
    x = "Threshold Probability (1-year risk)",
    y = "Net Benefit"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold")
  )

ggsave(file.path(OUT_DIR, "cox_dca.png"), p_dca, width = 10, height = 7, dpi = 300)


# =============================================================================
# STEP 14A — DCA SUMMARY TABLE (NB AT 1% / 2% / 3% / 4% / 5%)
# =============================================================================

cat("\n=== STEP 14A: DCA SUMMARY TABLE (NB at 1% / 2% / 3% / 4% / 5%) ===\n\n")

dca_long <- as_tibble(dca_result$dca)

strategy_labels <- c(
  "cox_score_calibrated" = "PAnCReAS (integer score, calibrated)",
  "cox_prob_1yr"         = "PAnCReAS (Cox 1-yr probability)",
  "clinical_rule"        = "Clinical rule",
  "all"                  = "Scan-All",
  "none"                 = "Scan-None"
)

dca_summary_table <- dca_long %>%
  filter(variable %in% names(strategy_labels)) %>%
  filter(sapply(threshold, function(t) any(abs(t - c(0.01, 0.02, 0.03, 0.04, 0.05)) < 1e-8))) %>%
  mutate(
    Strategy   = unname(strategy_labels[variable]),
    thresh_lbl = paste0("NB_at_", round(threshold * 100), "pct")
  ) %>%
  select(Strategy, thresh_lbl, net_benefit) %>%
  pivot_wider(names_from = thresh_lbl, values_from = net_benefit) %>%
  mutate(Strategy = factor(Strategy, levels = unname(strategy_labels))) %>%
  arrange(Strategy)

cat("Decision curve analysis - net benefit by strategy:\n")
print(dca_summary_table)

write.csv(dca_summary_table, file.path(OUT_DIR, "cox_dca_summary_table.csv"),
          row.names = FALSE)
cat("\nSaved: cox_dca_summary_table.csv\n\n")

# -----------------------------------------------------------------------------
# STEP 14B — TABLE S5: DCA (manuscript layout)
# PAnCReAS = integer score converted to predicted 1-year risk via S0(1yr).
# -----------------------------------------------------------------------------
table_s5 <- dca_summary_table %>%
  filter(Strategy %in% c("PAnCReAS (Cox 1-yr probability)",
                         "PAnCReAS (integer score, calibrated)", "Clinical rule",
                         "Scan-All", "Scan-None")) %>%
  mutate(Strategy = factor(as.character(Strategy),
                           levels = c("PAnCReAS (Cox 1-yr probability)",
                                      "PAnCReAS (integer score, calibrated)",
                                      "Clinical rule", "Scan-All", "Scan-None"))) %>%
  arrange(Strategy) %>%
  mutate(Strategy = recode(as.character(Strategy),
                           "PAnCReAS (integer score, calibrated)" = "PAnCReAS (integer score)"),
         across(starts_with("NB_at_"), ~ sprintf("%.6f", .x))) %>%
  rename_with(~ sub("NB_at_(\\d+)pct", "NB at \\1%", .x))

cat("=== TABLE S5: DECISION CURVE ANALYSIS (1-YEAR RISK) ===\n")
print(as.data.frame(table_s5), row.names = FALSE)
write.csv(table_s5, file.path(OUT_DIR, "Table_S5_decision_curve.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table_S5_decision_curve.csv\n\n")


# =============================================================================
# STEP 15 — FINAL SCORING CARD WITH PERFORMANCE TABLE
# =============================================================================

print_score_card("AP-CANCER COX RISK SCORE (1-YEAR MODEL) — FINAL SCORING CARD")
cat(sprintf("NOTE: Achievable score range %d to %d; no floor applied.\n",
            SCORE_MIN_ACHIEVABLE, SCORE_MAX_ACHIEVABLE))
cat(sprintf("Bands: %s | %s | %s\n", BAND_LABELS["Low"], BAND_LABELS["Moderate"], BAND_LABELS["High"]))
cat("----------------------------------------------------------------\n\n")
cat("SCORE   1-YR CANCER RATE   NNS   RECOMMENDATION\n")
cat("----------------------------------------------------------------\n")
for (i in seq_len(nrow(cancer_by_score))) {
  cat(sprintf("%5d   %13.2f%%  %5d   %s\n",
              cancer_by_score$cox_score[i],
              cancer_by_score$cancer_rate[i],
              cancer_by_score$NNS[i],
              cancer_by_score$Recommendation[i]))
}
cat("\n")
cat("CUTOFF PERFORMANCE (1-YEAR EVENTS):\n")
cat("----------------------------------------------------------------\n")
thresh_df %>%
  filter(Sensitivity >= 0.60) %>%
  select(Cutoff, Sensitivity, Specificity,
         Cancers_caught, Cancers_missed, Pct_flagged) %>%
  print()
cat(sprintf("\n1-yr Harrell's C, Cox model (apparent):            %.3f\n", c_stat$concordance))
cat(sprintf("1-yr Harrell's C, Cox model (optimism-corrected):  %.3f\n", corrected_C))
cat(sprintf("Optimism (Cox model):                              %.3f\n", optimism_C))
cat(sprintf("1-yr Harrell's C, integer score (apparent):        %.3f\n", app_score_c["C"]))
cat(sprintf("1-yr Harrell's C, integer score (corrected):       %.3f\n", corr_C1y_score))
cat(sprintf("Calibration slope (shrinkage):                     %.3f\n", shrinkage))
cat(sprintf("Integer score AUROC (1-year events):               %.3f\n", as.numeric(auc(roc_score))))
cat("----------------------------------------------------------------\n")
cat("NOTE: Follow-up truncated at 1 year; cancers after 1 year censored at 1 year.\n")
cat("NOTE: Same-day/same-admission cancers assigned surv_time = 1 day.\n")
cat("NOTE: In-sample performance. External validation required.\n")
cat("NOTE: Denominator = cohort after exclusions.\n")
cat(sprintf("NOTE: Excluded: chronic pancreatitis (nv_k86 == 1), first AP after %s, prior cancer, cancer with no date.\n",
            format(FU_CUTOFF)))
cat(sprintf("NOTE: Risk bands (1-year model): %s, %s, %s.\n",
            BAND_LABELS["Low"], BAND_LABELS["Moderate"], BAND_LABELS["High"]))
cat("================================================================\n\n")


# =============================================================================
# STEP 16 — SAVE
# =============================================================================

saveRDS(cox_model,             file.path(OUT_DIR, "cox_model.rds"))
saveRDS(point_map,             file.path(OUT_DIR, "cox_point_map.rds"))
saveRDS(score_derivation,      file.path(OUT_DIR, "cox_score_derivation.rds"))
saveRDS(s0_table,              file.path(OUT_DIR, "cox_baseline_survival.rds"))
# Count-bearing working objects go to the internal folder (unsuppressed)
saveRDS(restricted_ci_by_band, file.path(INTERNAL_DIR, "cox_1yr_restricted_incidence_by_PAnCReAS_band.rds"))
saveRDS(multi_horizon,         file.path(INTERNAL_DIR, "cox_multihorizon_incidence_by_PAnCReAS_band.rds"))

# Everything the VIC/TAS scripts need to score and calibrate, in one place.
# S0 is UNCENTRED (checked in Step 10A): P(cancer by 1yr) = 1 - S0_1yr ^ exp(score * ref_unit)
score_bundle <- list(
  point_map     = point_map,
  ref_unit      = ref_unit,
  ref_label     = ref_label,
  coefficients  = coef(cox_model),
  S0_1yr        = S0_1y,             # unrounded
  S0_centered   = FALSE,
  s0_table      = s0_table,
  score_min     = SCORE_MIN_ACHIEVABLE,
  score_max     = SCORE_MAX_ACHIEVABLE,
  band_low_max  = BAND_LOW_MAX,
  band_high_min = BAND_HIGH_MIN,
  horizon_yrs   = HORIZON_YRS,
  census_date   = CENSUS_DATE,
  fu_cutoff     = FU_CUTOFF,
  flowchart     = flow_sant,
  risk_formula  = "1 - S0_1yr ^ exp(score * ref_unit)",
  created       = Sys.time()
)
saveRDS(score_bundle, file.path(OUT_DIR, "PAnCReAS_score_bundle_1yr.rds"))

cat("All outputs saved to:", OUT_DIR, "\n\n")
cat("Files saved:\n")
cat("  cox_model.rds\n")
cat("  cox_point_map.rds\n")
cat("  cox_score_derivation.rds / .csv\n")
cat("  cox_baseline_survival.rds / .csv\n")
cat("  cox_bootstrap_validation.csv\n")
cat("  cox_1yr_harrell_C_optimism_corrected.csv / .rds\n")
cat("  Table3_risk_bands_1yr.csv / Table3_footnote.txt\n")
cat("  Table_S3_score_derivation.csv\n")
cat("  Table_S4_diagnostic_performance.csv\n")
cat("  Table_S5_decision_curve.csv\n")
cat("  cox_schoenfeld_PH_test.csv\n")
cat("  cox_score_performance.csv\n")
cat("  risk_by_point_1yr.csv\n")
cat("  band_cancer_distribution.csv\n")
cat("  cox_cutoff_performance.csv\n")
cat("  cox_1yr_restricted_incidence_by_PAnCReAS_band.csv / .rds\n")
cat("  cox_multihorizon_incidence_by_PAnCReAS_band.csv / .rds\n")
cat("  Figure1_KM_incidence_by_PAnCReAS_band.png\n")
cat("  cox_km_overall.png\n")
cat("  cox_km_by_score.png\n")
cat("  cox_roc.png\n")
cat("  cox_score_cancer_rate.png\n")
cat("  cox_score_sensitivity_specificity.png\n")
cat("  cox_dca.png\n")
cat("  cox_dca_summary_table.csv\n")
cat("  cox_schoenfeld.png\n")
cat("  flowchart_SANT_1yr.csv\n")
cat("  PAnCReAS_score_bundle_1yr.rds  (point map, ref unit, uncentred S0(1yr), bands)\n")
cat(sprintf("\nUnsuppressed copies (do not request release): %s\n\n", INTERNAL_DIR))
cat("=== 1-YEAR-FIT COX SCORE ANALYSIS COMPLETE ===\n")




















# =============================================================================
# TIME AP -> CANCER BY PAnCReAS BAND, 1-YEAR MODEL
# =============================================================================
library(tidyverse)
OUT_DIR <- path.expand("~/Desktop/PAnCReAS_16-10")
if (!exists("point_map")) point_map <- readRDS(file.path(OUT_DIR, "cox_point_map.rds"))
HORIZON_YRS   <- 1
BAND_LOW_MAX  <- 1
BAND_HIGH_MIN <- 9
census_date   <- as.Date("2022-01-01")

if (!exists("build_ap_cohort"))
  stop("Run Steps 0-0D at the top of this script first (cohort builder, small-cell helpers).")

wk_dat <- build_ap_cohort(cmbspine_pancreatitis) %>%
  mutate(
    age_50_59      = as.integer(nv_k85_firstage >= 50 & nv_k85_firstage < 60),
    age_60_69      = as.integer(nv_k85_firstage >= 60 & nv_k85_firstage < 70),
    age_70plus     = as.integer(nv_k85_firstage >= 70),
    male           = as.integer(sant_nv_sex == min(sant_nv_sex, na.rm = TRUE)),
    biliary_ap     = as.integer(nv_k851 == 1),
    alcohol_ap     = as.integer(nv_k852 == 1),
    unspecified_ap = as.integer(nv_k859 == 1),
    other_ap       = as.integer(coalesce(nv_k85own, 0L) == 1 | coalesce(nv_k85other, 0L) == 1 |
                                  coalesce(nv_k858, 0L) == 1   | coalesce(nv_k853, 0L) == 1)
  ) %>%
  drop_na(age_50_59, age_60_69, age_70plus, male,
          biliary_ap, alcohol_ap, unspecified_ap, other_ap) %>%
  mutate(
    score =
      age_50_59 * point_map["age_50_59"] + age_60_69 * point_map["age_60_69"] +
      age_70plus * point_map["age_70plus"] + male * point_map["male"] +
      biliary_ap * point_map["biliary_ap"] + alcohol_ap * point_map["alcohol_ap"] +
      unspecified_ap * point_map["unspecified_ap"] + other_ap * point_map["other_ap"],
    band = factor(case_when(score <= BAND_LOW_MAX  ~ "Low",
                            score >= BAND_HIGH_MIN ~ "High",
                            TRUE                   ~ "Moderate"),
                  levels = c("Low", "Moderate", "High"))
  )

cat("Cohort N:", nrow(wk_dat),
    "| cancers within 1 year:", sum(wk_dat$event),
    "| (full follow-up:", sum(wk_dat$event_full), ")\n")
if (exists("dat_cox")) {
  cat("Matches main pipeline (N / 1-yr events):",
      nrow(wk_dat) == nrow(dat_cox) && sum(wk_dat$event) == sum(dat_cox$event), "\n")
}
cat("\n")

# Two case sets:
#   1-year cancers  (event == 1)
#   ALL cancers     (event_full == 1, full follow-up)
cases_1yr <- wk_dat %>% filter(event == 1)      %>% select(band, diffweeks = k85_panctcancer_diffweeks)
cases_all <- wk_dat %>% filter(event_full == 1) %>% select(band, diffweeks = k85_panctcancer_diffweeks)

# Small-cell rule: counts 1-9 -> "<10"; groups with <10 cancers get no
# median/IQR (it would describe fewer than 10 people).
n_pct <- function(x) {
  k <- sum(x)
  if (sdc_small(k)) "<10" else sprintf("%d (%.1f%%)", k, mean(x) * 100)
}

summ <- function(d) d %>% summarise(
  n_cancers  = sdc_lab(n()),
  median_wk  = if (sdc_small(n())) NA_real_ else median(diffweeks),
  median_IQR = if (sdc_small(n())) "S" else
    sprintf("%.0f (%.0f–%.0f)", median(diffweeks),
            quantile(diffweeks, 0.25), quantile(diffweeks, 0.75)),
  median_mo  = if (sdc_small(n())) NA_real_ else round(median(diffweeks) * 7 / 30.4375, 1),
  by_4wk     = n_pct(diffweeks <= 4),
  by_12wk    = n_pct(diffweeks <= 12),
  by_26wk    = n_pct(diffweeks <= 26),
  by_36wk    = n_pct(diffweeks <= 36),
  by_52wk    = n_pct(diffweeks <= 52),
  after_52wk = n_pct(diffweeks > 52),
  .groups = "drop")

by_band <- function(cases) bind_rows(
  cases %>% group_by(band) %>% summ() %>% mutate(band = as.character(band)),
  cases %>% summ() %>% mutate(band = "Overall")
)

ttc_1yr <- by_band(cases_1yr)
ttc_all <- by_band(cases_all)

cat(sprintf("=== Time AP -> cancer: denominator = 1-YEAR cancers (N = %d) ===\n", nrow(cases_1yr)))
print(ttc_1yr %>% select(-after_52wk), width = Inf)
cat(sprintf("\n=== Time AP -> cancer: denominator = ALL cancers, full follow-up (N = %d) ===\n",
            nrow(cases_all)))
print(ttc_all, width = Inf)

# Combined manuscript table: each timepoint shown against both denominators
time_to_cancer <- ttc_1yr %>%
  select(band, n_1yr = n_cancers, median_IQR_1yr = median_IQR,
         by_4wk_1yr = by_4wk, by_12wk_1yr = by_12wk, by_26wk_1yr = by_26wk,
         by_36wk_1yr = by_36wk) %>%
  full_join(ttc_all %>%
              select(band, n_all = n_cancers, median_IQR_all = median_IQR,
                     by_4wk_all = by_4wk, by_12wk_all = by_12wk, by_26wk_all = by_26wk,
                     by_36wk_all = by_36wk, by_52wk_all = by_52wk, after_52wk_all = after_52wk),
            by = "band") %>%
  mutate(band = factor(band, levels = c("Low", "Moderate", "High", "Overall"))) %>%
  arrange(band) %>%
  relocate(band, n_1yr, n_all)

cat("\n=== Combined time-to-cancer table (1-year and all-cancer denominators) ===\n")
print(as.data.frame(time_to_cancer), row.names = FALSE)

kw_1yr <- kruskal.test(diffweeks ~ band, data = cases_1yr)
kw_all <- kruskal.test(diffweeks ~ band, data = cases_all)
cat(sprintf("\nKruskal-Wallis across bands: 1-year cancers p = %s | all cancers p = %s\n",
            format.pval(kw_1yr$p.value, digits = 2, eps = 0.001),
            format.pval(kw_all$p.value, digits = 2, eps = 0.001)))

write.csv(time_to_cancer, file.path(OUT_DIR, "time_to_cancer_weeks_by_band_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(ttc_1yr, file.path(OUT_DIR, "time_to_cancer_weeks_by_band_denom_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(ttc_all, file.path(OUT_DIR, "time_to_cancer_weeks_by_band_denom_all.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: time_to_cancer_weeks_by_band_1yr.csv (combined),\n")
cat("       time_to_cancer_weeks_by_band_denom_1yr.csv, time_to_cancer_weeks_by_band_denom_all.csv\n")












# =============================================================================
# TABLE 1 — SANT DEMOGRAPHICS: AP COHORT AND 1-YEAR CANCER COHORT
# Same cohort and event rules as Step 1. Columns:
#   SANT AP                (all AP patients)
#   SANT Cancer, 1-year    (cancer within 1 year of first AP = model outcome)
#   SANT Cancer, all       (cancer at any time after AP; for reference)
# Aetiology codes are not mutually exclusive (an episode can carry >1 code).
# Counts of 1-9 are shown as "<10" (data custodian rule); any secondary
# suppression ("S") to prevent back-calculation must be checked by hand.
# =============================================================================

library(tidyverse)

OUT_DIR     <- path.expand("~/Desktop/PAnCReAS_16-10")
HORIZON_YRS <- 1
census_date <- as.Date("2022-01-01")

if (!exists("build_ap_cohort"))
  stop("Run Steps 0-0D at the top of this script first (cohort builder, small-cell helpers).")

t1_dat <- build_ap_cohort(cmbspine_pancreatitis) %>%
  mutate(
    fu_years  = surv_years_full,
    event_1yr = event,
    age       = nv_k85_firstage,
    sex       = case_when(is.na(sant_nv_sex) ~ NA_character_,
                          sant_nv_sex == min(sant_nv_sex, na.rm = TRUE) ~ "Male",
                          TRUE ~ "Female"),
    biliary     = coalesce(nv_k851 == 1, FALSE),
    alcohol     = coalesce(nv_k852 == 1, FALSE),
    unspecified = coalesce(nv_k859 == 1, FALSE),
    other       = coalesce(nv_k85own, 0L) == 1 | coalesce(nv_k85other, 0L) == 1 |
      coalesce(nv_k858, 0L) == 1   | coalesce(nv_k853, 0L) == 1
  )

# ---- helpers ----------------------------------------------------------------
n_pct <- function(x, denom) {
  k <- sum(x, na.rm = TRUE)
  if (k > 0 && k < 10) return("<10")
  sprintf("%s (%.1f%%)", format(k, big.mark = ","), 100 * k / denom)
}
med_iqr <- function(x, dp = 1) {
  q <- quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf(paste0("%.", dp, "f (%.", dp, "f–%.", dp, "f)"), q[1], q[2], q[3])
}

t1_column <- function(d, fu_label_value) {
  N <- nrow(d)
  c(
    `Cancer events within 1 year, n (%)`  = n_pct(d$event_1yr == 1, N),
    `Cancer events, any time, n (%)`      = n_pct(d$event_full == 1, N),
    `Median age at first AP, years (IQR)` = med_iqr(d$age),
    `Age category`                        = "",
    `  <50 years`                         = n_pct(d$age < 50, N),
    `  50–64 years`                  = n_pct(d$age >= 50 & d$age < 65, N),
    `  ≥65 years`                    = n_pct(d$age >= 65, N),
    `Sex`                                 = "",
    `  Male`                              = n_pct(d$sex == "Male", N),
    `  Female`                            = n_pct(d$sex == "Female", N),
    `Aetiology of AP`                     = "",
    `  Biliary`                           = n_pct(d$biliary, N),
    `  Alcohol-induced`                   = n_pct(d$alcohol, N),
    `  Unspecified`                       = n_pct(d$unspecified, N),
    `  Other`                             = n_pct(d$other, N)
  )
}

ap      <- t1_dat
ca_1yr  <- t1_dat %>% filter(event_1yr == 1)
ca_all  <- t1_dat %>% filter(event_full == 1)

# (Median follow-up row removed: follow-up is capped at 1 year by design.)
# from first AP to cancer diagnosis.
col_ap     <- t1_column(ap,     med_iqr(ap$fu_years))
col_ca_1yr <- t1_column(ca_1yr, med_iqr(ca_1yr$fu_years, dp = 2))
col_ca_all <- t1_column(ca_all, med_iqr(ca_all$fu_years))

table1 <- tibble(
  Characteristic = names(col_ap),
  !!sprintf("SANT AP (N = %s)",                 format(nrow(ap),     big.mark = ",")) := unname(col_ap),
  !!sprintf("SANT Cancer, 1-year (N = %s)",     format(nrow(ca_1yr), big.mark = ",")) := unname(col_ca_1yr),
  !!sprintf("SANT Cancer, any time (N = %s)",   format(nrow(ca_all), big.mark = ",")) := unname(col_ca_all)
)

cat("=== TABLE 1: SANT AP AND CANCER COHORTS (1-YEAR) ===\n")
print(as.data.frame(table1), row.names = FALSE, right = FALSE)

# Missing data, for the footnote (not suppressed here - check before reporting)
cat(sprintf("\nMissing: sex n = %d | age n = %d (AP cohort)\n",
            sum(is.na(ap$sex)), sum(is.na(ap$age))))
cat(sprintf("AP patients with < 1 year follow-up and no 1-year cancer: %d\n",
            sum(ap$fu_years < HORIZON_YRS & ap$event_1yr == 0)))
cat("Check: SANT AP N should equal nrow(dat_cox); 1-year cancers should equal sum(dat_cox$event).\n")
if (exists("dat_cox"))
  cat(sprintf("  dat_cox N = %d, 1-year events = %d -> match: %s\n",
              nrow(dat_cox), sum(dat_cox$event),
              nrow(dat_cox) == nrow(ap) && sum(dat_cox$event) == nrow(ca_1yr)))

write.csv(table1, file.path(OUT_DIR, "Table1_SANT_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table1_SANT_1yr.csv\n")