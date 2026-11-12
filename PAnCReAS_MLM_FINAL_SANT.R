# =============================================================================
# MLM SANT (Development Dataset) — 1-YEAR OUTCOME (V2)
# =============================================================================

# =============================================================================
# 0. PACKAGE INSTALLATION AND LOADING
# =============================================================================

required_packages <- c(
  "tidyverse", "lubridate", "janitor",
  "caret", "glmnet", "randomForest", "xgboost",
  "themis",
  "ROSE",
  "pROC", "PRROC", "MLmetrics", "yardstick",
  "CalibrationCurves", "gbm",
  "SHAPforxgboost", "shapviz",
  "mice", "naniar", "visdat",
  "ggplot2", "patchwork", "scales", "viridis",
  "ggpubr", "gridExtra", "RColorBrewer",
  "survival", "survminer", "cmprsk", "officer", "flextable",
  "here", "tictoc"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[, "Package"])
]
if (length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

invisible(lapply(required_packages, library, character.only = TRUE))

set.seed(42)

cat("=== Packages loaded successfully ===\n\n")

OUT_DIR <- path.expand("~/Desktop/LASSO_19.10")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Output folder:", OUT_DIR, "\n\n")

HORIZON_YRS     <- 1        # outcome horizon (years)



# =============================================================================
# STEP 0C — COHORT DEFINITION WITH EXCLUSIONS (single source of truth)
# Used for every cohort built in this script, so all outputs share the same
# patients. Exclusions, applied in this order:
# Starts from the WHOLE cohort (no nv_k85 filter); patients without an AP
# episode drop out at the missing first-AP date step.
#   1. K86 == 1 (nv_k86) - removed
#      from the cohort entirely (neither event nor control)
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
# 1. DATA PREPARATION
# =============================================================================
# Entry criterion:  whole cohort with a first AP on/before FU_CUTOFF
#                   (CENSUS_DATE - 365 days): everyone can complete 1 year
# Exclusions:       chronic pancreatitis (nv_k86 == 1), prior pancreatic
#                   cancer, cancer with no diagnosis date, < 1 year possible
#                   follow-up (Step 0C, same as Cox script)
# Event definition: k85_panctcancer == 1 AND diffweeks >= 0
#                   AND diagnosed within HORIZON_YRS of first AP
# Controls:         everyone else (deaths within 1 year without cancer
#                   remain negatives)
# =============================================================================

census_date <- CENSUS_DATE
cat("Census date (end of follow-up):", format(census_date), "\n")
cat("Last eligible first-AP date:   ", format(FU_CUTOFF), "\n\n")

ml_base <- build_ap_cohort(cmbspine_pancreatitis)
flow_ml <- attr(ml_base, "flow")

dat_ml <- ml_base %>%
  mutate(
    age_cat = case_when(
      nv_k85_firstage < 50                          ~ 0,
      nv_k85_firstage >= 50 & nv_k85_firstage < 60 ~ 1,
      nv_k85_firstage >= 60 & nv_k85_firstage < 70 ~ 2,
      nv_k85_firstage >= 70                         ~ 3
    ),
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
  )

cat("=== COHORT SUMMARY ===\n")
print(as.data.frame(flow_ml), row.names = FALSE)
cat("\nTotal N (after exclusions):", nrow(dat_ml), "\n")
cat("Cancer events (all):       ", sum(dat_ml$event_full), "\n")
cat("Cancer events (1-year):    ", sum(dat_ml$event), " <- model outcome\n")
cat("Controls:                  ", sum(dat_ml$event == 0), "\n")
cat("Event rate:                ", round(mean(dat_ml$event) * 100, 2), "%\n\n")

dat_model <- dat_ml %>%
  mutate(
    outcome = factor(event,
                     levels = c(0, 1),
                     labels = c("No_Cancer", "Cancer"))
  ) %>%
  select(
    outcome,
    nv_k85_firstage,
    male,
    biliary_ap,
    alcohol_ap,
    unspecified_ap,
    other_ap
  )

cat("Feature matrix dimensions:", nrow(dat_model), "x",
    ncol(dat_model) - 1, "(features)\n\n")


# =============================================================================
# 2. MISSING DATA AUDIT
# =============================================================================

cat("=== MISSING DATA AUDIT ===\n\n")

missing_summary <- dat_model %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(),
               names_to  = "variable",
               values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(dat_model) * 100, 2)) %>%
  arrange(desc(pct_missing))

print(missing_summary)

p_miss <- vis_miss(dat_model, warn_large_data = FALSE) +
  labs(title = "Missing Data Pattern — Analytical Dataset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "missing_data_pattern.png"), p_miss, width = 10, height = 6, dpi = 300)

mar_table <- table(
  miss_age = is.na(dat_model$nv_k85_firstage),
  cancer   = dat_model$outcome == "Cancer"
)

cat("\n--- MAR Check ---\n")
print(mar_table)

if (nrow(mar_table) < 2 || ncol(mar_table) < 2) {
  cat("=> Skipping Fisher's test — no missing values in nv_k85_firstage",
      "(or all values missing). MAR assumption untestable.\n\n")
  mar_p <- NA
} else {
  mar_test <- fisher.test(mar_table)
  mar_p    <- mar_test$p.value
  cat(sprintf("Fisher's exact test p-value: %.4f\n", mar_p))
  if (mar_p < 0.05) {
    cat("=> WARNING: Possible MNAR — imputation results should be interpreted with caution.\n\n")
  } else {
    cat("=> No significant association (MAR assumption plausible).\n\n")
  }
}

max_miss <- max(missing_summary$pct_missing, na.rm = TRUE)
cat("Maximum missingness:", max_miss, "%\n")

if (max_miss < 5) {
  cat("=> Proceeding with COMPLETE CASE ANALYSIS.\n\n")
  dat_analysis <- dat_model %>% drop_na()
  imputed_flag <- FALSE
} else {
  cat("=> Applying MULTIPLE IMPUTATION (mice, m = 5, pmm).\n\n")
  imputed_flag <- TRUE
  imp <- mice(dat_model, m = 5, method = "pmm", seed = 42, printFlag = FALSE)
  dat_analysis <- complete(imp, 1)
  cat("Imputation complete. Using m = 1 for primary analysis.\n\n")
}

cat("Final analytical N:", nrow(dat_analysis), "\n")

flow_ml <- bind_rows(
  flow_ml,
  tibble(Step = c(if (imputed_flag) "Missing predictors imputed (mice); none excluded"
                  else "Excluded: missing predictor data (complete-case)",
                  "Analytic cohort",
                  "  Pancreatic cancer within 1 year (model outcome)"),
         N = c(nrow(dat_model) - nrow(dat_analysis), nrow(dat_analysis),
               sum(dat_analysis$outcome == "Cancer")))
)
write_release(flow_ml, list(N = character(0)), "flowchart_SANT_XGB_1yr.csv")
cat("Saved: flowchart_SANT_XGB_1yr.csv\n")
cat("Outcome distribution:\n")
print(table(dat_analysis$outcome))
cat("\n")


# =============================================================================
# 3. DATA SPLITTING — RANDOM 80/20 SPLIT (INTERNAL VALIDATION ONLY)
# =============================================================================

cat("=== DATA SPLITTING: Random 80/20 (internal validation) ===\n\n")

set.seed(42)
train_idx <- createDataPartition(dat_analysis$outcome,
                                 p     = 0.80,
                                 list  = FALSE,
                                 times = 1)

dat_train_val <- dat_analysis[ train_idx, ]
dat_test_val  <- dat_analysis[-train_idx, ]

cat(sprintf("Internal train (80%%): N = %d | Cancer = %d | No Cancer = %d\n",
            nrow(dat_train_val),
            sum(dat_train_val$outcome == "Cancer"),
            sum(dat_train_val$outcome == "No_Cancer")))
cat(sprintf("Internal test  (20%%): N = %d | Cancer = %d | No Cancer = %d\n\n",
            nrow(dat_test_val),
            sum(dat_test_val$outcome == "Cancer"),
            sum(dat_test_val$outcome == "No_Cancer")))

X_train_val <- dat_train_val %>% select(-outcome)
y_train_val <- dat_train_val$outcome
X_test_val  <- dat_test_val  %>% select(-outcome)
y_test_val  <- dat_test_val$outcome

y_train_val_bin <- ifelse(y_train_val == "Cancer", 1, 0)
y_test_val_bin  <- ifelse(y_test_val  == "Cancer", 1, 0)

n_pos            <- sum(y_train_val == "Cancer")
n_neg            <- sum(y_train_val == "No_Cancer")
scale_pos_weight <- n_neg / n_pos
cat("Internal train class ratio (neg/pos):", round(scale_pos_weight, 2), "\n\n")


# =============================================================================
# 4. CROSS-VALIDATION SETUP
# =============================================================================

cv_control <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling        = NULL
)

cat("=== 5-fold stratified cross-validation configured ===\n\n")


# =============================================================================
# 5. MODEL 1: LASSO LOGISTIC REGRESSION
# =============================================================================

cat("=== MODEL 1: LASSO Logistic Regression ===\n\n")
tic("LASSO")

X_train_val_mat <- model.matrix(~ . - 1, data = X_train_val)
X_test_val_mat  <- model.matrix(~ . - 1, data = X_test_val)

case_weights_train <- ifelse(y_train_val == "Cancer", n_neg / n_pos, 1)

set.seed(42)
cv_lasso <- cv.glmnet(
  x            = X_train_val_mat,
  y            = y_train_val_bin,
  family       = "binomial",
  alpha        = 1,
  weights      = case_weights_train,
  nfolds       = 5,
  type.measure = "auc",
  standardize  = TRUE
)

lambda_min  <- cv_lasso$lambda.min
lambda_1se  <- cv_lasso$lambda.1se
lambda_best <- lambda_min
cat("Lambda (1SE):", round(lambda_1se, 6),
    "| Lambda (min):", round(lambda_min, 6), "\n")
cat("Using lambda.min for final model.\n")

lasso_val <- glmnet(
  x           = X_train_val_mat,
  y           = y_train_val_bin,
  family      = "binomial",
  alpha       = 1,
  lambda      = lambda_best,
  weights     = case_weights_train,
  standardize = TRUE
)

lasso_coef <- coef(lasso_val) %>%
  as.matrix() %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  rename(Coefficient = s0) %>%
  filter(Coefficient != 0, Variable != "(Intercept)") %>%
  arrange(desc(abs(Coefficient)))

cat("Non-zero LASSO coefficients (lambda.min):\n")
print(lasso_coef)

lasso_prob_train_val <- predict(lasso_val, X_train_val_mat, type = "response")[, 1]
lasso_prob_test_val  <- predict(lasso_val, X_test_val_mat,  type = "response")[, 1]

toc(); cat("\n")


# =============================================================================
# 6. MODEL 2: RANDOM FOREST
# =============================================================================

cat("=== MODEL 2: Random Forest ===\n\n")
tic("Random Forest")

rf_classwt <- c("No_Cancer" = 1, "Cancer" = scale_pos_weight)

set.seed(42)
rf_tune <- train(
  x         = X_train_val,
  y         = y_train_val,
  method    = "rf",
  metric    = "ROC",
  trControl = cv_control,
  tuneGrid  = expand.grid(mtry = c(2, 3, 4, 5)),
  ntree     = 500,
  classwt   = rf_classwt
)

cat("Best RF mtry:", rf_tune$bestTune$mtry, "\n")

set.seed(42)
rf_val <- randomForest(
  x          = X_train_val,
  y          = y_train_val,
  ntree      = 500,
  mtry       = rf_tune$bestTune$mtry,
  classwt    = rf_classwt,
  importance = TRUE
)

rf_prob_train_val <- predict(rf_val, X_train_val, type = "prob")[, "Cancer"]
rf_prob_test_val  <- predict(rf_val, X_test_val,  type = "prob")[, "Cancer"]

toc(); cat("\n")


# =============================================================================
# 7. MODEL 3: XGBOOST — AUPRC-OPTIMISED
# =============================================================================

cat("=== MODEL 3: XGBoost (AUPRC-optimised) ===\n\n")
tic("XGBoost")

dtrain_val <- xgb.DMatrix(data = X_train_val_mat, label = y_train_val_bin)
dtest_val  <- xgb.DMatrix(data = X_test_val_mat,  label = y_test_val_bin)

xgb_grid <- expand.grid(
  nrounds          = c(200, 300, 500),
  max_depth        = c(2, 3, 4, 5, 6, 7, 8),
  eta              = c(0.01, 0.05, 0.1),
  gamma            = c(0, 1, 5),
  colsample_bytree = 0.8,
  min_child_weight = c(1, 5, 10, 20),
  subsample        = 0.8,
  spw_multiplier   = c(1.0, 1.5, 2.0, 3.0),
  stringsAsFactors = FALSE
)

set.seed(42)
xgb_grid_sample <- xgb_grid[sample(nrow(xgb_grid), 40), ]

cat("Tuning XGBoost via xgb.cv() (5-fold, stratified, metric = AUPRC)...\n")

xgb_cv_results <- purrr::pmap_dfr(xgb_grid_sample, function(
    nrounds, max_depth, eta, gamma,
    colsample_bytree, min_child_weight, subsample, spw_multiplier) {
  
  spw <- scale_pos_weight * spw_multiplier
  
  set.seed(42)
  cv <- xgb.cv(
    params = list(
      objective        = "binary:logistic",
      eval_metric      = "aucpr",
      eta              = eta,
      max_depth        = max_depth,
      gamma            = gamma,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      subsample        = subsample,
      scale_pos_weight = spw,
      nthread          = parallel::detectCores() - 1
    ),
    data                  = dtrain_val,
    nrounds               = nrounds,
    nfold                 = 5,
    stratified            = TRUE,
    verbose               = FALSE,
    showsd                = FALSE,
    early_stopping_rounds = 20
  )
  
  best_round <- if (!is.null(cv$best_iteration) && !is.na(cv$best_iteration)) {
    cv$best_iteration
  } else {
    which.max(cv$evaluation_log$test_aucpr_mean)
  }
  best_auprc <- max(cv$evaluation_log$test_aucpr_mean, na.rm = TRUE)
  
  tibble(
    nrounds          = best_round,
    max_depth        = max_depth,
    eta              = eta,
    gamma            = gamma,
    colsample_bytree = colsample_bytree,
    min_child_weight = min_child_weight,
    subsample        = subsample,
    spw_multiplier   = spw_multiplier,
    scale_pos_weight = spw,
    cv_auprc         = best_auprc
  )
})

best_params <- xgb_cv_results %>%
  slice_max(cv_auprc, n = 1, with_ties = FALSE)

cat("Best XGBoost parameters:\n")
print(best_params)
cat(sprintf("Best CV AUPRC: %.4f\n", best_params$cv_auprc))
cat(sprintf("scale_pos_weight used: %.2f (%.1fx natural ratio)\n",
            best_params$scale_pos_weight, best_params$spw_multiplier))

set.seed(42)
xgb_val <- xgb.train(
  params = list(
    objective        = "binary:logistic",
    eval_metric      = "aucpr",
    eta              = best_params$eta,
    max_depth        = best_params$max_depth,
    gamma            = best_params$gamma,
    colsample_bytree = best_params$colsample_bytree,
    min_child_weight = best_params$min_child_weight,
    subsample        = best_params$subsample,
    scale_pos_weight = best_params$scale_pos_weight,
    nthread          = parallel::detectCores() - 1
  ),
  data    = dtrain_val,
  nrounds = best_params$nrounds,
  evals   = list(train = dtrain_val, test = dtest_val),
  verbose = 0
)

xgb_prob_train_val <- predict(xgb_val, dtrain_val)
xgb_prob_test_val  <- predict(xgb_val, dtest_val)

toc(); cat("\n")


# =============================================================================
# 8. POST-HOC CALIBRATION — PLATT SCALING
# =============================================================================

cat("=== POST-HOC CALIBRATION: Platt Scaling ===\n\n")

calibrate_platt <- function(prob_train, prob_test, y_train_bin, model_name) {
  cal_df   <- data.frame(prob = prob_train, y = y_train_bin)
  platt_m  <- glm(y ~ prob, data = cal_df, family = binomial)
  prob_cal <- predict(platt_m,
                      newdata = data.frame(prob = prob_test),
                      type    = "response")
  cat(model_name, "- Platt scaling applied.\n")
  list(model = platt_m, prob_cal = prob_cal)
}

platt_lasso_val <- calibrate_platt(lasso_prob_train_val, lasso_prob_test_val,
                                   y_train_val_bin, "LASSO")
platt_rf_val    <- calibrate_platt(rf_prob_train_val,    rf_prob_test_val,
                                   y_train_val_bin, "Random Forest")
platt_xgb_val   <- calibrate_platt(xgb_prob_train_val,  xgb_prob_test_val,
                                   y_train_val_bin, "XGBoost")

lasso_prob_test_val_cal <- platt_lasso_val$prob_cal
rf_prob_test_val_cal    <- platt_rf_val$prob_cal
xgb_prob_test_val_cal   <- platt_xgb_val$prob_cal

cat("\n")


# =============================================================================
# 9. THRESHOLD OPTIMISATION (Precision-Recall Curve on validation set)
# =============================================================================

cat("=== THRESHOLD OPTIMISATION via Precision-Recall ===\n\n")

optimise_threshold <- function(probs, labels_bin, model_name) {
  thresholds <- seq(0.01, 0.99, by = 0.01)
  results <- map_dfr(thresholds, function(thresh) {
    pred_class <- ifelse(probs >= thresh, 1, 0)
    tp <- sum(pred_class == 1 & labels_bin == 1)
    fp <- sum(pred_class == 1 & labels_bin == 0)
    fn <- sum(pred_class == 0 & labels_bin == 1)
    tn <- sum(pred_class == 0 & labels_bin == 0)
    precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
    recall    <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
    f1        <- ifelse((precision + recall) == 0, 0,
                        2 * precision * recall / (precision + recall))
    tibble(threshold = thresh, precision = precision,
           recall = recall, f1 = f1, tp = tp, fp = fp, fn = fn, tn = tn)
  })
  best <- results %>% slice_max(f1, n = 1, with_ties = FALSE)
  cat(model_name, "- Optimal threshold:", best$threshold,
      "| F1:", round(best$f1, 4),
      "| Sensitivity:", round(best$recall, 4),
      "| Specificity:", round(best$tn / (best$tn + best$fp), 4),
      "| PPV:", round(best$precision, 4), "\n")
  list(threshold = best$threshold, results = results)
}

thresh_lasso <- optimise_threshold(lasso_prob_train_val, y_train_val_bin, "LASSO")
thresh_rf    <- optimise_threshold(rf_prob_train_val,    y_train_val_bin, "Random Forest")
thresh_xgb   <- optimise_threshold(xgb_prob_train_val,  y_train_val_bin, "XGBoost")
cat("\n")


# =============================================================================
# 10. PERFORMANCE EVALUATION FUNCTION
# =============================================================================

evaluate_model <- function(probs, labels_bin, threshold, model_name,
                           dataset_label = "Validation") {
  pred_class <- factor(ifelse(probs >= threshold, "Cancer", "No_Cancer"),
                       levels = c("No_Cancer", "Cancer"))
  true_class <- factor(ifelse(labels_bin == 1, "Cancer", "No_Cancer"),
                       levels = c("No_Cancer", "Cancer"))
  
  cm      <- confusionMatrix(pred_class, true_class, positive = "Cancer")
  roc_obj <- roc(labels_bin, probs, quiet = TRUE)
  auroc   <- as.numeric(auc(roc_obj))
  pr_obj  <- pr.curve(scores.class0 = probs[labels_bin == 1],
                      scores.class1 = probs[labels_bin == 0],
                      curve = TRUE)
  auprc   <- pr_obj$auc.integral
  brier   <- mean((probs - labels_bin)^2)
  
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  ppv  <- cm$byClass["Pos Pred Value"]
  f1   <- cm$byClass["F1"]
  
  cat("\n---", model_name, "|", dataset_label, "---\n")
  cat(sprintf("  AUROC:        %.4f\n", auroc))
  cat(sprintf("  AUPRC:        %.4f\n", auprc))
  cat(sprintf("  Sensitivity:  %.4f\n", sens))
  cat(sprintf("  Specificity:  %.4f\n", spec))
  cat(sprintf("  F1 Score:     %.4f\n", f1))
  cat(sprintf("  PPV:          %.4f\n", ppv))
  cat(sprintf("  Brier Score:  %.4f\n", brier))
  
  list(model = model_name, dataset = dataset_label,
       auroc = auroc, auprc = auprc,
       sensitivity = sens, specificity = spec,
       f1 = f1, ppv = ppv, brier = brier,
       roc_obj = roc_obj, pr_obj = pr_obj, cm = cm,
       pred_class = pred_class, true_class = true_class,
       probs = probs, labels_bin = labels_bin)
}


# =============================================================================
# 11. RUN EVALUATIONS — INTERNAL VALIDATION SET
# =============================================================================

cat("\n=== MODEL PERFORMANCE — INTERNAL VALIDATION SET (20%) ===\n")

results_lasso <- evaluate_model(lasso_prob_test_val, y_test_val_bin,
                                thresh_lasso$threshold, "LASSO")
results_rf    <- evaluate_model(rf_prob_test_val,    y_test_val_bin,
                                thresh_rf$threshold,    "Random Forest")
results_xgb   <- evaluate_model(xgb_prob_test_val,  y_test_val_bin,
                                thresh_xgb$threshold,   "XGBoost")

brier_cal_lasso <- mean((lasso_prob_test_val_cal - y_test_val_bin)^2)
brier_cal_rf    <- mean((rf_prob_test_val_cal    - y_test_val_bin)^2)
brier_cal_xgb   <- mean((xgb_prob_test_val_cal  - y_test_val_bin)^2)

brier_null        <- mean(y_test_val_bin) * (1 - mean(y_test_val_bin))
brier_skill_lasso <- 1 - (brier_cal_lasso / brier_null)
brier_skill_rf    <- 1 - (brier_cal_rf    / brier_null)
brier_skill_xgb   <- 1 - (brier_cal_xgb  / brier_null)

cat("\n--- Calibrated Brier Scores & Skill Scores ---\n")
cat(sprintf("  Null model Brier:  %.4f\n", brier_null))
cat(sprintf("  LASSO         — Brier: %.4f  |  Skill: %.4f\n", brier_cal_lasso, brier_skill_lasso))
cat(sprintf("  Random Forest — Brier: %.4f  |  Skill: %.4f\n", brier_cal_rf,    brier_skill_rf))
cat(sprintf("  XGBoost       — Brier: %.4f  |  Skill: %.4f\n", brier_cal_xgb,   brier_skill_xgb))

performance_table <- bind_rows(
  tibble(Model = "LASSO",
         AUROC        = results_lasso$auroc,
         AUPRC        = results_lasso$auprc,
         Sensitivity  = results_lasso$sensitivity,
         Specificity  = results_lasso$specificity,
         F1           = results_lasso$f1,
         PPV          = results_lasso$ppv,
         Brier_raw    = results_lasso$brier,
         Brier_cal    = brier_cal_lasso,
         Brier_skill  = brier_skill_lasso),
  tibble(Model = "Random Forest",
         AUROC        = results_rf$auroc,
         AUPRC        = results_rf$auprc,
         Sensitivity  = results_rf$sensitivity,
         Specificity  = results_rf$specificity,
         F1           = results_rf$f1,
         PPV          = results_rf$ppv,
         Brier_raw    = results_rf$brier,
         Brier_cal    = brier_cal_rf,
         Brier_skill  = brier_skill_rf),
  tibble(Model = "XGBoost",
         AUROC        = results_xgb$auroc,
         AUPRC        = results_xgb$auprc,
         Sensitivity  = results_xgb$sensitivity,
         Specificity  = results_xgb$specificity,
         F1           = results_xgb$f1,
         PPV          = results_xgb$ppv,
         Brier_raw    = results_xgb$brier,
         Brier_cal    = brier_cal_xgb,
         Brier_skill  = brier_skill_xgb)
) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

cat("\n=== PERFORMANCE SUMMARY TABLE (Internal Validation) ===\n")
print(performance_table)

best_model_name <- performance_table %>%
  slice_max(AUROC, n = 1, with_ties = FALSE) %>%
  pull(Model)
cat("\n=> Best model (internal validation):", best_model_name, "\n\n")


# =============================================================================
# 12. VISUALISATIONS (INTERNAL VALIDATION)
# =============================================================================

model_colours <- c("LASSO" = "#2166AC",
                   "Random Forest" = "#4DAF4A",
                   "XGBoost" = "#E41A1C")

plot_roc <- function(results_list, title = "ROC Curves — Internal Validation") {
  roc_data <- map_dfr(results_list, function(res) {
    data.frame(
      specificity = res$roc_obj$specificities,
      sensitivity = res$roc_obj$sensitivities,
      model       = paste0(res$model, " (AUC = ", round(res$auroc, 3), ")")
    )
  })
  ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, colour = model)) +
    geom_line(linewidth = 1.2) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = setNames(
      c("#2166AC", "#4DAF4A", "#E41A1C"), unique(roc_data$model)
    )) +
    labs(title  = title,
         x      = "1 - Specificity (False Positive Rate)",
         y      = "Sensitivity (True Positive Rate)",
         colour = "Model") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

p_roc <- plot_roc(list(results_lasso, results_rf, results_xgb))
ggsave(file.path(OUT_DIR, "roc_curves.png"), p_roc, width = 8, height = 7, dpi = 300)

plot_prc <- function(results_list, title = "Precision-Recall Curves — Internal Validation") {
  prc_data <- map_dfr(results_list, function(res) {
    data.frame(recall    = res$pr_obj$curve[, 1],
               precision = res$pr_obj$curve[, 2],
               model     = paste0(res$model, " (AUPRC = ", round(res$auprc, 3), ")"))
  })
  prevalence <- mean(results_lasso$labels_bin)
  ggplot(prc_data, aes(x = recall, y = precision, colour = model)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = prevalence, linetype = "dashed",
               colour = "grey50", linewidth = 0.8) +
    annotate("text", x = 0.8, y = prevalence + 0.005,
             label = paste0("Baseline = ", round(prevalence, 3)),
             size = 3.5, colour = "grey40") +
    scale_colour_manual(values = setNames(
      c("#2166AC", "#4DAF4A", "#E41A1C"), unique(prc_data$model)
    )) +
    labs(title  = title, x = "Recall (Sensitivity)",
         y = "Precision (PPV)", colour = "Model") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

p_prc <- plot_prc(list(results_lasso, results_rf, results_xgb))
ggsave(file.path(OUT_DIR, "precision_recall_curves.png"), p_prc, width = 8, height = 7, dpi = 300)

plot_cm <- function(result) {
  cm_table <- as.data.frame(result$cm$table) %>%
    group_by(Reference) %>%
    mutate(lab = if (any(sdc_small(Freq))) ifelse(sdc_small(Freq), "<10", "S")
           else format(Freq, big.mark = ",")) %>%
    ungroup()
  ggplot(cm_table, aes(x = Reference, y = Prediction, fill = Freq)) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = lab), size = 8, fontface = "bold") +
    scale_fill_gradient(low = "#DEEBF7", high = "#2166AC", name = "Count") +
    labs(title = paste("Confusion Matrix —", result$model),
         x = "True Class", y = "Predicted Class") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

p_cm_all <- plot_cm(results_lasso) + plot_cm(results_rf) + plot_cm(results_xgb) +
  plot_annotation(title = "Confusion Matrices — Internal Validation",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(OUT_DIR, "confusion_matrices.png"), p_cm_all, width = 15, height = 5, dpi = 300)

# Quantile bins that tolerate tied predictions. With six mostly binary
# predictors the calibrated probabilities take few distinct values, so the
# decile cut-points repeat and cut() fails ("'breaks' are not unique").
# Duplicate cut-points are dropped, so some bins hold more patients.
quantile_bins <- function(p, n_bins = 10) {
  brks <- unique(quantile(p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(brks) < 2) return(factor(rep("all", length(p))))
  cut(p, breaks = brks, include.lowest = TRUE)
}

plot_calibration <- function(probs, labels_bin, brier, model_name, n_bins = 10) {
  brier_null  <- mean(labels_bin) * (1 - mean(labels_bin))
  brier_skill <- round(1 - (brier / brier_null), 3)
  calib_df <- tibble(prob = probs, label = labels_bin) %>%
    mutate(bin = quantile_bins(prob, n_bins)) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(prob), mean_obs = mean(label), n = n(),
              n_events = sum(label), .groups = "drop") %>%
    filter(!is.na(bin))
  if (any(sdc_small(calib_df$n_events)))
    warning(sprintf("Calibration plot %s: %d bin(s) with 1-9 cancers - observed rate x bin N ",
                    model_name, sum(sdc_small(calib_df$n_events))),
            "would reveal the count. Internal figure only, or check with the custodian.")
  ggplot(calib_df, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(linetype = "dashed", colour = "grey50", linewidth = 1) +
    geom_line(colour = "#2166AC", linewidth = 1) +
    geom_point(aes(size = n), colour = "#2166AC", alpha = 0.8) +
    scale_size_continuous(name = "N in bin", range = c(2, 8)) +
    annotate("text", x = max(calib_df$mean_pred) * 0.05, y = 0.9,
             label = paste0("Brier = ", round(brier, 4),
                            "\nSkill = ", brier_skill),
             hjust = 0, size = 3.5) +
    labs(title = paste("Calibration —", model_name),
         x = "Mean Predicted Probability (quantile bins)",
         y = "Observed Event Rate") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
}

p_cal_all <- plot_calibration(lasso_prob_test_val_cal, y_test_val_bin,
                              brier_cal_lasso, "LASSO (Platt-calibrated)") +
  plot_calibration(rf_prob_test_val_cal,    y_test_val_bin,
                   brier_cal_rf,    "Random Forest (Platt-calibrated)") +
  plot_calibration(xgb_prob_test_val_cal,   y_test_val_bin,
                   brier_cal_xgb,   "XGBoost (Platt-calibrated)") +
  plot_annotation(title = "Calibration Curves — Internal Validation (Platt-scaled)",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))
ggsave(file.path(OUT_DIR, "calibration_curves.png"), p_cal_all, width = 15, height = 5, dpi = 300)

rf_imp_df <- importance(rf_val) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  arrange(desc(MeanDecreaseGini))

p_rf_imp <- ggplot(rf_imp_df, aes(x = reorder(Feature, MeanDecreaseGini),
                                  y = MeanDecreaseGini)) +
  geom_col(fill = "#4DAF4A", alpha = 0.85) + coord_flip() +
  labs(title = "Random Forest — Feature Importance (Mean Decrease Gini)",
       x = NULL, y = "Mean Decrease in Gini") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "feature_importance_rf.png"), p_rf_imp, width = 9, height = 6, dpi = 300)

xgb_imp    <- xgb.importance(model = xgb_val, feature_names = colnames(X_train_val_mat))
xgb_imp_df <- as.data.frame(xgb_imp) %>% arrange(desc(Gain))

p_xgb_imp <- ggplot(xgb_imp_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_col(fill = "#E41A1C", alpha = 0.85) + coord_flip() +
  labs(title = "XGBoost — Feature Importance (Gain)", x = NULL, y = "Gain") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "feature_importance_xgb.png"), p_xgb_imp, width = 9, height = 6, dpi = 300)

p_lasso_imp <- ggplot(lasso_coef,
                      aes(x = reorder(Variable, abs(Coefficient)),
                          y = Coefficient, fill = Coefficient > 0)) +
  geom_col(alpha = 0.85) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#D73027"),
                    labels = c("Negative", "Positive"), name = "Direction") +
  labs(title = "LASSO — Non-Zero Coefficients", x = NULL, y = "Coefficient") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "feature_importance_lasso.png"), p_lasso_imp, width = 9, height = 6, dpi = 300)

cat("=== Computing SHAP Values (XGBoost) ===\n\n")
shap       <- shapviz(xgb_val, X_pred = X_test_val_mat)
p_shap_sum <- sv_importance(shap, kind = "beeswarm") +
  labs(title = "SHAP Summary Plot — XGBoost") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "shap_summary_xgb.png"), p_shap_sum, width = 9, height = 7, dpi = 300)

p_shap_bar <- sv_importance(shap, kind = "bar") +
  labs(title = "SHAP Feature Importance — XGBoost", x = "Mean |SHAP Value|") +
  theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUT_DIR, "shap_importance_bar_xgb.png"), p_shap_bar, width = 9, height = 6, dpi = 300)


# =============================================================================
# 13. SENSITIVITY ANALYSIS — SMOTE
# =============================================================================

cat("=== SENSITIVITY ANALYSIS: SMOTE ===\n\n")
set.seed(42)
smote_train <- ROSE::ovun.sample(
  outcome ~ ., data = dat_train_val, method = "over", seed = 42
)$data

cat("SMOTE training set class distribution:\n")
print(table(smote_train$outcome))

X_smote_mat <- model.matrix(~ . - outcome - 1, data = smote_train)
y_smote_bin <- ifelse(smote_train$outcome == "Cancer", 1, 0)
dsmote      <- xgb.DMatrix(data = X_smote_mat, label = y_smote_bin)

set.seed(42)
xgb_smote <- xgb.train(
  params  = list(objective = "binary:logistic", eval_metric = "auc",
                 eta = best_params$eta, max_depth = best_params$max_depth,
                 nthread = parallel::detectCores() - 1),
  data    = dsmote,
  nrounds = best_params$nrounds,
  verbose = 0
)

xgb_smote_prob_val <- predict(xgb_smote, dtest_val)
results_xgb_smote  <- evaluate_model(xgb_smote_prob_val, y_test_val_bin,
                                     thresh_xgb$threshold,
                                     "XGBoost (SMOTE)", "Validation (SMOTE Sensitivity)")


# =============================================================================
# 14. RE-FIT BEST MODEL ON FULL SANT DATASET
# =============================================================================

# Override model selection — use LASSO for full dataset re-fit and transfer
# (LASSO had the highest internal-validation AUROC; XGBoost tuning settled on
# very few trees). The XGBoost and RF branches below are kept but not run.
best_model_name <- "LASSO"
cat(sprintf("=> Model override: using %s for full dataset re-fit and transfer bundle.\n\n",
            best_model_name))

cat("\n=== RE-FITTING BEST MODEL ON FULL SANT DATASET ===\n\n")
cat("Best model:", best_model_name, "\n")
cat("Full SANT N:", nrow(dat_analysis), "\n\n")

X_full     <- dat_analysis %>% select(-outcome)
y_full     <- dat_analysis$outcome
y_full_bin <- ifelse(y_full == "Cancer", 1, 0)

n_pos_full    <- sum(y_full == "Cancer")
n_neg_full    <- sum(y_full == "No_Cancer")
scale_pw_full <- n_neg_full / n_pos_full

X_full_mat <- model.matrix(~ . - 1, data = X_full)
dfull      <- xgb.DMatrix(data = X_full_mat, label = y_full_bin)

calibrate_platt_return <- function(prob_train, y_train_bin, model_name) {
  cal_df  <- data.frame(prob = prob_train, y = y_train_bin)
  platt_m <- glm(y ~ prob, data = cal_df, family = binomial)
  cat(model_name, "- Platt model fitted on full SANT set.\n")
  platt_m
}

if (best_model_name == "XGBoost") {
  
  set.seed(42)
  model_stage1 <- xgb.train(
    params = list(
      objective        = "binary:logistic",
      eval_metric      = "aucpr",
      eta              = best_params$eta,
      max_depth        = best_params$max_depth,
      gamma            = best_params$gamma,
      colsample_bytree = best_params$colsample_bytree,
      min_child_weight = best_params$min_child_weight,
      subsample        = best_params$subsample,
      scale_pos_weight = best_params$scale_pos_weight,
      nthread          = parallel::detectCores() - 1
    ),
    data    = dfull,
    nrounds = best_params$nrounds,
    verbose = 0
  )
  
  cat("\n--- Stage 2: max_depth refinement on full SANT cohort ---\n\n")
  
  max_depth_grid <- c(2, 3, 4, 5, 6, 7, 8, 9, 10)
  
  depth_cv_results <- purrr::map_dfr(max_depth_grid, function(depth) {
    set.seed(42)
    cv <- xgb.cv(
      params = list(
        objective        = "binary:logistic",
        eval_metric      = "aucpr",
        eta              = best_params$eta,
        max_depth        = depth,
        gamma            = best_params$gamma,
        colsample_bytree = best_params$colsample_bytree,
        min_child_weight = best_params$min_child_weight,
        subsample        = best_params$subsample,
        scale_pos_weight = best_params$scale_pos_weight,
        nthread          = parallel::detectCores() - 1
      ),
      data                  = dfull,
      nrounds               = best_params$nrounds * 2,
      nfold                 = 5,
      stratified            = TRUE,
      verbose               = FALSE,
      showsd                = TRUE,
      early_stopping_rounds = 20
    )
    
    best_round <- if (!is.null(cv$best_iteration) && !is.na(cv$best_iteration)) {
      cv$best_iteration
    } else {
      which.max(cv$evaluation_log$test_aucpr_mean)
    }
    
    tibble(
      max_depth     = depth,
      best_nrounds  = best_round,
      cv_auprc_mean = max(cv$evaluation_log$test_aucpr_mean, na.rm = TRUE),
      cv_auprc_std  = cv$evaluation_log$test_aucpr_std[best_round]
    )
  })
  
  cat("max_depth CV results (full SANT, AUPRC):\n")
  print(depth_cv_results)
  
  best_depth_row <- depth_cv_results %>%
    slice_max(cv_auprc_mean, n = 1, with_ties = FALSE)
  
  cat(sprintf(
    "\n=> Best max_depth: %d  |  CV AUPRC: %.4f (±%.4f)  |  nrounds: %d\n\n",
    best_depth_row$max_depth,
    best_depth_row$cv_auprc_mean,
    best_depth_row$cv_auprc_std,
    best_depth_row$best_nrounds
  ))
  
  stage1_auprc <- max(predict(model_stage1, dfull) %>%
                        { pr.curve(scores.class0 = .[y_full_bin == 1],
                                   scores.class1 = .[y_full_bin == 0])$auc.integral })
  
  cat(sprintf("Stage 1 AUPRC (best_params max_depth = %d): %.4f\n",
              best_params$max_depth, stage1_auprc))
  cat(sprintf("Stage 2 CV AUPRC (refined max_depth = %d):  %.4f\n",
              best_depth_row$max_depth, best_depth_row$cv_auprc_mean))
  
  final_max_depth <- best_depth_row$max_depth
  final_nrounds   <- best_depth_row$best_nrounds
  
  set.seed(42)
  model_final <- xgb.train(
    params = list(
      objective        = "binary:logistic",
      eval_metric      = "aucpr",
      eta              = best_params$eta,
      max_depth        = final_max_depth,
      gamma            = best_params$gamma,
      colsample_bytree = best_params$colsample_bytree,
      min_child_weight = best_params$min_child_weight,
      subsample        = best_params$subsample,
      scale_pos_weight = best_params$scale_pos_weight,
      nthread          = parallel::detectCores() - 1
    ),
    data    = dfull,
    nrounds = final_nrounds,
    verbose = 0
  )
  
  prob_full_train    <- predict(model_final, dfull)
  platt_model_final  <- calibrate_platt_return(prob_full_train, y_full_bin, "XGBoost")
  feature_names_used <- colnames(X_full_mat)
  
  best_params <- best_params %>%
    mutate(max_depth = final_max_depth,
           nrounds   = final_nrounds)
  
  cat(sprintf("\nFinal model: max_depth = %d | nrounds = %d\n\n",
              final_max_depth, final_nrounds))
  
  p_depth <- ggplot(depth_cv_results,
                    aes(x = factor(max_depth), y = cv_auprc_mean)) +
    geom_col(fill = "#E41A1C", alpha = 0.85) +
    geom_errorbar(aes(ymin = cv_auprc_mean - cv_auprc_std,
                      ymax = cv_auprc_mean + cv_auprc_std),
                  width = 0.2) +
    geom_point(size = 3, colour = "black") +
    labs(title = "Stage 2: max_depth Refinement on Full SANT Cohort (CV AUPRC)",
         x = "max_depth", y = "CV AUPRC (mean ± SD)") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUT_DIR, "depth_refinement_full_cohort.png"), p_depth, width = 8, height = 6, dpi = 300)
  
} else if (best_model_name == "Random Forest") {
  
  rf_classwt_full <- c("No_Cancer" = 1, "Cancer" = scale_pw_full)
  set.seed(42)
  model_final <- randomForest(
    x          = X_full,
    y          = y_full,
    ntree      = 500,
    mtry       = rf_tune$bestTune$mtry,
    classwt    = rf_classwt_full,
    importance = TRUE
  )
  prob_full_train    <- predict(model_final, X_full, type = "prob")[, "Cancer"]
  platt_model_final  <- calibrate_platt_return(prob_full_train, y_full_bin, "Random Forest")
  feature_names_used <- names(X_full)
  
} else {  # LASSO
  
  case_wts_full_lasso <- ifelse(y_full == "Cancer", scale_pw_full, 1)
  set.seed(42)
  model_final <- glmnet(
    x           = X_full_mat,
    y           = y_full_bin,
    family      = "binomial",
    alpha       = 1,
    lambda      = lambda_best,
    weights     = case_wts_full_lasso,
    standardize = TRUE
  )
  prob_full_train    <- predict(model_final, X_full_mat, type = "response")[, 1]
  platt_model_final  <- calibrate_platt_return(prob_full_train, y_full_bin, "LASSO")
  feature_names_used <- colnames(X_full_mat)
}

cat("\nFull-dataset re-fit complete.\n\n")


# =============================================================================
# 14B. 5-FOLD CROSS-VALIDATION ON FULL SANT DATASET
# =============================================================================

cat("\n=== STEP 14B: 5-FOLD CV ON FULL SANT DATASET ===\n\n")
cat("Model type:", best_model_name, "\n")
cat("Full SANT N:", nrow(dat_analysis), "\n\n")

set.seed(42)
cv_fold_idx <- createFolds(dat_analysis$outcome, k = 5, list = TRUE, returnTrain = FALSE)

oof_probs  <- numeric(nrow(dat_analysis))
oof_labels <- numeric(nrow(dat_analysis))
fold_models        <- vector("list", 5)
fold_platt_models  <- vector("list", 5)
fold_auroc         <- numeric(5)
fold_auprc         <- numeric(5)
fold_sensitivity   <- numeric(5)
fold_specificity   <- numeric(5)
fold_f1            <- numeric(5)
fold_brier         <- numeric(5)

for (fold in 1:5) {
  
  cat(sprintf("--- Fold %d / 5 ---\n", fold))
  
  test_idx     <- cv_fold_idx[[fold]]
  train_idx_cv <- setdiff(seq_len(nrow(dat_analysis)), test_idx)
  
  fold_train <- dat_analysis[ train_idx_cv, ]
  fold_test  <- dat_analysis[ test_idx,     ]
  
  X_fold_train <- fold_train %>% select(-outcome)
  y_fold_train <- fold_train$outcome
  X_fold_test  <- fold_test  %>% select(-outcome)
  y_fold_test  <- fold_test$outcome
  
  y_fold_train_bin <- ifelse(y_fold_train == "Cancer", 1, 0)
  y_fold_test_bin  <- ifelse(y_fold_test  == "Cancer", 1, 0)
  
  n_pos_fold <- sum(y_fold_train_bin)
  n_neg_fold <- sum(1 - y_fold_train_bin)
  spw_fold   <- n_neg_fold / max(n_pos_fold, 1)
  
  if (best_model_name == "XGBoost") {
    
    X_fold_train_mat <- model.matrix(~ . - 1, data = X_fold_train)
    X_fold_test_mat  <- model.matrix(~ . - 1, data = X_fold_test)
    dtrain_fold      <- xgb.DMatrix(data = X_fold_train_mat, label = y_fold_train_bin)
    dtest_fold       <- xgb.DMatrix(data = X_fold_test_mat,  label = y_fold_test_bin)
    
    set.seed(42)
    fold_model <- xgb.train(
      params = list(
        objective        = "binary:logistic",
        eval_metric      = "aucpr",
        eta              = best_params$eta,
        max_depth        = best_params$max_depth,
        gamma            = best_params$gamma,
        colsample_bytree = best_params$colsample_bytree,
        min_child_weight = best_params$min_child_weight,
        subsample        = best_params$subsample,
        scale_pos_weight = spw_fold,
        nthread          = parallel::detectCores() - 1
      ),
      data    = dtrain_fold,
      nrounds = best_params$nrounds,
      verbose = 0
    )
    
    prob_fold_train <- predict(fold_model, dtrain_fold)
    prob_fold_test  <- predict(fold_model, dtest_fold)
    
  } else if (best_model_name == "Random Forest") {
    
    rf_wt_fold <- c("No_Cancer" = 1, "Cancer" = spw_fold)
    set.seed(42)
    fold_model <- randomForest(
      x       = X_fold_train,
      y       = y_fold_train,
      ntree   = 500,
      mtry    = rf_tune$bestTune$mtry,
      classwt = rf_wt_fold
    )
    prob_fold_train <- predict(fold_model, X_fold_train, type = "prob")[, "Cancer"]
    prob_fold_test  <- predict(fold_model, X_fold_test,  type = "prob")[, "Cancer"]
    
  } else {  # LASSO
    
    X_fold_train_mat <- model.matrix(~ . - 1, data = X_fold_train)
    X_fold_test_mat  <- model.matrix(~ . - 1, data = X_fold_test)
    wts_fold         <- ifelse(y_fold_train == "Cancer", spw_fold, 1)
    
    set.seed(42)
    fold_model <- glmnet(
      x           = X_fold_train_mat,
      y           = y_fold_train_bin,
      family      = "binomial",
      alpha       = 1,
      lambda      = lambda_best,
      weights     = wts_fold,
      standardize = TRUE
    )
    prob_fold_train <- predict(fold_model, X_fold_train_mat, type = "response")[, 1]
    prob_fold_test  <- predict(fold_model, X_fold_test_mat,  type = "response")[, 1]
  }
  
  cal_df_fold   <- data.frame(prob = prob_fold_train, y = y_fold_train_bin)
  platt_fold    <- glm(y ~ prob, data = cal_df_fold, family = binomial)
  prob_fold_cal <- predict(platt_fold,
                           newdata = data.frame(prob = prob_fold_test),
                           type    = "response")
  
  oof_probs[test_idx]  <- prob_fold_cal
  oof_labels[test_idx] <- y_fold_test_bin
  
  roc_fold <- roc(y_fold_test_bin, prob_fold_cal, quiet = TRUE)
  auroc_f  <- as.numeric(auc(roc_fold))
  pr_fold  <- pr.curve(scores.class0 = prob_fold_cal[y_fold_test_bin == 1],
                       scores.class1 = prob_fold_cal[y_fold_test_bin == 0],
                       curve = FALSE)
  auprc_f  <- pr_fold$auc.integral
  brier_f  <- mean((prob_fold_cal - y_fold_test_bin)^2)
  
  thresh_f <- optimise_threshold(prob_fold_train, y_fold_train_bin,
                                 paste0("Fold ", fold))$threshold
  pred_f   <- ifelse(prob_fold_cal >= thresh_f, 1, 0)
  tp_f <- sum(pred_f == 1 & y_fold_test_bin == 1)
  fp_f <- sum(pred_f == 1 & y_fold_test_bin == 0)
  fn_f <- sum(pred_f == 0 & y_fold_test_bin == 1)
  tn_f <- sum(pred_f == 0 & y_fold_test_bin == 0)
  sens_f <- tp_f / max(tp_f + fn_f, 1)
  spec_f <- tn_f / max(tn_f + fp_f, 1)
  prec_f <- ifelse((tp_f + fp_f) == 0, 0, tp_f / (tp_f + fp_f))
  f1_f   <- ifelse((sens_f + prec_f) == 0, 0,
                   2 * sens_f * prec_f / (sens_f + prec_f))
  
  fold_models[[fold]]       <- fold_model
  fold_platt_models[[fold]] <- platt_fold
  fold_auroc[[fold]]        <- auroc_f
  fold_auprc[[fold]]        <- auprc_f
  fold_sensitivity[[fold]]  <- sens_f
  fold_specificity[[fold]]  <- spec_f
  fold_f1[[fold]]           <- f1_f
  fold_brier[[fold]]        <- brier_f
  
  cat(sprintf("  AUROC: %.4f | AUPRC: %.4f | Sensitivity: %.4f | Specificity: %.4f | F1: %.4f | Brier: %.4f\n",
              auroc_f, auprc_f, sens_f, spec_f, f1_f, brier_f))
}

oof_roc   <- roc(oof_labels, oof_probs, quiet = TRUE)
oof_auroc <- as.numeric(auc(oof_roc))
oof_pr    <- pr.curve(scores.class0 = oof_probs[oof_labels == 1],
                      scores.class1 = oof_probs[oof_labels == 0],
                      curve = FALSE)
oof_auprc <- oof_pr$auc.integral
oof_brier <- mean((oof_probs - oof_labels)^2)

cv_full_summary <- data.frame(
  Fold        = c(1:5, "OOF aggregate"),
  AUROC       = round(c(fold_auroc,       oof_auroc),  4),
  AUPRC       = round(c(fold_auprc,       oof_auprc),  4),
  Sensitivity = round(c(fold_sensitivity, NA),         4),
  Specificity = round(c(fold_specificity, NA),         4),
  F1          = round(c(fold_f1,          NA),         4),
  Brier       = round(c(fold_brier,       oof_brier),  4)
)

cat("\n=== 5-FOLD CV RESULTS (Full SANT Dataset) ===\n")
print(cv_full_summary)

cat(sprintf("\nMean CV AUROC:       %.4f (SD = %.4f)\n", mean(fold_auroc), sd(fold_auroc)))
cat(sprintf("Mean CV AUPRC:       %.4f (SD = %.4f)\n", mean(fold_auprc), sd(fold_auprc)))
cat(sprintf("Mean CV Sensitivity: %.4f (SD = %.4f)\n", mean(fold_sensitivity), sd(fold_sensitivity)))
cat(sprintf("Mean CV Specificity: %.4f (SD = %.4f)\n", mean(fold_specificity), sd(fold_specificity)))
cat(sprintf("Mean CV F1:          %.4f (SD = %.4f)\n", mean(fold_f1), sd(fold_f1)))
cat(sprintf("OOF aggregate AUROC: %.4f\n", oof_auroc))
cat(sprintf("OOF aggregate AUPRC: %.4f\n\n", oof_auprc))

write.csv(cv_full_summary, file.path(OUT_DIR, "cv_full_sant_performance.csv"), row.names = FALSE)

best_fold_idx <- which.max(fold_auroc)
model_cv_best <- fold_models[[best_fold_idx]]
platt_cv_best <- fold_platt_models[[best_fold_idx]]

cat(sprintf("Best fold: Fold %d (AUROC = %.4f) — saved to transfer bundle.\n\n",
            best_fold_idx, fold_auroc[best_fold_idx]))

cv_plot_df <- data.frame(
  Fold  = factor(1:5),
  AUROC = fold_auroc,
  AUPRC = fold_auprc
) %>%
  pivot_longer(cols = c(AUROC, AUPRC), names_to = "Metric", values_to = "Value")

p_cv <- ggplot(cv_plot_df, aes(x = Fold, y = Value, fill = Metric)) +
  geom_col(position = "dodge", alpha = 0.85) +
  geom_hline(yintercept = oof_auroc, linetype = "dashed",
             colour = "#2166AC", linewidth = 0.8) +
  geom_hline(yintercept = oof_auprc, linetype = "dotted",
             colour = "#E41A1C", linewidth = 0.8) +
  scale_fill_manual(values = c("AUROC" = "#2166AC", "AUPRC" = "#E41A1C")) +
  labs(title   = paste0("5-Fold CV Performance — Full SANT Dataset (", best_model_name, ")"),
       x       = "Fold",
       y       = "Metric Value",
       caption = "Dashed = OOF AUROC. Dotted = OOF AUPRC.") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(OUT_DIR, "cv_full_sant_performance.png"), p_cv, width = 9, height = 6, dpi = 300)


# =============================================================================
# 14C. SANT 1-YEAR OUTPUTS: TABLE S6 (SANT), FIGURE S2 (SHAP), FIGURE S3(A)
# =============================================================================

cat("\n=== STEP 14C: SANT 1-YEAR OUTPUTS ===\n\n")

# Platt-calibrated predictions for the full SANT cohort (final model)
prob_full_cal <- as.numeric(predict(platt_model_final,
                                    newdata = data.frame(prob = prob_full_train),
                                    type    = "response"))

roc_full_cal <- roc(y_full_bin, prob_full_cal, quiet = TRUE)
sant_youden  <- as.numeric(coords(roc_full_cal, "best", best.method = "youden",
                                  ret = "threshold", transpose = FALSE)[1, 1])
cat(sprintf("SANT Youden threshold (calibrated 1-year probability): %.4f\n", sant_youden))
cat("\n")

# ---- metrics at a threshold (Table S6 layout) ------------------------------
s6_metrics <- function(probs, y, thr, platt_label) {
  flag <- probs >= thr
  tp <- sum(flag & y == 1); fp <- sum(flag & y == 0)
  fn <- sum(!flag & y == 1); tn <- sum(!flag & y == 0)
  roc_o <- roc(y, probs, quiet = TRUE)
  ci    <- as.numeric(ci.auc(roc_o, method = "delong"))
  auprc <- pr.curve(scores.class0 = probs[y == 1],
                    scores.class1 = probs[y == 0])$auc.integral
  brier <- mean((probs - y)^2)
  bss   <- 1 - brier / (mean(y) * (1 - mean(y)))
  ppv   <- tp / max(tp + fp, 1)
  # small-cell rule: any metric built from a 1-9 cell is shown as "S"
  s <- function(cells, txt) if (any(sdc_small(cells))) "S" else txt
  c(`Platt calibration model`  = platt_label,
    `Threshold applied`        = sprintf("%.4f", thr),
    `AUROC (95% CI)`           = sprintf("%.3f (%.3f–%.3f)", ci[2], ci[1], ci[3]),
    `AUPRC`                    = sprintf("%.4f", auprc),
    `Sensitivity`              = s(c(tp, fn), sprintf("%.1f%%", 100 * tp / (tp + fn))),
    `Specificity`              = s(c(tn, fp), sprintf("%.1f%%", 100 * tn / (tn + fp))),
    `NNI`                      = s(c(tp, fp, fn), sprintf("%.1f", 1 / ppv)),
    `PPV`                      = s(c(tp, fp, fn), sprintf("%.2f%%", 100 * ppv)),
    `NPV`                      = s(c(tn, fn, tp), sprintf("%.2f%%", 100 * tn / max(tn + fn, 1))),
    `Cancers caught / total`   = sprintf("%s / %d",
                                         if (sdc_small(tp)) "<10" else if (sdc_small(fn)) "S" else tp,
                                         tp + fn),
    `Patients flagged, n (%)`  = if (sdc_small(tp + fp)) "<10" else if (sdc_small(fp)) "S" else
      sprintf("%s (%.1f%%)", format(tp + fp, big.mark = ","), 100 * (tp + fp) / length(y)),
    `Brier score (calibrated)` = sprintf("%.4f", brier),
    `Brier skill score`        = sprintf("%.4f", bss))
}

s6_cols <- list(
  `SANT Derivation (Youden)` =
    s6_metrics(prob_full_cal, y_full_bin, sant_youden,     "SANT-derived"),
  `SANT 5-fold CV (out-of-fold)` =
    s6_metrics(oof_probs, oof_labels, sant_youden, "Fold-specific Platt")
)

# ---- Cohort reconciliation with the Cox pipeline (dat_cox) -----------------
cat("=== COHORT CHECK vs Cox pipeline ===\n")
cat(sprintf("XGBoost analytic N: %d | 1-year events: %d | (before complete-case: N = %d)\n",
            nrow(dat_analysis), sum(y_full_bin), nrow(dat_ml)))
if (exists("dat_cox")) {
  same <- nrow(dat_analysis) == nrow(dat_cox) && sum(y_full_bin) == sum(dat_cox$event)
  cat(sprintf("Cox dat_cox N: %d | 1-year events: %d  ->  identical cohort: %s\n\n",
              nrow(dat_cox), sum(dat_cox$event), same))
  if (!same) warning("XGBoost and Cox cohorts differ - check missing-data handling.")
} else cat("(dat_cox not in session - run the Cox pipeline first to compare)\n\n")
if (imputed_flag) warning("MICE imputation was used: XGBoost keeps rows the Cox model drops.")

table_s6_sant <- bind_rows(
  tibble(Metric = c("N (AP)", "Cancer events within 1 year, n (%)",
                    "Cancer events, full follow-up, n (%)"),
         Value  = c(format(nrow(dat_analysis), big.mark = ","),
                    sprintf("%d (%.2f%%)", sum(y_full_bin), 100 * mean(y_full_bin)),
                    sprintf("%d (%.2f%%)", sum(dat_ml$event_full), 100 * mean(dat_ml$event_full)))) %>%
    mutate(`SANT Derivation (Youden)` = Value,
           `SANT 5-fold CV (out-of-fold)` = Value) %>%
    select(-Value),
  tibble(Metric = names(s6_cols[[1]]),
         !!!lapply(s6_cols, unname))
)

cat("=== TABLE S6 (SANT columns, 1-year outcome) ===\n")
print(as.data.frame(table_s6_sant), row.names = FALSE)
write.csv(table_s6_sant, file.path(OUT_DIR, "Table_S6_SANT_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table_S6_SANT_1yr.csv\n\n")

# ---- TABLE S7: model selection comparison (SANT 20% internal holdout) ------
s7_col <- function(res, brier_cal, bss) {
  ct <- res$cm$table
  tp <- ct["Cancer", "Cancer"];    fn <- ct["No_Cancer", "Cancer"]
  fp <- ct["Cancer", "No_Cancer"]; tn <- ct["No_Cancer", "No_Cancer"]
  s  <- function(cells, txt) if (any(sdc_small(cells))) "S" else txt
  c(AUROC                 = sprintf("%.3f", res$auroc),
    AUPRC                 = sprintf("%.3f", res$auprc),
    Sensitivity           = s(c(tp, fn), sprintf("%.1f%%", 100 * res$sensitivity)),
    Specificity           = s(c(tn, fp), sprintf("%.1f%%", 100 * res$specificity)),
    PPV                   = s(c(tp, fp, fn), sprintf("%.1f%%", 100 * res$ppv)),
    `F1 score`            = s(c(tp, fp, fn), sprintf("%.3f", res$f1)),
    `Brier (calibrated)`  = sprintf("%.4f", brier_cal),
    `Brier skill score`   = sprintf("%.3f", bss))
}
s7_list <- list(
  LASSO           = s7_col(results_lasso, brier_cal_lasso, brier_skill_lasso),
  `Random Forest` = s7_col(results_rf,    brier_cal_rf,    brier_skill_rf),
  XGBoost         = s7_col(results_xgb,   brier_cal_xgb,   brier_skill_xgb)
)
table_s7 <- tibble(Metric = names(s7_list[[1]]), !!!lapply(s7_list, unname))

cat(sprintf("=== TABLE S7: MODEL SELECTION - SANT 20%% HOLDOUT (N = %d, %d 1-year cancers) ===\n",
            length(y_test_val_bin), sum(y_test_val_bin)))
print(as.data.frame(table_s7), row.names = FALSE)
write.csv(table_s7, file.path(OUT_DIR, "Table_S7_model_selection_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table_S7_model_selection_1yr.csv\n\n")

# ---- TABLE S8: final model (LASSO) 5-fold CV (full SANT cohort) -------------
final_auprc <- pr.curve(scores.class0 = prob_full_cal[y_full_bin == 1],
                        scores.class1 = prob_full_cal[y_full_bin == 0])$auc.integral
final_brier <- mean((prob_full_cal - y_full_bin)^2)

table_s8 <- tibble(
  Fold  = c(as.character(1:5), "OOF Aggregate", "Full cohort (final model)"),
  AUROC = sprintf("%.3f", c(fold_auroc, oof_auroc, as.numeric(auc(roc_full_cal)))),
  AUPRC = sprintf("%.3f", c(fold_auprc, oof_auprc, final_auprc)),
  Brier = sprintf("%.4f", c(fold_brier, oof_brier, final_brier))
)
table_s8_footnote <- sprintf(paste0(
  "Platt calibration applied within each fold. OOF aggregate AUROC = %.3f; ",
  "OOF AUPRC = %.3f. Outcome: pancreatic cancer within 1 year of AP ",
  "(%d events, N = %s)."),
  oof_auroc, oof_auprc, sum(y_full_bin), format(length(y_full_bin), big.mark = ","))

cat(sprintf("=== TABLE S8: %s 5-FOLD CV (FULL SANT COHORT) ===\n", toupper(best_model_name)))
print(as.data.frame(table_s8), row.names = FALSE)
cat("Footnote:", table_s8_footnote, "\n")
write.csv(table_s8, file.path(OUT_DIR, sprintf("Table_S8_%s_5fold_CV_1yr.csv", best_model_name)),
          row.names = FALSE, fileEncoding = "UTF-8")
writeLines(table_s8_footnote, file.path(OUT_DIR, "Table_S8_footnote.txt"))
cat(sprintf("Saved: Table_S8_%s_5fold_CV_1yr.csv / Table_S8_footnote.txt\n\n", best_model_name))

# ---- AUROC summary (all models + CV) ---------------------------------------
auc_ci_str <- function(y, p) {
  ci <- as.numeric(ci.auc(roc(y, p, quiet = TRUE), method = "delong"))
  sprintf("%.3f (%.3f–%.3f)", ci[2], ci[1], ci[3])
}
auroc_summary <- tibble(
  Analysis = c("LASSO - internal validation (20%)",
               "Random Forest - internal validation (20%)",
               "XGBoost - internal validation (20%)",
               paste(best_model_name, "- 5-fold CV, out-of-fold (full SANT)"),
               paste(best_model_name, "- final model, apparent (full SANT)")),
  `AUROC (95% CI)` = c(auc_ci_str(y_test_val_bin, lasso_prob_test_val),
                       auc_ci_str(y_test_val_bin, rf_prob_test_val),
                       auc_ci_str(y_test_val_bin, xgb_prob_test_val),
                       auc_ci_str(oof_labels, oof_probs),
                       auc_ci_str(y_full_bin, prob_full_cal)),
  Events = c(rep(sum(y_test_val_bin), 3), sum(oof_labels), sum(y_full_bin))
)
cat("=== AUROC SUMMARY (1-year outcome) ===\n")
print(as.data.frame(auroc_summary), row.names = FALSE)
write.csv(auroc_summary, file.path(OUT_DIR, "AUROC_summary_1yr.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: AUROC_summary_1yr.csv\n\n")

# ---- Figure S2: SHAP beeswarm, development cohort (final model) -----------
feat_labels <- c(nv_k85_firstage = "First age of AP", male = "Male",
                 biliary_ap = "Biliary AP", alcohol_ap = "Alcohol AP",
                 unspecified_ap = "Unspecified AP", other_ap = "Other AP")

if (best_model_name == "LASSO") {
  # Exact SHAP for a linear model (log-odds scale): beta_j * (x_j - mean(x_j)).
  # glmnet returns coefficients on the original scale even with standardize = TRUE.
  beta   <- as.matrix(coef(model_final))[colnames(X_full_mat), 1]
  X_cent <- sweep(X_full_mat, 2, colMeans(X_full_mat))
  S_lin  <- sweep(X_cent, 2, beta, `*`)
  base_lp <- as.numeric(coef(model_final)[1, 1] + sum(colMeans(X_full_mat) * beta))
  shap_full <- shapviz(S_lin, X = as.data.frame(X_full_mat), baseline = base_lp)
} else {
  shap_full <- shapviz(model_final, X_pred = X_full_mat)
}
lab <- unname(feat_labels[colnames(shap_full$S)])
colnames(shap_full$S) <- lab
colnames(shap_full$X) <- lab

p_s2 <- sv_importance(shap_full, kind = "beeswarm",
                      viridis_args = list(option = "inferno", begin = 0.2, end = 0.85)) +
  labs(x = "SHAP value") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "FigureS2_SHAP_development_1yr.png"),
       p_s2, width = 9, height = 6.5, dpi = 400)
cat("Saved: FigureS2_SHAP_development_1yr.png\n")

# ---- Figure S3(A): SANT calibration (Platt-calibrated, quantile bins) ------
calib_panel <- function(probs, y, colour, title, n_bins = 10) {
  df <- tibble(prob = probs, y = y) %>%
    mutate(bin = quantile_bins(prob, n_bins)) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(prob), mean_obs = mean(y), n = n(),
              n_events = sum(y), .groups = "drop")
  if (any(sdc_small(df$n_events)))
    warning(sprintf("Calibration plot %s: %d bin(s) with 1-9 cancers - observed rate x bin N ",
                    title, sum(sdc_small(df$n_events))),
            "would reveal the count. Use fewer bins or check with the custodian.")
  ggplot(df, aes(mean_pred, mean_obs)) +
    geom_abline(linetype = "dashed", colour = "grey45", linewidth = 0.9) +
    geom_line(colour = colour, linewidth = 1) +
    geom_point(aes(size = n), colour = colour) +
    scale_size_continuous(name = "N in bin", range = c(2, 6)) +
    labs(title = title, x = "Mean Predicted Probability (quantile bins)",
         y = "Observed Event Rate") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
}

p_s3a <- calib_panel(prob_full_cal, y_full_bin, "grey50", "(A)")
ggsave(file.path(OUT_DIR, "FigureS3A_calibration_SANT_1yr.png"),
       p_s3a, width = 7, height = 6, dpi = 400)
cat("Saved: FigureS3A_calibration_SANT_1yr.png\n")

# Out-of-fold version (honest internal calibration) - optional alternative
p_s3a_oof <- calib_panel(oof_probs, oof_labels, "grey50", "(A) SANT, out-of-fold")
ggsave(file.path(OUT_DIR, "FigureS3A_calibration_SANT_OOF_1yr.png"),
       p_s3a_oof, width = 7, height = 6, dpi = 400)
cat("Saved: FigureS3A_calibration_SANT_OOF_1yr.png\n\n")


# =============================================================================
# 15. SAVE TRANSFER BUNDLE
# =============================================================================

cat("=== SAVING TRANSFER BUNDLE ===\n\n")

best_threshold <- switch(best_model_name,
                         "LASSO"         = thresh_lasso$threshold,
                         "Random Forest" = thresh_rf$threshold,
                         "XGBoost"       = thresh_xgb$threshold)

# glm() objects keep a full copy of the data they were fitted on (every
# patient's outcome and predicted probability). Strip them so no
# patient-level SANT data leave in the bundle; predict(type = "response")
# still works on the stripped object.
strip_glm <- function(m) {
  if (is.null(m) || !inherits(m, "glm")) return(m)
  for (f in c("data", "model", "y", "residuals", "fitted.values", "effects",
              "linear.predictors", "weights", "prior.weights", "na.action"))
    m[[f]] <- NULL
  m$qr$qr <- NULL
  # terms/formula carry the fitting function's environment (which holds the
  # training vectors) - point them at the global env so it is not saved
  attr(m$terms, ".Environment") <- globalenv()
  if (!is.null(m$formula)) environment(m$formula) <- globalenv()
  m
}
platt_model_final_tx <- strip_glm(platt_model_final)
platt_cv_best_tx     <- strip_glm(platt_cv_best)

transfer_bundle <- list(
  best_model_name       = best_model_name,
  best_threshold        = sant_youden,      # Youden on Platt-calibrated 1-year probs
  # (= Table S6 "SANT Derivation" threshold)
  threshold_scale       = "Platt-calibrated probability from platt_model_final",
  f1_threshold_internal = best_threshold,   # old F1 threshold (raw XGB scale, 80/20 split) - not for validation
  model_final           = model_final,
  platt_model_final     = platt_model_final_tx,
  model_cv_best         = model_cv_best,
  platt_cv_best         = platt_cv_best_tx,
  best_fold_idx         = best_fold_idx,
  best_fold_auroc       = fold_auroc[best_fold_idx],
  cv_full_summary       = cv_full_summary,
  cv_oof_auroc          = oof_auroc,
  cv_oof_auprc          = oof_auprc,
  cv_mean_auroc         = mean(fold_auroc),
  cv_sd_auroc           = sd(fold_auroc),
  feature_names         = feature_names_used,
  xgb_best_params       = if (best_model_name == "XGBoost") best_params else NULL,
  lasso_lambda          = if (best_model_name == "LASSO") lambda_best else NULL,
  performance_table     = performance_table,
  imputed_flag          = imputed_flag,
  sant_prevalence       = mean(y_full_bin),
  scale_pos_weight_full = scale_pw_full,
  outcome_horizon_yrs   = HORIZON_YRS,
  sant_youden_threshold = sant_youden,
  sant_final_auroc      = as.numeric(auc(roc_full_cal)),   # SANT apparent AUROC (Table S6)
  census_date           = CENSUS_DATE,
  fu_cutoff             = FU_CUTOFF,                      # first AP must be <= this date
  created               = Sys.time()
)

saveRDS(transfer_bundle, file.path(OUT_DIR, "SANT_transfer_bundle_1yr.rds"))

# Check: stripped Platt model gives identical predictions, and nothing in
# the bundle has one row per SANT patient.
stopifnot(isTRUE(all.equal(
  as.numeric(predict(platt_model_final_tx, data.frame(prob = c(0.01, 0.1, 0.5)), type = "response")),
  as.numeric(predict(platt_model_final,    data.frame(prob = c(0.01, 0.1, 0.5)), type = "response")))))
n_sant  <- nrow(dat_analysis)
row_hit <- function(x) {
  if (is.data.frame(x) || is.matrix(x)) return(NROW(x) == n_sant)
  if (is.atomic(x)) return(length(x) == n_sant)
  if (is.list(x))   return(any(vapply(x, row_hit, logical(1))))
  FALSE
}
if (row_hit(transfer_bundle))
  stop("Transfer bundle still contains an object with one entry per SANT patient - do not send.")
cat(sprintf("Transfer bundle size: %.1f KB - no patient-level data found.\n",
            file.size(file.path(OUT_DIR, "SANT_transfer_bundle_1yr.rds")) / 1024))

cat("Transfer bundle saved: SANT_transfer_bundle_1yr.rds\n")
cat("Transfer this file to the Tasmanian computer along with the testing pipeline script.\n\n")


# =============================================================================
# 16. SAVE PERFORMANCE TABLE AND SUMMARY
# =============================================================================

write.csv(performance_table, file.path(OUT_DIR, "training_performance_summary.csv"), row.names = FALSE)

cat("=== TRAINING PIPELINE COMPLETE ===\n\n")
cat("Files generated:\n")
cat("  - SANT_transfer_bundle_1yr.rds      (TRANSFER TO TAS COMPUTER)\n")
cat("  - flowchart_SANT_XGB_1yr.csv\n")
cat(sprintf("  - Unsuppressed copies (do not request release): %s\n", INTERNAL_DIR))
cat("  - training_performance_summary.csv\n")
cat("  - cv_full_sant_performance.csv\n")
cat("  - cv_full_sant_performance.png\n")
cat("  - roc_curves.png\n")
cat("  - precision_recall_curves.png\n")
cat("  - confusion_matrices.png\n")
cat("  - calibration_curves.png\n")
cat("  - feature_importance_rf.png\n")
cat("  - feature_importance_xgb.png\n")
cat("  - feature_importance_lasso.png\n")
cat("  - shap_summary_xgb.png\n")
cat("  - shap_importance_bar_xgb.png\n")
cat("  - missing_data_pattern.png\n")
cat("  - Table_S6_SANT_1yr.csv\n")
cat("  - Table_S7_model_selection_1yr.csv\n")
cat(sprintf("  - Table_S8_%s_5fold_CV_1yr.csv / Table_S8_footnote.txt\n", best_model_name))
cat("  - AUROC_summary_1yr.csv\n")
cat("  - FigureS2_SHAP_development_1yr.png\n")
cat("  - FigureS3A_calibration_SANT_1yr.png (+ _OOF version)\n\n")
cat("All files in:", OUT_DIR, "\n\n")
cat("Transfer bundle contains:\n")
cat("  - Best model (full re-fit):    model_final + platt_model_final\n")
cat("  - Best CV fold model:          model_cv_best + platt_cv_best\n")
cat("  - All 5 fold models:           fold_models + fold_platt_models\n")
cat("  - CV summary:                  cv_full_summary\n")
cat(sprintf("  - OOF AUROC: %.4f | OOF AUPRC: %.4f\n", oof_auroc, oof_auprc))
cat("  - Best model:", best_model_name, "\n")
cat(sprintf("  - Validation threshold (SANT Youden, calibrated scale): %.4f\n", sant_youden))
if (best_model_name == "LASSO") {
  cat("  - Apply as: X <- as.matrix(X[, feature_names])\n")
  cat("              raw <- predict(model_final, newx = X, type = \"response\")[, 1]\n")
  cat("              p <- predict(platt_model_final, data.frame(prob = raw), type = \"response\"); flag if p >= best_threshold\n")
} else {
  cat("  - Apply as: p = predict(platt_model_final, data.frame(prob = predict(model_final, X)), type = \"response\"); flag if p >= best_threshold\n")
}
cat("  - Feature names:", paste(feature_names_used, collapse = ", "), "\n")
cat("\n=== READY FOR EXTERNAL VALIDATION ===\n")