# =================================================
# =================================================
# SANT Supplement and additional analysis
# =================================================
# =================================================

# =============================================================================
# PAnCReAS — TABLE S8: TIMING OF PANCREATIC CANCER DIAGNOSIS AFTER INDEX AP
#            BY PAnCReAS RISK BAND (SANT), 1-YEAR CANCERS
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
})

# -----------------------------------------------------------------------------
# SETTINGS — check these before running
# -----------------------------------------------------------------------------
DATA_PATH    <- NULL     # e.g. "C:/Data/cmbspine_pancreatitis.rds"; NULL = use object in session
OUT_DIR      <- path.expand("~/Desktop/PAnCReAS_TableS8")

# Index-admission discharge (separation) date column, if the dataset has one.
# Leave NULL if not: the first interval is then "<=30 days" only, and the
# script prints a reminder to edit the column header / footnote.
SEP_DATE_COL <- NULL     # e.g. "ap_sepdate_clean"

CENSUS_DATE   <- as.Date("2022-01-01")
HORIZON_YRS   <- 1
BAND_LOW_MAX  <- 1       # Low      : score <= 1
BAND_HIGH_MIN <- 9       # High     : score >= 9  (Moderate 2-8)
DAYS_1M       <- 30      # upper limit of interval 1
DAYS_3M       <- 91      # upper limit of interval 2 (3 months = 91 days)
SDC_MIN       <- 10      # custodian rule: counts 1-9 -> "<10"

# Expected 1-year cancers per band (from the main Cox pipeline / Table S8).
# Used only as a check that this pipeline reproduces the same patients.
EXPECTED_N <- c(Low = 16, Moderate = 88, High = 47)

INTERNAL_DIR <- file.path(OUT_DIR, "internal_unsuppressed_DO_NOT_RELEASE")
dir.create(INTERNAL_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
if (!is.null(DATA_PATH)) {
  cmbspine_pancreatitis <- if (grepl("\\.rds$", DATA_PATH, ignore.case = TRUE))
    readRDS(DATA_PATH) else read.csv(DATA_PATH, stringsAsFactors = FALSE)
}
if (!exists("cmbspine_pancreatitis"))
  stop("cmbspine_pancreatitis is not loaded. Load it, or set DATA_PATH.")

required_cols <- c(
  "ap_date_clean", "pc_date_clean", "death_date_clean",
  "k85_panctcancer", "k85_panctcancer_diffweeks", "nv_k86",
  "sant_deaths_flag", "nv_k85_firstage", "sant_nv_sex",
  "nv_k851", "nv_k852", "nv_k859", "nv_k85own", "nv_k85other", "nv_k858", "nv_k853",
  SEP_DATE_COL
)
missing_cols <- setdiff(required_cols, names(cmbspine_pancreatitis))
if (length(missing_cols) > 0)
  stop("Missing column(s): ", paste(missing_cols, collapse = ", "))


# -----------------------------------------------------------------------------
# 2. COHORT (same exclusions and outcome definition as the main Cox script)
#    1. nv_k86 == 1   2. missing first-AP date   3. cancer before first AP
#    4. cancer with no diagnosis date   5. first AP after CENSUS_DATE - 365
# -----------------------------------------------------------------------------
cutoff <- CENSUS_DATE - 365

cohort <- cmbspine_pancreatitis %>%
  mutate(
    ap_date    = as.Date(ap_date_clean,    format = "%Y-%m-%d"),
    pc_date    = as.Date(pc_date_clean,    format = "%Y-%m-%d"),
    death_date = as.Date(death_date_clean, format = "%Y-%m-%d"),
    cp_flag    = coalesce(nv_k86 == 1, FALSE),
    ca_flag    = coalesce(k85_panctcancer == 1, FALSE),
    prior_ca   = ca_flag & coalesce(k85_panctcancer_diffweeks < 0, FALSE),
    ca_nodate  = ca_flag & !prior_ca &
      (is.na(k85_panctcancer_diffweeks) | is.na(pc_date))
  ) %>%
  filter(!cp_flag, !is.na(ap_date), !prior_ca, !ca_nodate, ap_date <= cutoff) %>%
  mutate(
    event_full = as.integer(ca_flag & !is.na(pc_date) & k85_panctcancer_diffweeks >= 0),
    end_date   = case_when(
      event_full == 1                            ~ pc_date,
      sant_deaths_flag == 1 & !is.na(death_date) ~ death_date,
      TRUE                                       ~ CENSUS_DATE),
    surv_years_full = pmax(as.numeric(end_date - ap_date) / 365.25, 1 / 365.25),
    event      = as.integer(event_full == 1 & surv_years_full <= HORIZON_YRS),
    surv_years = pmin(surv_years_full, HORIZON_YRS),
    # predictors
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
  drop_na(event, surv_years, event_full, age_50_59, age_60_69, age_70plus, male,
          biliary_ap, alcohol_ap, unspecified_ap, other_ap)

cat(sprintf("Analytic cohort N = %d | 1-year cancers = %d\n\n",
            nrow(cohort), sum(cohort$event)))


# -----------------------------------------------------------------------------
# 3. REFIT 1-YEAR COX MODEL AND RE-DERIVE INTEGER POINTS
#    Reference unit = smallest |log-HR| > 0.10; points = round(log-HR / unit)
# -----------------------------------------------------------------------------
cox_fit <- coxph(Surv(surv_years, event) ~ age_50_59 + age_60_69 + age_70plus + male +
                   biliary_ap + alcohol_ap + unspecified_ap + other_ap, data = cohort)
b         <- coef(cox_fit)
ref_unit  <- min(abs(b[abs(b) > 0.10]))
point_map <- setNames(as.integer(round(b / ref_unit)), names(b))

cat("Point map (re-derived):\n"); print(point_map)

# Optional cross-check against the point map saved by the main pipeline
saved_pm <- file.path(path.expand("~/Desktop/PAnCReAS_16-10"), "cox_point_map.rds")
if (file.exists(saved_pm)) {
  pm0 <- readRDS(saved_pm)
  cat("Matches saved cox_point_map.rds:", isTRUE(all(pm0[names(point_map)] == point_map)), "\n")
}

X <- as.matrix(cohort[, names(point_map)])
cohort <- cohort %>%
  mutate(score = as.integer(X %*% point_map),
         band  = factor(case_when(score <= BAND_LOW_MAX  ~ "Low",
                                  score >= BAND_HIGH_MIN ~ "High",
                                  TRUE                   ~ "Moderate"),
                        levels = c("Low", "Moderate", "High")))


# -----------------------------------------------------------------------------
# 4. TIME FROM INDEX AP TO DIAGNOSIS, 1-YEAR CANCERS
# -----------------------------------------------------------------------------
cases <- cohort %>%
  filter(event == 1) %>%
  mutate(days = as.numeric(pc_date - ap_date))

if (!is.null(SEP_DATE_COL)) {
  cases$sep_date  <- as.Date(cases[[SEP_DATE_COL]], format = "%Y-%m-%d")
  cases$in_index  <- coalesce(cases$pc_date <= cases$sep_date, FALSE)
} else {
  cases$in_index  <- FALSE
}

cases <- cases %>%
  mutate(interval = factor(case_when(
    days <= DAYS_1M | in_index ~ "early",
    days <= DAYS_3M            ~ "mid",
    TRUE                       ~ "late"),
    levels = c("early", "mid", "late")))

# Checks
stopifnot(all(cases$days >= 0), all(cases$days <= 366))
wk_gap <- abs(cases$days / 7 - cases$k85_panctcancer_diffweeks)
if (any(wk_gap > 1, na.rm = TRUE))
  warning(sprintf("%d case(s): pc_date - ap_date differs from diffweeks by > 1 week. ",
                  sum(wk_gap > 1, na.rm = TRUE)),
          "Table uses the dates; check which source is correct.")
if (!is.null(SEP_DATE_COL))
  cat(sprintf("Cancers diagnosed in index admission but > %d days after AP: %d\n",
              DAYS_1M, sum(cases$in_index & cases$days > DAYS_1M)))

got_n <- table(cases$band)
cat("\n1-year cancers by band:\n"); print(got_n)
if (!all(got_n[names(EXPECTED_N)] == EXPECTED_N))
  warning("Band counts differ from EXPECTED_N (16 / 88 / 47). The cohort or point map ",
          "does not match the main pipeline - resolve before using this table.")


# -----------------------------------------------------------------------------
# 5. COUNTS, THEN SMALL-CELL SUPPRESSION
#    Primary:   counts 1-9 -> "<10"
#    Secondary: if a row or column has exactly one hidden cell, its value
#               could be back-calculated from the total, so the next-smallest
#               cell in that row/column is also hidden ("S"). Repeats until
#               no row or column has a single hidden cell.
# -----------------------------------------------------------------------------
cnt <- as.matrix(table(cases$band, cases$interval))
cnt <- rbind(cnt, Total = colSums(cnt))          # rows: Low, Moderate, High, Total
tot <- rowSums(cnt)                              # 1-year cancers, n

prim <- cnt > 0 & cnt < SDC_MIN
hide <- prim

repeat {
  changed <- FALSE
  # rows: 3 intervals sum to the shown row total
  for (i in seq_len(nrow(cnt))) {
    if (sum(hide[i, ]) == 1) {
      cand <- which(!hide[i, ]); j <- cand[which.min(cnt[i, cand])]
      hide[i, j] <- TRUE; changed <- TRUE
    }
  }
  # columns: 3 bands sum to the Total row
  for (j in seq_len(ncol(cnt))) {
    if (sum(hide[, j]) == 1) {
      cand <- which(!hide[, j]); i <- cand[which.min(cnt[cand, j])]
      hide[i, j] <- TRUE; changed <- TRUE
    }
  }
  if (!changed) break
}

fmt_cell <- function(k, n, p, h) {
  if (p) return("<10")
  if (h) return("S")
  sprintf("%d (%.1f)", k, 100 * k / n)
}

build_table <- function(suppress) {
  out <- data.frame(
    `Risk band` = c(sprintf("Low (<=%d)", BAND_LOW_MAX),
                    sprintf("Moderate (%d-%d)", BAND_LOW_MAX + 1, BAND_HIGH_MIN - 1),
                    sprintf("High (>=%d)", BAND_HIGH_MIN),
                    "Total"),
    `1-year cancers, n` = ifelse(suppress & tot > 0 & tot < SDC_MIN, "<10", as.character(tot)),
    check.names = FALSE)
  labs <- c(if (is.null(SEP_DATE_COL)) "<=30 days, n (%)" else "<=30 days or index admission, n (%)",
            ">30 days-3 months, n (%)",
            ">3-12 months, n (%)")
  for (j in 1:3)
    out[[labs[j]]] <- mapply(fmt_cell, cnt[, j], tot,
                             if (suppress) prim[, j] else FALSE,
                             if (suppress) hide[, j] & !prim[, j] else FALSE)
  out
}

table_s8_internal <- build_table(suppress = FALSE)
table_s8_release  <- build_table(suppress = TRUE)

cat("\n=== TABLE S8 (RELEASE, suppressed) ===\n")
print(table_s8_release, row.names = FALSE, right = FALSE)
cat(sprintf("\n[SDC] %d primary cell(s) <10, %d secondary cell(s) S.\n",
            sum(prim), sum(hide & !prim)))
cat("Check by hand that no hidden cell can still be derived across rows AND columns together.\n")

if (is.null(SEP_DATE_COL))
  cat("\nNOTE: no index-admission discharge date set (SEP_DATE_COL = NULL), so column 1 is\n",
      "      '<=30 days' only. Either set SEP_DATE_COL, or change the header to '<=30 days'.\n", sep = "")

write.csv(table_s8_release,  file.path(OUT_DIR, "Table_S8_timing_by_band_SANT.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table_s8_internal, file.path(INTERNAL_DIR, "Table_S8_timing_by_band_SANT_UNSUPPRESSED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\nSaved: Table_S8_timing_by_band_SANT.csv (release) and unsuppressed copy in",
    basename(INTERNAL_DIR), "\n")




























# =============================================================================
# PAnCReAS — KM CUMULATIVE INCIDENCE BY RISK BAND (SANT), STEP-FUNCTION CURVES
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
DATA_PATH     <- NULL      # e.g. "C:/Data/cmbspine_pancreatitis.rds"; NULL = object in session
OUT_DIR       <- path.expand("~/Desktop/PAnCReAS_Figure_KM")
CENSUS_DATE   <- as.Date("2022-01-01")
HORIZON_YRS   <- 1
BAND_LOW_MAX  <- 1         # Low      : score <= 1
BAND_HIGH_MIN <- 9         # High     : score >= 9  (Moderate 2-8)
SHOW_CI       <- TRUE      # 95% CI bands
RISK_TIMES_MO <- c(0, 3, 6, 9, 12)
DPI           <- 300
FIG_W_IN      <- 7.5
FIG_H_IN      <- 6

EXPECTED_N <- c(Low = 16, Moderate = 88, High = 47)   # 1-year cancers, main pipeline

BAND_LABELS  <- c(Low      = sprintf("Low (\u2264%d)", BAND_LOW_MAX),
                  Moderate = sprintf("Moderate (%d\u2013%d)", BAND_LOW_MAX + 1, BAND_HIGH_MIN - 1),
                  High     = sprintf("High (\u2265%d)", BAND_HIGH_MIN))
BAND_COLOURS <- c(Low = "#0D0887", Moderate = "#CC4778", High = "#F0A227")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
if (!is.null(DATA_PATH)) {
  cmbspine_pancreatitis <- if (grepl("\\.rds$", DATA_PATH, ignore.case = TRUE))
    readRDS(DATA_PATH) else read.csv(DATA_PATH, stringsAsFactors = FALSE)
}
if (!exists("cmbspine_pancreatitis"))
  stop("cmbspine_pancreatitis is not loaded. Load it, or set DATA_PATH.")

required_cols <- c("ap_date_clean", "pc_date_clean", "death_date_clean",
                   "k85_panctcancer", "k85_panctcancer_diffweeks", "nv_k86",
                   "sant_deaths_flag", "nv_k85_firstage", "sant_nv_sex",
                   "nv_k851", "nv_k852", "nv_k859", "nv_k85own", "nv_k85other",
                   "nv_k858", "nv_k853")
missing_cols <- setdiff(required_cols, names(cmbspine_pancreatitis))
if (length(missing_cols) > 0)
  stop("Missing column(s): ", paste(missing_cols, collapse = ", "),
       "\nCheck: grep('panctcancer', names(cmbspine_pancreatitis), value = TRUE)")


# -----------------------------------------------------------------------------
# 2. COHORT — same exclusions and outcome as the main pipeline (Step 0C/1)
# -----------------------------------------------------------------------------
cutoff <- CENSUS_DATE - 365

dat <- cmbspine_pancreatitis %>%
  mutate(
    ap_date    = as.Date(ap_date_clean,    format = "%Y-%m-%d"),
    pc_date    = as.Date(pc_date_clean,    format = "%Y-%m-%d"),
    death_date = as.Date(death_date_clean, format = "%Y-%m-%d"),
    cp_flag    = coalesce(nv_k86 == 1, FALSE),
    ca_flag    = coalesce(k85_panctcancer == 1, FALSE),
    prior_ca   = ca_flag & coalesce(k85_panctcancer_diffweeks < 0, FALSE),
    ca_nodate  = ca_flag & !prior_ca &
      (is.na(k85_panctcancer_diffweeks) | is.na(pc_date))
  ) %>%
  filter(!cp_flag, !is.na(ap_date), !prior_ca, !ca_nodate, ap_date <= cutoff) %>%
  mutate(
    event_full = as.integer(ca_flag & !is.na(pc_date) & k85_panctcancer_diffweeks >= 0),
    end_date   = case_when(
      event_full == 1                            ~ pc_date,
      sant_deaths_flag == 1 & !is.na(death_date) ~ death_date,
      TRUE                                       ~ CENSUS_DATE),
    surv_years_full = pmax(as.numeric(end_date - ap_date) / 365.25, 1 / 365.25),
    event      = as.integer(event_full == 1 & surv_years_full <= HORIZON_YRS),
    surv_years = pmin(surv_years_full, HORIZON_YRS),
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
  drop_na(event, surv_years, age_50_59, age_60_69, age_70plus, male,
          biliary_ap, alcohol_ap, unspecified_ap, other_ap)

cat(sprintf("Analytic cohort N = %d | 1-year cancers = %d\n", nrow(dat), sum(dat$event)))


# -----------------------------------------------------------------------------
# 3. REFIT COX, RE-DERIVE POINTS AND BANDS (Steps 3, 5, 7 rules)
# -----------------------------------------------------------------------------
pred_vars <- c("age_50_59", "age_60_69", "age_70plus", "male",
               "biliary_ap", "alcohol_ap", "unspecified_ap", "other_ap")
cox_fit   <- coxph(as.formula(paste("Surv(surv_years, event) ~",
                                    paste(pred_vars, collapse = " + "))), data = dat)
b         <- coef(cox_fit)
ref_unit  <- min(abs(b[abs(b) > 0.10]))
point_map <- setNames(as.integer(round(b / ref_unit)), names(b))
cat("Point map:\n"); print(point_map)

dat <- dat %>%
  mutate(score = as.integer(as.matrix(dat[, pred_vars]) %*% point_map[pred_vars]),
         band  = factor(case_when(score <= BAND_LOW_MAX  ~ "Low",
                                  score >= BAND_HIGH_MIN ~ "High",
                                  TRUE                   ~ "Moderate"),
                        levels = c("Low", "Moderate", "High")))

got_n <- tapply(dat$event, dat$band, sum)
cat("\n1-year cancers by band:\n"); print(got_n)
if (!all(got_n[names(EXPECTED_N)] == EXPECTED_N))
  warning("Band event counts differ from 16 / 88 / 47 - the cohort or point map ",
          "does not match the main pipeline. Resolve before using the figure.")


# -----------------------------------------------------------------------------
# 4. KM BY BAND -> STEP-FUNCTION DATA
#    Each band's curve starts at (0, 0) and steps at every event time; CI
#    bands step with it. Time is plotted in months.
# -----------------------------------------------------------------------------
km <- survfit(Surv(surv_years, event) ~ band, data = dat)

km_df <- bind_rows(lapply(names(km$strata), function(s) {
  idx  <- if (s == names(km$strata)[1]) seq_len(km$strata[1]) else
    (sum(km$strata[seq_len(which(names(km$strata) == s) - 1)]) + 1):
    sum(km$strata[seq_len(which(names(km$strata) == s))])
  band <- sub("^band=", "", s)
  tibble(band  = band,
         time  = c(0, km$time[idx]),
         cuminc = 100 * (1 - c(1, km$surv[idx])),
         lo     = 100 * (1 - c(1, km$upper[idx])),
         hi     = 100 * (1 - c(1, km$lower[idx])))
})) %>%
  mutate(lo = ifelse(is.na(lo), cuminc, lo),
         hi = ifelse(is.na(hi), cuminc, hi)) %>%
  group_by(band) %>%
  # extend each curve flat to the end of follow-up
  group_modify(~ bind_rows(.x, tail(.x, 1) %>% mutate(time = HORIZON_YRS))) %>%
  ungroup() %>%
  mutate(months   = time * 12,
         band_lab = factor(unname(BAND_LABELS[band]), levels = unname(BAND_LABELS)))

# Step-shaped CI ribbon: duplicate each point so the ribbon is flat between events
ribbon_df <- km_df %>%
  group_by(band_lab) %>%
  arrange(months, .by_group = TRUE) %>%
  reframe(months = c(rbind(months, lead(months)))[-2 * n()],
          lo     = rep(lo, each = 2)[-2 * n()],
          hi     = rep(hi, each = 2)[-2 * n()])

# Check: curve at 3/6/9/12 months must match summary.survfit
chk <- summary(km, times = RISK_TIMES_MO[-1] / 12, extend = TRUE)
cat("\nKM cumulative incidence (%) at 3/6/9/12 months (check against Figure/Table 3):\n")
print(data.frame(band = sub("band=", "", chk$strata), months = round(chk$time * 12),
                 cuminc = round(100 * (1 - chk$surv), 2),
                 lo95 = round(100 * (1 - chk$upper), 2), hi95 = round(100 * (1 - chk$lower), 2)),
      row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. NUMBERS AT RISK
# -----------------------------------------------------------------------------
nar <- summary(km, times = RISK_TIMES_MO / 12, extend = TRUE)
risk_tab <- tibble(band     = sub("band=", "", nar$strata),
                   months   = round(nar$time * 12),
                   n_risk   = nar$n.risk) %>%
  mutate(band_lab = factor(unname(BAND_LABELS[band]), levels = rev(unname(BAND_LABELS))),
         lab      = ifelse(n_risk > 0 & n_risk < 10, "<10", format(n_risk, big.mark = ",")))
write.csv(risk_tab %>% select(band, months, n_risk = lab),
          file.path(OUT_DIR, "Figure_KM_step_numbers_at_risk.csv"), row.names = FALSE)


# -----------------------------------------------------------------------------
# 6. PLOT (no title, no caption; single legend at bottom)
# -----------------------------------------------------------------------------
cols <- setNames(unname(BAND_COLOURS), unname(BAND_LABELS))
ymax <- max(km_df$hi[km_df$months <= 12], na.rm = TRUE)
ymax <- ceiling(ymax)

p_km <- ggplot(km_df, aes(x = months, y = cuminc, colour = band_lab)) +
  { if (SHOW_CI) geom_ribbon(data = ribbon_df,
                             aes(x = months, ymin = lo, ymax = hi, fill = band_lab),
                             inherit.aes = FALSE, alpha = 0.15, colour = NA) } +
  geom_step(linewidth = 0.9, direction = "hv") +
  scale_colour_manual(values = cols, name = "PAnCReAS band") +
  scale_fill_manual(values = cols, name = "PAnCReAS band") +
  scale_x_continuous(breaks = RISK_TIMES_MO,
                     expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = c(0, 12)) +
  scale_y_continuous(limits = c(0, ymax), breaks = scales::pretty_breaks(6),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Months since index AP",
       y = "Cumulative incidence of pancreatic cancer") +
  theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        legend.justification = "center",
        panel.grid.minor = element_blank())

p_risk <- ggplot(risk_tab, aes(x = months, y = band_lab, label = lab, colour = band_lab)) +
  geom_text(size = 3.4, show.legend = FALSE) +
  # band names in their own column, left of month 0
  geom_text(data = distinct(risk_tab, band_lab),
            aes(x = -0.9, y = band_lab, label = band_lab, colour = band_lab),
            hjust = 1, size = 3.4, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
  scale_colour_manual(values = setNames(unname(BAND_COLOURS), unname(BAND_LABELS))) +
  scale_x_continuous(breaks = RISK_TIMES_MO,
                     expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = NULL, y = NULL, subtitle = "Number at risk") +
  coord_cartesian(xlim = c(0, 12), clip = "off") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        plot.subtitle = element_text(size = 10, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(2, 22, 2, 115))

p_km <- p_km + theme(plot.margin = margin(5, 22, 2, 5))
fig  <- (p_km / p_risk) +
  plot_layout(heights = c(4, 1), guides = "collect") &
  theme(legend.position = "bottom", legend.justification = "center")

ggsave(file.path(OUT_DIR, "Figure_KM_step_by_PAnCReAS_band.png"), fig,
       width = FIG_W_IN, height = FIG_H_IN, dpi = DPI, bg = "white")
ggsave(file.path(OUT_DIR, "Figure_KM_step_by_PAnCReAS_band.tiff"), fig,
       width = FIG_W_IN, height = FIG_H_IN, dpi = DPI, bg = "white", compression = "lzw")
cat("\nSaved: Figure_KM_step_by_PAnCReAS_band.png / .tiff (", DPI, " dpi) to ", OUT_DIR, "\n", sep = "")

# Small-cell note: each step in the curve is one or more cancers. The Low band
# has 16 one-year cancers in total, so its early steps represent counts of 1-9.
cat("NOTE: step curves show individual event times; confirm with the custodian that\n",
    "      KM curves are acceptable for release (standard practice, but check the Low band).\n")























# =============================================================================
# TABLE 3 — LAST COLUMN: "Diagnosed <=30 d / <=3 mo, %; deaths <=1 y, n (%)"
# Same cohort, score and bands as Steps 0C-7. Paste at the end of
# PAnCReAS_COX_FINAL_SANT.R and run after Steps 0-0D (cohort builder and
# small-cell helpers). Uses point_map from Step 6, or cox_point_map.rds.
# =============================================================================

library(tidyverse)
OUT_DIR <- path.expand("~/Desktop/PAnCReAS_16-10")
if (!exists("point_map")) point_map <- readRDS(file.path(OUT_DIR, "cox_point_map.rds"))
HORIZON_YRS   <- 1
BAND_LOW_MAX  <- 1
BAND_HIGH_MIN <- 9
DAYS_30       <- 30
DAYS_3MO      <- 91
EXPECTED_N    <- c(Low = 16, Moderate = 88, High = 47)   # 1-year cancers (Table 3)

if (!exists("build_ap_cohort") || !exists("sdc_small"))
  stop("Run Steps 0-0D at the top of this script first (cohort builder, small-cell helpers).")

t3_dat <- build_ap_cohort(cmbspine_pancreatitis) %>%
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
                  levels = c("Low", "Moderate", "High")),
    days_to_ca  = ifelse(event == 1, as.numeric(pc_date - ap_date), NA_real_),
    died_1yr    = as.integer(coalesce(sant_deaths_flag == 1, FALSE) & !is.na(death_date) &
                               as.numeric(death_date - ap_date) <= 365)
  )

# ---- checks against the main pipeline ---------------------------------------
got_n <- tapply(t3_dat$event, t3_dat$band, sum)
cat("Cohort N:", nrow(t3_dat), "| 1-year cancers by band:\n"); print(got_n)
if (!all(got_n[names(EXPECTED_N)] == EXPECTED_N))
  warning("1-year cancers by band differ from 16 / 88 / 47 - cohort or point map does not match Table 3.")
if (exists("dat_cox"))
  cat("Matches dat_cox (N / 1-yr events):",
      nrow(t3_dat) == nrow(dat_cox) && sum(t3_dat$event) == sum(dat_cox$event), "\n")
stopifnot(all(t3_dat$days_to_ca[t3_dat$event == 1] >= 0),
          all(t3_dat$days_to_ca[t3_dat$event == 1] <= 366))

# ---- counts (unsuppressed) ---------------------------------------------------
t3_counts <- t3_dat %>%
  group_by(band) %>%
  summarise(N          = n(),
            ca_1yr     = sum(event),
            ca_le30    = sum(event == 1 & days_to_ca <= DAYS_30),
            ca_le3mo   = sum(event == 1 & days_to_ca <= DAYS_3MO),
            deaths_1yr = sum(died_1yr),
            deaths_1yr_no_ca = sum(died_1yr == 1 & event == 0),
            .groups = "drop")

cat("\n=== Unsuppressed counts (internal only) ===\n")
print(as.data.frame(t3_counts), row.names = FALSE)

# ---- suppression -------------------------------------------------------------
# For each count column: primary = 1-9; complement (band total minus count)
# 1-9 also reveals it; then column-wise secondary suppression because the
# column totals are published in Table S8 (single hidden band -> hide the
# next smallest band too).
suppress_col <- function(k, n_total) {
  hide <- sdc_small(k) | sdc_small(n_total - k)
  prim <- sdc_small(k)
  if (sum(hide) == 1) {
    cand <- which(!hide)
    hide[cand[which.min(k[cand])]] <- TRUE
  }
  list(prim = prim, hide = hide)
}

s30  <- suppress_col(t3_counts$ca_le30,    t3_counts$ca_1yr)
s3mo <- suppress_col(t3_counts$ca_le3mo,   t3_counts$ca_1yr)
sdth <- suppress_col(t3_counts$deaths_1yr, t3_counts$N)

pct_cell <- function(k, n, s) ifelse(s$prim, "<10", ifelse(s$hide, "S", sprintf("%.1f", 100 * k / n)))
npct_cell <- function(k, n, s) ifelse(s$prim, "<10", ifelse(s$hide, "S",
                                                            sprintf("%s (%.1f)", format(k, big.mark = ",", trim = TRUE), 100 * k / n)))

table3_lastcol <- t3_counts %>%
  transmute(
    Band          = as.character(band),
    le30_pct      = pct_cell(ca_le30,  ca_1yr, s30),
    le3mo_pct     = pct_cell(ca_le3mo, ca_1yr, s3mo),
    deaths_1yr    = npct_cell(deaths_1yr, N, sdth),
    Table3_cell =
      sprintf("%s / %s; %s", le30_pct, le3mo_pct, deaths_1yr)
  )

cat("\n=== TABLE 3 - LAST COLUMN (release, suppressed) ===\n")
print(as.data.frame(table3_lastcol), row.names = FALSE, right = FALSE)
cat("\nDeaths are all-cause within 1 year of the index AP (including patients who\n",
    "had a cancer diagnosis first). Deaths without a 1-year cancer, by band:\n", sep = "")
print(t3_counts %>% transmute(band, deaths_1yr_no_ca = sdc_lab(deaths_1yr_no_ca)) %>% as.data.frame(),
      row.names = FALSE)
cat("\nCheck by hand: the <=30 d and <=3 mo cells must stay consistent with the\n",
    "suppression used in Table S8 (e.g. High <=30 d is 'S' there too).\n", sep = "")

write.csv(t3_counts, file.path(INTERNAL_DIR, "Table3_timing_deaths_UNSUPPRESSED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table3_lastcol, file.path(OUT_DIR, "Table3_timing_deaths_column.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("Saved: Table3_timing_deaths_column.csv (release) and the unsuppressed copy in",
    basename(INTERNAL_DIR), "\n")






































# =============================================================================
# PAnCReAS — TABLE S7 (SANT ROWS): DIAGNOSTIC PERFORMANCE OF MODERATE + HIGH
#            AND HIGH-ONLY THRESHOLDS, CANCER WITHIN 1 YEAR — WILSON 95% CIs
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
})

# -----------------------------------------------------------------------------
# SETTINGS — check these before running
# -----------------------------------------------------------------------------
DATA_PATH     <- NULL      # e.g. "C:/Data/cmbspine_pancreatitis.rds"; NULL = object in session
OUT_DIR       <- path.expand("~/Desktop/PAnCReAS_TableS7")
MAIN_DIR      <- path.expand("~/Desktop/PAnCReAS_16-10")   # main pipeline outputs (for cross-checks)

CENSUS_DATE   <- as.Date("2022-01-01")
HORIZON_YRS   <- 1
BAND_LOW_MAX  <- 1         # Low      : score <= 1
BAND_HIGH_MIN <- 9         # High     : score >= 9  (Moderate 2-8)
SDC_MIN       <- 10        # custodian rule: counts 1-9 -> "<10"
Z             <- qnorm(0.975)

EXPECTED_N    <- c(Low = 16, Moderate = 88, High = 47)   # 1-year cancers by band (main pipeline)

# Current Table S7 SANT point estimates (%), Clopper-Pearson version.
# Point estimates must not change when only the CI method changes.
EXPECTED_PT <- rbind(
  `Moderate + High` = c(Sens = 89.4, Spec = 51.6, PPV = 2.0, NPV = 99.8),
  `High only`       = c(Sens = 31.1, Spec = 92.3, PPV = 4.3, NPV = 99.2)
)

INTERNAL_DIR <- file.path(OUT_DIR, "internal_unsuppressed_DO_NOT_RELEASE")
dir.create(INTERNAL_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. LOAD DATA
# -----------------------------------------------------------------------------
if (!is.null(DATA_PATH)) {
  cmbspine_pancreatitis <- if (grepl("\\.rds$", DATA_PATH, ignore.case = TRUE))
    readRDS(DATA_PATH) else read.csv(DATA_PATH, stringsAsFactors = FALSE)
}
if (!exists("cmbspine_pancreatitis"))
  stop("cmbspine_pancreatitis is not loaded. Load it, or set DATA_PATH.")

required_cols <- c("ap_date_clean", "pc_date_clean", "death_date_clean",
                   "k85_panctcancer", "k85_panctcancer_diffweeks", "nv_k86",
                   "sant_deaths_flag", "nv_k85_firstage", "sant_nv_sex",
                   "nv_k851", "nv_k852", "nv_k859", "nv_k85own", "nv_k85other",
                   "nv_k858", "nv_k853")
missing_cols <- setdiff(required_cols, names(cmbspine_pancreatitis))
if (length(missing_cols) > 0)
  stop("Missing column(s): ", paste(missing_cols, collapse = ", "))


# -----------------------------------------------------------------------------
# 2. COHORT — same exclusions and outcome as the main pipeline (Step 0C / 1)
#    1. nv_k86 == 1   2. missing first-AP date   3. cancer before first AP
#    4. cancer with no diagnosis date   5. first AP after CENSUS_DATE - 365
#    then complete-case on predictors
# -----------------------------------------------------------------------------
cutoff    <- CENSUS_DATE - 365
pred_vars <- c("age_50_59", "age_60_69", "age_70plus", "male",
               "biliary_ap", "alcohol_ap", "unspecified_ap", "other_ap")

dat <- cmbspine_pancreatitis %>%
  mutate(
    ap_date    = as.Date(ap_date_clean,    format = "%Y-%m-%d"),
    pc_date    = as.Date(pc_date_clean,    format = "%Y-%m-%d"),
    death_date = as.Date(death_date_clean, format = "%Y-%m-%d"),
    cp_flag    = coalesce(nv_k86 == 1, FALSE),
    ca_flag    = coalesce(k85_panctcancer == 1, FALSE),
    prior_ca   = ca_flag & coalesce(k85_panctcancer_diffweeks < 0, FALSE),
    ca_nodate  = ca_flag & !prior_ca &
      (is.na(k85_panctcancer_diffweeks) | is.na(pc_date))
  ) %>%
  filter(!cp_flag, !is.na(ap_date), !prior_ca, !ca_nodate, ap_date <= cutoff) %>%
  mutate(
    event_full = as.integer(ca_flag & !is.na(pc_date) & k85_panctcancer_diffweeks >= 0),
    end_date   = case_when(
      event_full == 1                            ~ pc_date,
      sant_deaths_flag == 1 & !is.na(death_date) ~ death_date,
      TRUE                                       ~ CENSUS_DATE),
    surv_years_full = pmax(as.numeric(end_date - ap_date) / 365.25, 1 / 365.25),
    event      = as.integer(event_full == 1 & surv_years_full <= HORIZON_YRS),
    surv_years = pmin(surv_years_full, HORIZON_YRS),
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
  drop_na(event, surv_years, event_full, all_of(pred_vars))

cat(sprintf("Analytic cohort N = %d | 1-year cancers = %d\n\n", nrow(dat), sum(dat$event)))


# -----------------------------------------------------------------------------
# 3. REFIT 1-YEAR COX, RE-DERIVE INTEGER POINTS, ASSIGN BANDS
#    Reference unit = smallest |log-HR| > 0.10; points = round(log-HR / unit)
# -----------------------------------------------------------------------------
cox_fit   <- coxph(as.formula(paste("Surv(surv_years, event) ~",
                                    paste(pred_vars, collapse = " + "))), data = dat)
b         <- coef(cox_fit)
ref_unit  <- min(abs(b[abs(b) > 0.10]))
point_map <- setNames(as.integer(round(b / ref_unit)), names(b))
cat("Point map (re-derived):\n"); print(point_map)

saved_pm <- file.path(MAIN_DIR, "cox_point_map.rds")
if (file.exists(saved_pm)) {
  pm0 <- readRDS(saved_pm)
  ok  <- isTRUE(all(pm0[names(point_map)] == point_map))
  cat("Matches saved cox_point_map.rds:", ok, "\n")
  if (!ok) warning("Point map differs from the main pipeline - resolve before using this table.")
}

dat <- dat %>%
  mutate(score = as.integer(as.matrix(dat[, pred_vars]) %*% point_map[pred_vars]),
         band  = factor(case_when(score <= BAND_LOW_MAX  ~ "Low",
                                  score >= BAND_HIGH_MIN ~ "High",
                                  TRUE                   ~ "Moderate"),
                        levels = c("Low", "Moderate", "High")))

got_n <- tapply(dat$event, dat$band, sum)
cat("\n1-year cancers by band:\n"); print(got_n)
if (!all(got_n[names(EXPECTED_N)] == EXPECTED_N))
  warning("Band event counts differ from 16 / 88 / 47 - cohort or point map does not ",
          "match the main pipeline. Resolve before using this table.")


# -----------------------------------------------------------------------------
# 4. WILSON SCORE INTERVAL (no continuity correction)
#    Equivalent to prop.test(x, n, correct = FALSE)$conf.int; checked below.
# -----------------------------------------------------------------------------
wilson_ci <- function(x, n, z = Z) {
  p     <- x / n
  denom <- 1 + z^2 / n
  ctr   <- (p + z^2 / (2 * n)) / denom
  half  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(lo = max(0, ctr - half), hi = min(1, ctr + half))
}
stopifnot(isTRUE(all.equal(unname(wilson_ci(17, 40)),
                           as.numeric(prop.test(17, 40, correct = FALSE)$conf.int),
                           tolerance = 1e-10)))

sdc_small <- function(x) !is.na(x) & x > 0 & x < SDC_MIN


# -----------------------------------------------------------------------------
# 5. DIAGNOSTIC PERFORMANCE BY THRESHOLD
#    Patients censored (death) before 1 year without cancer are classified
#    negative, as in the main pipeline Step 9A.
#    Suppression rules identical to Step 9A:
#      sens/spec -> "S" if numerator or its complement is 1-9
#      PPV -> "S" if FN is 1-9; NPV -> "S" if TP is 1-9
# -----------------------------------------------------------------------------
y <- dat$event

perf_counts <- function(label, flag) {
  tibble(Threshold = label,
         TP = sum(flag & y == 1), FN = sum(!flag & y == 1),
         FP = sum(flag & y == 0), TN = sum(!flag & y == 0))
}

counts <- bind_rows(
  perf_counts(sprintf("SANT: Moderate + High (\u2265%d)", BAND_LOW_MAX + 1),
              dat$score >= BAND_LOW_MAX + 1),
  perf_counts(sprintf("SANT: High only (\u2265%d)", BAND_HIGH_MIN),
              dat$score >= BAND_HIGH_MIN)
)

fmt_est <- function(x, n, suppress = FALSE) {
  if (suppress) return("S")
  ci <- wilson_ci(x, n) * 100
  sprintf("%.1f (%.1f\u2013%.1f)", 100 * x / n, ci["lo"], ci["hi"])
}

build_s7 <- function(cnt, apply_sdc) {
  cnt %>%
    rowwise() %>%
    transmute(
      Threshold,
      `Sensitivity (95% CI)` = fmt_est(TP, TP + FN, apply_sdc && (sdc_small(TP) || sdc_small(FN))),
      `Specificity (95% CI)` = fmt_est(TN, TN + FP, apply_sdc && (sdc_small(TN) || sdc_small(FP))),
      `PPV (95% CI)`         = fmt_est(TP, TP + FP, apply_sdc && (sdc_small(FN) || sdc_small(TP) || sdc_small(FP))),
      `NPV (95% CI)`         = fmt_est(TN, TN + FN, apply_sdc && (sdc_small(TP) || sdc_small(FN) || sdc_small(TN)))
    ) %>%
    ungroup()
}

table_s7_release  <- build_s7(counts, apply_sdc = TRUE)
table_s7_internal <- build_s7(counts, apply_sdc = FALSE) %>%
  bind_cols(counts %>% select(TP, FN, FP, TN))


# -----------------------------------------------------------------------------
# 6. CHECKS
#    (a) point estimates unchanged from the current Table S7
#    (b) side-by-side Wilson vs Clopper-Pearson, for the record
# -----------------------------------------------------------------------------
pt <- counts %>%
  transmute(Sens = round(100 * TP / (TP + FN), 1), Spec = round(100 * TN / (TN + FP), 1),
            PPV  = round(100 * TP / (TP + FP), 1), NPV  = round(100 * TN / (TN + FN), 1)) %>%
  as.matrix()
pt_ok <- isTRUE(all(abs(pt - EXPECTED_PT) < 0.051))
cat("\nPoint estimates match current Table S7:", pt_ok, "\n")
if (!pt_ok) {
  print(cbind(new = pt, current = EXPECTED_PT))
  warning("Point estimates differ from the current Table S7 - cohort or score has changed.")
}

cp_vs_wilson <- bind_rows(lapply(seq_len(nrow(counts)), function(i) {
  r <- counts[i, ]
  m <- list(Sensitivity = c(r$TP, r$TP + r$FN), Specificity = c(r$TN, r$TN + r$FP),
            PPV = c(r$TP, r$TP + r$FP),         NPV = c(r$TN, r$TN + r$FN))
  bind_rows(lapply(names(m), function(k) {
    x <- m[[k]][1]; n <- m[[k]][2]
    cp <- binom.test(x, n)$conf.int * 100; wi <- wilson_ci(x, n) * 100
    tibble(Threshold = r$Threshold, Metric = k, Est = round(100 * x / n, 1),
           CP = sprintf("%.1f\u2013%.1f", cp[1], cp[2]),
           Wilson = sprintf("%.1f\u2013%.1f", wi["lo"], wi["hi"]))
  }))
}))
cat("\nClopper-Pearson vs Wilson (internal check, %):\n")
print(as.data.frame(cp_vs_wilson), row.names = FALSE)


# -----------------------------------------------------------------------------
# 7. OUTPUT
# -----------------------------------------------------------------------------
cat("\n=== TABLE S7, SANT ROWS (RELEASE, Wilson 95% CI) ===\n")
print(as.data.frame(table_s7_release), row.names = FALSE, right = FALSE)

write.csv(table_s7_release, file.path(OUT_DIR, "Table_S7_SANT_Wilson.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(table_s7_internal, file.path(INTERNAL_DIR, "Table_S7_SANT_Wilson_UNSUPPRESSED.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cp_vs_wilson, file.path(INTERNAL_DIR, "Table_S7_SANT_CP_vs_Wilson.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\nSaved: Table_S7_SANT_Wilson.csv (release) to", OUT_DIR, "\n")
cat("       unsuppressed copy and CP-vs-Wilson check in", basename(INTERNAL_DIR), "\n")
cat("\nFootnote edit: replace 'exact (Clopper-Pearson) 95% CIs in SANT and Wilson score\n",
    "95% CIs in VIC and TAS' with 'Wilson score 95% CIs'.\n", sep = "")



















