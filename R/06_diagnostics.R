# ---------------------------------------------------------------------------
# 06_diagnostics.R
#
# Assumption checks and sensitivity analyses.
#
# Three things are done here:
#
#   1. The proportional hazards assumption is tested formally (scaled Schoenfeld
#      residuals) and inspected graphically, rather than only by eye.
#
#   2. The primary design is compared with two alternatives that differ only in
#      how time zero and the unexposed group are defined. On synthetic data,
#      where the true effect is known, this shows what each design recovers.
#      On real data it shows how much the answer depends on those choices,
#      which is the honest thing to report.
#
#   3. The sensitivity of the estimate to unmeasured confounding is quantified
#      with an E-value (VanderWeele & Ding, Ann Intern Med 2017): the minimum
#      strength of association an unmeasured confounder would need with both
#      vaccination and infection to explain away the result.
#
# Input:  data/derived/matched_pairs.rds
#         data/derived/cohort_person_time.rds
#         data/derived/cohort_baseline.rds
# Output: data/derived/diagnostics.rds
#         figures/ph_loglog.png, figures/ph_schoenfeld.png,
#         figures/design_comparison.png
# ---------------------------------------------------------------------------

source("R/00_functions.R")
library(survival)

OUT_DIR     <- "data/derived"
FIG_DIR     <- "figures"
STUDY_END   <- as.Date("2022-06-30")
MAX_FOLLOWUP <- 360

m  <- readRDS(file.path(OUT_DIR, "matched_pairs.rds"))
pt <- readRDS(file.path(OUT_DIR, "cohort_person_time.rds"))
d  <- readRDS(file.path(OUT_DIR, "cohort_baseline.rds"))

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# --- 1. Proportional hazards -----------------------------------------------

ph_fit <- coxph(Surv(followup_days, event) ~ exposed + cluster(pair), data = m)
ph_test <- cox.zph(ph_fit)

cat("Proportional hazards test (scaled Schoenfeld residuals)\n")
print(ph_test)

png(file.path(FIG_DIR, "ph_schoenfeld.png"), width = 1400, height = 900, res = 170)
plot(ph_test, main = "Scaled Schoenfeld residuals for the vaccination term")
abline(h = coef(ph_fit)[["exposed"]], col = "grey40", lty = 2)
invisible(dev.off())

# Complementary log-log plot: parallel curves support proportional hazards
km <- survfit(Surv(followup_days, event) ~ exposed, data = m)
png(file.path(FIG_DIR, "ph_loglog.png"), width = 1400, height = 900, res = 170)
plot(km, fun = "cloglog", lty = c(1, 2), lwd = 2,
     xlab = "Days since matching (log scale)",
     ylab = "log(-log(S(t)))",
     main = "Complementary log-log plot by treatment strategy")
legend("bottomright", c("Unvaccinated", "Complete primary schedule"),
       lty = c(1, 2), lwd = 2, bty = "n")
invisible(dev.off())

# --- 2. Comparison of designs ----------------------------------------------

designs <- list()

# (a) Primary: risk-set matched pairs, stratified by pair
f_matched <- coxph(Surv(followup_days, event) ~ exposed + strata(pair), data = m)
v <- ve_from_cox(f_matched, "exposed")
designs[["Risk-set matched (primary)"]] <- v

# (b) Whole cohort, time-varying exposure, calendar time scale
f_tv <- coxph(
  Surv(tstart, tstop, event) ~ exposed + age_group + sex + comorbidity +
    city + regimen + cluster(id),
  data = pt
)
designs[["Whole cohort, time-varying, calendar time"]] <- ve_from_cox(f_tv, "exposed")

# (c) Whole cohort, time-varying exposure, but time measured from each
#     individual's own origin: ignores that exposed person-time falls later in
#     the epidemic
f_origin <- coxph(
  Surv(followup_days, event) ~ exposed + age_group + sex + comorbidity +
    city + regimen + cluster(id),
  data = pt
)
designs[["Whole cohort, time since individual origin"]] <- ve_from_cox(f_origin, "exposed")

# (d) The original thesis design: unexposed defined as never vaccinated during
#     the whole study period, time measured from each individual's origin
nv <- d[d$analysis_set %in% c("exposed", "never vaccinated"), ]
nv$exposed <- as.integer(nv$analysis_set == "exposed")
nv$t0 <- as.Date(ifelse(nv$exposed == 1,
                        as.character(nv$exposure_start),
                        as.character(nv$eligibility_date)))
nv$event_flag <- as.integer(!is.na(nv$event_date) & nv$event_date > nv$t0)
nv$time <- ifelse(nv$event_flag == 1,
                  as.numeric(nv$event_date - nv$t0),
                  as.numeric(STUDY_END - nv$t0))
nv <- nv[!is.na(nv$time) & nv$time > 0, ]
nv$event_flag[nv$time > MAX_FOLLOWUP] <- 0
nv$time <- pmin(nv$time, MAX_FOLLOWUP)

f_naive <- coxph(
  Surv(time, event_flag) ~ exposed + age_group + sex + comorbidity + city + regimen,
  data = nv
)
designs[["Never-vaccinated comparison (original)"]] <- ve_from_cox(f_naive, "exposed")

comparison <- do.call(rbind, lapply(names(designs), function(nm) {
  cbind(design = nm, designs[[nm]][, c("hr", "ve", "ve_low", "ve_high")])
}))

cat("\nComparison of designs\n---------------------\n")
for (i in seq_len(nrow(comparison))) {
  cat(sprintf("%-44s VE %5.1f%% (95%% CI %5.1f to %5.1f)\n",
              comparison$design[i], 100 * comparison$ve[i],
              100 * comparison$ve_low[i], 100 * comparison$ve_high[i]))
}
cat("\nThe first two designs align exposed and unexposed person-time in calendar\n",
    "time; the last two do not, and are shown to quantify what that costs.\n", sep = "")

png(file.path(FIG_DIR, "design_comparison.png"),
    width = 1600, height = 170 + 80 * nrow(comparison), res = 170)
op <- par(mar = c(5, 20, 4, 2))
yy <- rev(seq_len(nrow(comparison)))
plot(NA, xlim = c(40, 80), ylim = c(0.5, nrow(comparison) + 0.5),
     yaxt = "n", ylab = "", xlab = "Vaccine effectiveness (%)",
     main = "Effectiveness by analytic design")
abline(v = seq(40, 80, 5), col = "grey92")
segments(100 * comparison$ve_low, yy, 100 * comparison$ve_high, yy, lwd = 2)
points(100 * comparison$ve, yy, pch = 19, cex = 1.1)
axis(2, at = yy, labels = comparison$design, las = 1, cex.axis = 0.7)
par(op); invisible(dev.off())

# --- 3. Sensitivity to unmeasured confounding -------------------------------
# E-value on the hazard ratio scale, using the square-root transformation
# recommended for common outcomes.

evalue <- function(hr, hr_low, hr_high) {
  # Strip names: confint() returns a named vector and the names would
  # propagate into the result, breaking lookup by element name.
  hr <- unname(hr); hr_low <- unname(hr_low); hr_high <- unname(hr_high)

  # Square-root transformation of the hazard ratio, recommended when the
  # outcome is common.
  rr <- sqrt(hr)

  ev_one <- function(x) {
    if (x == 1) return(1)
    b <- if (x < 1) 1 / x else x
    b + sqrt(b * (b - 1))
  }

  # The confidence limit closest to the null is the one that matters.
  bound <- if (rr < 1) sqrt(hr_high) else sqrt(hr_low)
  ci <- if ((rr < 1 && bound >= 1) || (rr > 1 && bound <= 1)) 1 else ev_one(bound)

  c(point = ev_one(rr), ci = ci)
}

primary_hr <- exp(coef(f_matched)[["exposed"]])
primary_ci <- exp(confint(f_matched)["exposed", ])
ev <- evalue(primary_hr, primary_ci[1], primary_ci[2])

cat(sprintf(
"\nSensitivity to unmeasured confounding
-------------------------------------
Primary hazard ratio                      %.3f (95%% CI %.3f to %.3f)
E-value for the point estimate            %.2f
E-value for the confidence limit          %.2f

An unmeasured confounder would have to be associated with both vaccination and
infection by a risk ratio of at least %.2f, above and beyond the matched
covariates, to explain away the estimate entirely.
",
primary_hr, primary_ci[1], primary_ci[2], ev[["point"]], ev[["ci"]], ev[["point"]]))

# --- Save -------------------------------------------------------------------

diagnostics <- list(
  ph_test    = ph_test,
  comparison = comparison,
  evalue     = ev
)
saveRDS(diagnostics, file.path(OUT_DIR, "diagnostics.rds"))

cat(sprintf("\nDiagnostics saved to %s\nFigures written to %s\n",
            file.path(OUT_DIR, "diagnostics.rds"), FIG_DIR))
