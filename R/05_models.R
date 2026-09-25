# ---------------------------------------------------------------------------
# 05_models.R
#
# Vaccine effectiveness estimates from the matched cohort.
#
# Three corrections relative to the original thesis code:
#
#   1. Confidence intervals are computed from the variance matrix of the model
#      being reported. The original code estimated a crude model, stored its
#      variance matrix, and then reused that matrix for the adjusted model, so
#      the adjusted intervals belonged to a different model.
#
#   2. Effectiveness by follow-up window truncates person-time at the end of
#      the window, not only the event indicator. Redefining the event as
#      "infection before day t" while leaving follow-up at 360 days changes the
#      numerator without changing the denominator, which is not a risk over
#      that window.
#
#   3. Subgroup models are fitted by looping over factor levels with explicit
#      bookkeeping, and estimates from subgroups with too few events are
#      reported as such rather than silently returned.
#
# Because both members of a matched pair enter on the same calendar day, time
# since matching and calendar time are aligned within pair; models are
# stratified by pair, which gives each pair its own baseline hazard.
#
# Input:  data/derived/matched_pairs.rds, data/derived/cohort_baseline.rds
# Output: data/derived/ve_results.rds
#         figures/ve_by_window.png, figures/ve_subgroups.png
# ---------------------------------------------------------------------------

source("R/00_functions.R")
library(survival)

OUT_DIR <- "data/derived"
FIG_DIR <- "figures"
MIN_EVENTS <- 20      # below this, an estimate is reported as unstable

m <- readRDS(file.path(OUT_DIR, "matched_pairs.rds"))
d <- readRDS(file.path(OUT_DIR, "cohort_baseline.rds"))

# --- Severe outcomes --------------------------------------------------------
# Hospitalisation and death are defined among those infected, dated at symptom
# onset. An individual contributes an event only if the infection that led to
# the outcome occurred inside their follow-up interval.

sev <- d[, c("id", "hospitalised", "died")]
names(sev) <- c("id", "hosp_flag", "death_flag")
m <- merge(m, sev, by = "id", all.x = TRUE, suffixes = c("", ".base"))

m$event_infection     <- m$event
m$event_hospitalised  <- as.integer(m$event == 1 & m$hosp_flag  == 1)
m$event_death         <- as.integer(m$event == 1 & m$death_flag == 1)

# The vaccine product is recorded for the vaccinated member only; propagate it
# to the pair so that subgroup analyses by product are possible.
pair_vaccine <- tapply(m$vaccine[m$exposed == 1], m$pair[m$exposed == 1],
                       function(x) x[1])
m$pair_vaccine <- pair_vaccine[as.character(m$pair)]

# --- Helper: fit one model and return a tidy row ---------------------------

fit_ve <- function(data, outcome, label, group = "Overall") {
  n_events <- sum(data[[outcome]], na.rm = TRUE)
  if (n_events < MIN_EVENTS || length(unique(data$exposed)) < 2) {
    return(data.frame(outcome = outcome, analysis = label, group = group,
                      n = nrow(data), events = n_events,
                      ve = NA, ve_low = NA, ve_high = NA,
                      note = "too few events", stringsAsFactors = FALSE))
  }
  f <- try(coxph(
    Surv(followup_days, data[[outcome]]) ~ exposed + strata(pair),
    data = data
  ), silent = TRUE)
  if (inherits(f, "try-error") || !"exposed" %in% names(coef(f))) {
    return(data.frame(outcome = outcome, analysis = label, group = group,
                      n = nrow(data), events = n_events,
                      ve = NA, ve_low = NA, ve_high = NA,
                      note = "model did not converge", stringsAsFactors = FALSE))
  }
  v <- ve_from_cox(f, "exposed")
  data.frame(outcome = outcome, analysis = label, group = group,
             n = nrow(data), events = n_events,
             ve = v$ve, ve_low = v$ve_low, ve_high = v$ve_high,
             note = "", stringsAsFactors = FALSE)
}

results <- list()

# --- Primary analysis: three outcomes --------------------------------------

for (out in c("event_infection", "event_hospitalised", "event_death")) {
  results[[length(results) + 1]] <- fit_ve(m, out, "Primary")
}

# --- Effectiveness by follow-up window -------------------------------------
# Person-time is truncated at the end of each window; an individual who was
# infected after the window closes is treated as event-free within it.

windows <- c(60, 120, 180, 240, 300, 360)

for (w in windows) {
  mw <- m
  inside <- mw$followup_days <= w
  mw$followup_days <- pmin(mw$followup_days, w)
  mw$event_infection <- as.integer(mw$event_infection == 1 & inside)
  mw <- mw[mw$followup_days > 0, ]
  results[[length(results) + 1]] <-
    fit_ve(mw, "event_infection", "By follow-up window",
           group = sprintf("0-%d days", w))
}

# --- Subgroups --------------------------------------------------------------

subgroup_vars <- list(
  "Age group"         = "age_group",
  "Sex"               = "sex",
  "Comorbidity group" = "comorbidity",
  "City"              = "city",
  "Vaccine product"   = "pair_vaccine"
)

for (nm in names(subgroup_vars)) {
  v <- subgroup_vars[[nm]]
  levs <- sort(unique(as.character(m[[v]])))
  levs <- levs[!is.na(levs)]
  for (l in levs) {
    sub <- m[!is.na(m[[v]]) & as.character(m[[v]]) == l, ]
    results[[length(results) + 1]] <-
      fit_ve(sub, "event_infection", nm, group = l)
  }
}

ve_results <- do.call(rbind, results)
saveRDS(ve_results, file.path(OUT_DIR, "ve_results.rds"))

# --- Figures ----------------------------------------------------------------

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

forest <- function(df, title, file) {
  df <- df[!is.na(df$ve), ]
  if (!nrow(df)) return(invisible(NULL))
  png(file, width = 1500, height = 170 + 72 * nrow(df), res = 170)
  op <- par(mar = c(5, 14, 4, 2))
  yy <- rev(seq_len(nrow(df)))
  xlo <- min(0, floor(min(100 * df$ve_low) / 10) * 10)
  plot(NA, xlim = c(xlo, 100), ylim = c(0.5, nrow(df) + 0.5),
       yaxt = "n", ylab = "", xlab = "Vaccine effectiveness (%)", main = title)
  abline(v = seq(xlo, 100, 20), col = "grey90"); abline(v = 0, col = "grey50")
  segments(100 * df$ve_low, yy, 100 * df$ve_high, yy, lwd = 2)
  points(100 * df$ve, yy, pch = 19, cex = 1.1)
  axis(2, at = yy, labels = df$group, las = 1, cex.axis = 0.85)
  par(op); invisible(dev.off())
}

forest(ve_results[ve_results$analysis == "By follow-up window", ],
       "Effectiveness against infection, by follow-up window",
       file.path(FIG_DIR, "ve_by_window.png"))

sub_rows <- ve_results[ve_results$analysis %in% names(subgroup_vars), ]
sub_rows$group <- paste0(sub_rows$analysis, ": ", sub_rows$group)
forest(sub_rows, "Effectiveness against infection, by subgroup",
       file.path(FIG_DIR, "ve_subgroups.png"))

# --- Report -----------------------------------------------------------------

pretty <- function(r) {
  if (is.na(r$ve)) return(sprintf("%-28s  %6d events  %s", r$group, r$events, r$note))
  sprintf("%-28s  %6d events  VE %5.1f%% (95%% CI %5.1f to %5.1f)",
          r$group, r$events, 100 * r$ve, 100 * r$ve_low, 100 * r$ve_high)
}

cat("\nPrimary analysis\n----------------\n")
prim <- ve_results[ve_results$analysis == "Primary", ]
labels <- c(event_infection = "Confirmed infection",
            event_hospitalised = "Hospitalisation",
            event_death = "Death")
for (i in seq_len(nrow(prim))) {
  r <- prim[i, ]; r$group <- labels[[r$outcome]]
  cat(pretty(r), "\n")
}

for (sec in c("By follow-up window", names(subgroup_vars))) {
  rows <- ve_results[ve_results$analysis == sec, ]
  if (!nrow(rows)) next
  cat(sprintf("\n%s\n%s\n", sec, strrep("-", nchar(sec))))
  for (i in seq_len(nrow(rows))) cat(pretty(rows[i, ]), "\n")
}

cat(sprintf("\nResults saved to %s\nFigures written to %s\n",
            file.path(OUT_DIR, "ve_results.rds"), FIG_DIR))
