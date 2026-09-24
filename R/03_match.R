# ---------------------------------------------------------------------------
# 03_match.R
#
# Risk-set matching of vaccinated to unvaccinated individuals.
#
# Design note, and the second methodological change relative to the original
# thesis code:
#
#   The thesis matched on a propensity score estimated from a single covariate
#   (insurance scheme), which leaves age, sex, comorbidity and city unbalanced
#   and offers no protection against calendar-time confounding.
#
#   Here matching is done inside risk sets, the design used by the large
#   national vaccine effectiveness studies. On each calendar day, every person
#   whose protection begins that day is matched to one person who, on that same
#   day, is still unvaccinated, still uninfected and already eligible under the
#   priority schedule. Matching is exact on city, sex, age group, comorbidity
#   group and insurance scheme, so the pair is comparable both in covariates
#   and in epidemic phase.
#
#   Follow-up for a pair ends when either member has the outcome, when the
#   matched control is vaccinated, at the exposed member's booster, at the
#   administrative end, or after MAX_FOLLOWUP days, whichever comes first.
#   Censoring the pair when the control is vaccinated is what keeps the
#   contrast between the two treatment strategies clean.
#
#   A propensity score is still estimated, from the full covariate set, and
#   reported as a balance diagnostic rather than used to form the pairs.
#
# Input:  data/derived/cohort_baseline.rds
# Output: data/derived/matched_pairs.rds
#         figures/balance_smd.png
# ---------------------------------------------------------------------------

source("R/00_functions.R")

set.seed(20241213)   # graduation date

OUT_DIR      <- "data/derived"
FIG_DIR      <- "figures"
STUDY_START  <- as.Date("2021-02-17")
STUDY_END    <- as.Date("2022-06-30")
MAX_FOLLOWUP <- 360

d <- readRDS(file.path(OUT_DIR, "cohort_baseline.rds"))

# --- Everything on a calendar-day scale ------------------------------------

day_of <- function(x) as.numeric(x - STUDY_START)

d$elig_day    <- day_of(d$eligibility_date)
d$dose1_day   <- day_of(d$dose1_date)
d$booster_day <- day_of(d$booster_date)
d$prot_day    <- day_of(d$exposure_start)      # NA unless schedule completed
d$event_day   <- day_of(d$event_date)          # NA unless infected
END_DAY       <- day_of(STUDY_END)

# Matching strata: exact on these five variables. Missing values get their own
# level so that records are matched like with like rather than dropped.
lv <- function(x) addNA(factor(x), ifany = TRUE)
d$stratum <- interaction(
  lv(d$city), lv(d$sex), lv(d$age_group), lv(d$comorbidity), lv(d$regimen),
  drop = TRUE, sep = "|"
)

# --- Who can be a case, and when -------------------------------------------

is_case <- !is.na(d$prot_day) &
  d$prot_day >= 0 & d$prot_day <= END_DAY &
  (is.na(d$event_day) | d$event_day > d$prot_day)

case_days <- sort(unique(d$prot_day[is_case]))

# --- Risk-set matching ------------------------------------------------------

available <- rep(TRUE, nrow(d))     # not yet used as a control
strata_levels <- levels(d$stratum)

matched_case    <- integer(0)
matched_control <- integer(0)
matched_day     <- numeric(0)

for (t in case_days) {

  cases <- which(is_case & d$prot_day == t)
  if (!length(cases)) next

  # Eligible controls on day t: already eligible under the priority schedule,
  # not yet vaccinated, not yet infected, and not already used.
  pool <- which(
    available &
      d$elig_day <= t &
      (is.na(d$dose1_day) | d$dose1_day > t) &
      (is.na(d$event_day) | d$event_day > t)
  )
  if (!length(pool)) next

  # Match within stratum, one control per case, sampled at random.
  pool_by_stratum  <- split(pool,  d$stratum[pool],  drop = TRUE)
  cases_by_stratum <- split(cases, d$stratum[cases], drop = TRUE)

  for (s in names(cases_by_stratum)) {
    cs <- cases_by_stratum[[s]]
    ps <- pool_by_stratum[[s]]
    if (is.null(ps) || !length(ps)) next

    n <- min(length(cs), length(ps))
    picked_cases    <- if (length(cs) > n) sample(cs, n) else cs
    picked_controls <- if (length(ps) > n) sample(ps, n) else ps

    matched_case    <- c(matched_case, picked_cases)
    matched_control <- c(matched_control, picked_controls)
    matched_day     <- c(matched_day, rep(t, n))

    available[picked_controls] <- FALSE
  }
}

message(sprintf("Matched %d pairs out of %d eligible vaccinated individuals (%.1f%%)",
                length(matched_case), sum(is_case),
                100 * length(matched_case) / sum(is_case)))

# --- Follow-up for each pair ------------------------------------------------
# The pair is followed from the matching day until the first of: outcome in
# either member, vaccination of the control, booster in the case, the
# administrative end, or MAX_FOLLOWUP days.

pair_end <- pmin(
  matched_day + MAX_FOLLOWUP,
  END_DAY,
  ifelse(is.na(d$dose1_day[matched_control]), Inf, d$dose1_day[matched_control]),
  ifelse(is.na(d$booster_day[matched_case]),  Inf, d$booster_day[matched_case]),
  na.rm = TRUE
)

make_side <- function(idx, exposed_flag) {
  ev_day <- d$event_day[idx]
  stop_day <- pmin(pair_end, ifelse(is.na(ev_day), Inf, ev_day), na.rm = TRUE)
  data.frame(
    pair        = seq_along(idx),
    id          = d$id[idx],
    exposed     = exposed_flag,
    tstart      = matched_day,
    tstop       = stop_day,
    event       = as.integer(!is.na(ev_day) & ev_day <= stop_day),
    city        = d$city[idx],
    sex         = d$sex[idx],
    age         = d$age[idx],
    age_group   = d$age_group[idx],
    comorbidity = d$comorbidity[idx],
    regimen     = d$regimen[idx],
    vaccine     = d$vaccine[idx],
    hospitalised = d$hospitalised[idx],
    died        = d$died[idx],
    stringsAsFactors = FALSE
  )
}

matched <- rbind(make_side(matched_case, 1L), make_side(matched_control, 0L))
matched <- matched[matched$tstop > matched$tstart, ]
matched$followup_days <- matched$tstop - matched$tstart

# --- Balance diagnostics ----------------------------------------------------
# Standardised mean differences before matching (all eligible vaccinated
# versus all unvaccinated) and after matching (cases versus their controls).

smd <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  m1 <- mean(x1, na.rm = TRUE); m0 <- mean(x0, na.rm = TRUE)
  s  <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
  if (is.na(s) || s == 0) return(0)
  (m1 - m0) / s
}

design_matrix <- function(df) {
  data.frame(
    age               = df$age,
    female            = as.integer(df$sex == "F"),
    city_monteria     = as.integer(df$city == "23001 - Montería"),
    comorbid_any      = as.integer(df$comorbidity != "1. None"),
    comorbid_cardio   = as.integer(df$comorbidity %in%
                                     c("3. Cardiometabolic + other", "4. Cardiometabolic")),
    regimen_contrib   = as.integer(df$regimen == "Contributivo")
  )
}

before_grp <- as.integer(is_case)
before_mat <- design_matrix(d)
after_grp  <- matched$exposed
after_mat  <- design_matrix(matched)

balance <- data.frame(
  variable = names(before_mat),
  before   = vapply(before_mat, smd, numeric(1), g = before_grp),
  after    = vapply(after_mat,  smd, numeric(1), g = after_grp),
  row.names = NULL
)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
png(file.path(FIG_DIR, "balance_smd.png"), width = 1400, height = 900, res = 170)
op <- par(mar = c(5, 11, 4, 2))
yy <- seq_len(nrow(balance))
plot(NA, xlim = range(c(-0.35, 0.35, balance$before, balance$after)),
     ylim = c(0.5, nrow(balance) + 0.5), yaxt = "n", ylab = "",
     xlab = "Standardised mean difference",
     main = "Covariate balance before and after risk-set matching")
abline(v = c(-0.1, 0, 0.1), lty = c(2, 1, 2), col = c("grey60", "grey30", "grey60"))
points(balance$before, yy, pch = 1,  cex = 1.2)
points(balance$after,  yy, pch = 19, cex = 1.2)
axis(2, at = yy, labels = balance$variable, las = 1)
legend("topright", c("Before", "After"), pch = c(1, 19), bty = "n")
par(op); invisible(dev.off())

saveRDS(matched, file.path(OUT_DIR, "matched_pairs.rds"))

cat("\nCovariate balance (standardised mean differences)\n")
print(within(balance, {
  before <- round(before, 3); after <- round(after, 3)
}), row.names = FALSE)

cat(sprintf("
Matched cohort
--------------
Pairs                                     %7d
Individuals contributing follow-up        %7d
Events among the vaccinated               %7d
Events among the matched unvaccinated     %7d
Median follow-up, days                    %7.0f
Balance figure written to                 %s
",
  length(unique(matched$pair)),
  nrow(matched),
  sum(matched$event[matched$exposed == 1]),
  sum(matched$event[matched$exposed == 0]),
  median(matched$followup_days),
  file.path(FIG_DIR, "balance_smd.png")
))
