# ---------------------------------------------------------------------------
# 00_functions.R
#
# Shared helper functions. Sourced by every other script.
# No side effects: this file only defines functions.
# ---------------------------------------------------------------------------

#' Eligibility date under Colombia's National Vaccination Plan
#'
#' Returns the calendar date on which each individual's priority group was
#' opened for vaccination. This is the date used as time zero for person-time
#' that has not yet been exposed, which is what makes the emulation possible:
#' it is determined by age and documented comorbidity at baseline and therefore
#' does not depend on the individual's own later decisions.
#'
#' @param age Numeric vector of ages in years.
#' @param any_comorbidity Integer vector, 1 if any documented comorbidity.
#' @return A Date vector.
pnv_eligibility_date <- function(age, any_comorbidity) {
  stopifnot(length(age) == length(any_comorbidity))
  out <- rep(NA_character_, length(age))
  out[age >= 12]                      <- "2021-07-17"
  out[age <  12]                      <- "2021-10-29"
  out[age >= 40]                      <- "2021-06-17"
  out[any_comorbidity == 1 & age < 50] <- "2021-05-22"
  out[age >= 50]                      <- "2021-05-22"
  out[age >= 60]                      <- "2021-03-08"
  out[age >= 80]                      <- "2021-02-17"
  as.Date(out)
}

#' Recode the three flavours of missing sex found in the registry extract
clean_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("INDEFINIDO", "NO DEFINIDO", "NULL", "")] <- NA
  factor(x, levels = c("F", "M"))
}

#' Collapse the six comorbidity flags into the four-level grouping used in the
#' original analysis.
comorbidity_group <- function(hta, diabetes, vih, peh, cancer, artritis) {
  cardiometabolic <- (hta == 1 | diabetes == 1)
  other           <- (vih == 1 | peh == 1 | cancer == 1 | artritis == 1)
  out <- rep(NA_character_, length(hta))
  out[!cardiometabolic & !other] <- "1. None"
  out[!cardiometabolic &  other] <- "2. Other"
  out[ cardiometabolic &  other] <- "3. Cardiometabolic + other"
  out[ cardiometabolic & !other] <- "4. Cardiometabolic"
  factor(out, levels = c("1. None", "2. Other",
                         "3. Cardiometabolic + other", "4. Cardiometabolic"))
}

#' Age groups used throughout the analysis
age_group <- function(age) {
  cut(age,
      breaks = c(0, 12, 26, 40, 50, 60, Inf),
      labels = c("3-11", "12-25", "26-39", "40-49", "50-59", ">=60"),
      right  = FALSE)
}

#' Vaccine effectiveness with confidence interval from a fitted Cox model
#'
#' VE = 1 - HR. The standard error is always taken from the model passed in,
#' which avoids the easy mistake of reusing a variance matrix from a previous,
#' differently specified model.
#'
#' @param fit A coxph object.
#' @param term Name of the coefficient to convert, e.g. "exposed".
#' @param conf_level Confidence level, default 0.95.
#' @return A one-row data frame with the estimate and confidence limits.
ve_from_cox <- function(fit, term, conf_level = 0.95) {
  coefs <- coef(fit)
  if (!term %in% names(coefs)) {
    stop("Term '", term, "' not found. Available: ",
         paste(names(coefs), collapse = ", "))
  }
  beta <- coefs[[term]]
  se   <- sqrt(diag(vcov(fit)))[[term]]
  z    <- qnorm(1 - (1 - conf_level) / 2)

  data.frame(
    term     = term,
    hr       = exp(beta),
    hr_low   = exp(beta - z * se),
    hr_high  = exp(beta + z * se),
    ve       = 1 - exp(beta),
    ve_low   = 1 - exp(beta + z * se),   # note the deliberate swap:
    ve_high  = 1 - exp(beta - z * se),   # VE is decreasing in the HR
    row.names = NULL
  )
}

#' Pretty-print a vaccine effectiveness row as a percentage
format_ve <- function(x, digits = 1) {
  sprintf("%.*f%% (95%% CI %.*f to %.*f)",
          digits, 100 * x$ve, digits, 100 * x$ve_low, digits, 100 * x$ve_high)
}
