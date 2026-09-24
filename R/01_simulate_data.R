# ---------------------------------------------------------------------------
# 01_simulate_data.R
#
# Generates a synthetic dataset that mirrors the structure of the linked
# national records used in the original analysis (immunisation registry,
# surveillance system and vital statistics) WITHOUT reproducing any real
# record. Variable names, categories, date ranges and missingness patterns
# match the original extract, so the rest of the pipeline runs unchanged.
#
# Infections are simulated day by day in calendar time, from a baseline hazard
# with two epidemic waves. Vaccination acts as a TIME-VARYING exposure: an
# individual's hazard is multiplied by (1 - TRUE_VE) only from 14 days after
# the second dose onwards. This is what lets 02_build_cohort.R and
# 05_models.R be validated: a correct time-varying analysis recovers TRUE_VE,
# and analyses that mishandle time zero do not.
#
# Confounding is built in as well: uptake depends on age, comorbidity and
# insurance scheme, which also predict infection risk.
#
# Base R only: no packages are required to generate the data.
#
# Output: data/synthetic/synthetic_cohort.csv
# ---------------------------------------------------------------------------

set.seed(20240813)  # thesis defence date, for luck and reproducibility

# --- Parameters ------------------------------------------------------------

N            <- 60000                      # individuals
STUDY_START  <- as.Date("2021-02-17")      # national vaccination plan start
STUDY_END    <- as.Date("2022-06-30")      # administrative end of follow-up
IMMUNITY_LAG <- 14                         # days from 2nd dose to protection
TRUE_VE      <- 0.55                       # true vaccine effectiveness
OUT_DIR      <- "data/synthetic"
OUT_FILE     <- file.path(OUT_DIR, "synthetic_cohort.csv")

CITIES   <- c("76001 - Cali", "23001 - Montería")
VACCINES <- c("PFIZER", "ASTRAZENECA", "JANSSEN", "MODERNA", "SINOVAC")
REGIMENS <- c("Contributivo", "Subsidiado", "Especial", "Excepcion")

# --- Demographics ----------------------------------------------------------

city <- sample(CITIES, N, replace = TRUE, prob = c(0.72, 0.28))

# Age: gamma mixture giving an adult-heavy distribution with a paediatric tail
age <- round(pmin(pmax(rgamma(N, shape = 4.2, scale = 9.5), 3), 99))

sex <- sample(c("M", "F"), N, replace = TRUE, prob = c(0.47, 0.53))
# The original extract carries three flavours of missing sex; keep them so the
# cleaning code in 02_build_cohort.R is exercised.
junk <- sample(seq_len(N), size = round(0.004 * N))
sex[junk] <- sample(c("INDEFINIDO", "NO DEFINIDO", "NULL"), length(junk), replace = TRUE)

regimen <- sample(REGIMENS, N, replace = TRUE, prob = c(0.52, 0.41, 0.05, 0.02))
regimen[sample(seq_len(N), round(0.02 * N))] <- NA  # unlinked records

# --- Comorbidities (prevalence increases with age) -------------------------

p_age <- plogis((age - 55) / 14)

cac_hta      <- rbinom(N, 1, 0.02  + 0.42 * p_age)
cac_diabetes <- rbinom(N, 1, 0.01  + 0.20 * p_age)
cac_cancer   <- rbinom(N, 1, 0.002 + 0.05 * p_age)
cac_artritis <- rbinom(N, 1, 0.003 + 0.08 * p_age)
cac_peh      <- rbinom(N, 1, 0.004 + 0.03 * p_age)   # chronic kidney disease
cac_vih      <- rbinom(N, 1, 0.004)

any_comorbidity <- as.integer(
  cac_hta | cac_diabetes | cac_cancer | cac_artritis | cac_peh | cac_vih
)

# --- Eligibility date under the National Vaccination Plan ------------------
# Reproduces the priority schedule: the date on which each person's priority
# group was opened. Determined by age and comorbidity at baseline, so it does
# not depend on the individual's own later decisions.

eligibility_date <- as.Date(ifelse(
  age >= 80, "2021-02-17",
  ifelse(age >= 60, "2021-03-08",
  ifelse(age >= 50, "2021-05-22",
  ifelse(any_comorbidity == 1, "2021-05-22",
  ifelse(age >= 40, "2021-06-17",
  ifelse(age >= 12, "2021-07-17",
         "2021-10-29")))))))

# --- Vaccination ------------------------------------------------------------
# Uptake depends on age, comorbidity and insurance scheme: this is the
# confounding structure the propensity score is meant to address.

lp_vax <- -0.35 +
  1.25 * p_age +
  0.45 * any_comorbidity +
  0.35 * (regimen %in% c("Contributivo", "Especial", "Excepcion"))
lp_vax[is.na(lp_vax)] <- -0.35
vaccinated <- rbinom(N, 1, plogis(lp_vax))

# First dose: between eligibility and the end of the roll-out period, sooner
# for the older priority groups.
max_delay <- as.numeric(STUDY_END - eligibility_date) - 40
max_delay[max_delay < 1] <- 1
delay_1d <- round(runif(N) * max_delay * (0.35 + 0.65 * (1 - p_age)))

fec_vac_1d <- as.Date(rep(NA, N))
fec_vac_1d[vaccinated == 1] <- eligibility_date[vaccinated == 1] + delay_1d[vaccinated == 1]

# Second dose: most within the 14-28 day protocol window, some outside it
# (those follow neither treatment strategy).
gap <- round(rnorm(N, mean = 21, sd = 3))
off_protocol <- sample(seq_len(N), round(0.08 * N))
gap[off_protocol] <- round(runif(length(off_protocol), 30, 90))
fec_vac_2d <- fec_vac_1d + gap

# A small share never completed the schedule
dropout <- rbinom(N, 1, 0.05) == 1
fec_vac_2d[dropout] <- as.Date(NA)

nom_vac_1d <- rep(NA_character_, N)
nom_vac_1d[!is.na(fec_vac_1d)] <- sample(
  VACCINES, sum(!is.na(fec_vac_1d)), replace = TRUE,
  prob = c(0.34, 0.24, 0.13, 0.11, 0.18)
)
nom_vac_2d <- ifelse(is.na(fec_vac_2d), NA_character_, nom_vac_1d)

# Booster from late 2021 onwards, for a subset of those who completed the schedule
booster <- !is.na(fec_vac_2d) & rbinom(N, 1, 0.28) == 1
fec_vac_rd <- as.Date(rep(NA, N))
fec_vac_rd[booster] <- fec_vac_2d[booster] + round(runif(sum(booster), 120, 240))
fec_vac_rd[!is.na(fec_vac_rd) & fec_vac_rd > STUDY_END] <- as.Date(NA)
nom_vac_rd <- rep(NA_character_, N)
nom_vac_rd[!is.na(fec_vac_rd)] <- sample(VACCINES, sum(!is.na(fec_vac_rd)), replace = TRUE)

# Date from which protection applies (NA for those who never complete the
# schedule within the protocol window)
complete_schedule <- !is.na(fec_vac_2d) &
  (fec_vac_2d - fec_vac_1d) >= 14 & (fec_vac_2d - fec_vac_1d) <= 28
protection_start <- as.Date(rep(NA, N))
protection_start[complete_schedule] <- fec_vac_2d[complete_schedule] + IMMUNITY_LAG

# --- Infections: daily hazard in calendar time -----------------------------
# Baseline hazard: a low constant plus two Gaussian epidemic waves, roughly
# matching the Colombian third wave (mid-2021) and the Omicron wave (early
# 2022). Individual hazard is multiplied by covariate effects and, from
# protection_start onwards, by (1 - TRUE_VE).

days <- seq(STUDY_START, STUDY_END, by = "day")

wave <- function(d, peak, sd, amp) amp * exp(-0.5 * ((as.numeric(d - peak)) / sd)^2)
baseline_daily <- 0.00018 +
  wave(days, as.Date("2021-06-20"), 32, 0.0022) +
  wave(days, as.Date("2022-01-12"), 20, 0.0038)

indiv_multiplier <- exp(0.15 * (city == CITIES[2]) - 0.20 * p_age)
hr_vaccinated    <- 1 - TRUE_VE

infection_day <- rep(NA_real_, N)          # days since STUDY_START
still_at_risk  <- rep(TRUE, N)
elig_num       <- as.numeric(eligibility_date - STUDY_START)
prot_num       <- as.numeric(protection_start - STUDY_START)

for (t in seq_along(days)) {
  idx <- which(still_at_risk & elig_num <= (t - 1))
  if (!length(idx)) next

  protected <- !is.na(prot_num[idx]) & prot_num[idx] <= (t - 1)
  p <- baseline_daily[t] * indiv_multiplier[idx] * ifelse(protected, hr_vaccinated, 1)
  p[p > 0.9] <- 0.9

  hit <- idx[runif(length(idx)) < p]
  if (length(hit)) {
    infection_day[hit] <- t - 1
    still_at_risk[hit] <- FALSE
  }
}

confirmado    <- as.integer(!is.na(infection_day))
symptom_onset <- STUDY_START + infection_day

# --- Severity, conditional on infection ------------------------------------

was_protected <- as.integer(!is.na(prot_num) & !is.na(infection_day) &
                              prot_num <= infection_day)

p_hosp <- plogis(-4.4 + 3.1 * p_age + 0.55 * any_comorbidity - 0.85 * was_protected)
hospitalised <- ifelse(confirmado == 1, rbinom(N, 1, p_hosp), 0)

p_death <- plogis(-2.0 + 2.2 * p_age + 0.45 * any_comorbidity - 0.70 * was_protected)
died <- ifelse(hospitalised == 1, rbinom(N, 1, p_death), 0)

# --- Assemble ---------------------------------------------------------------
# Column names deliberately match the original registry extract.

synthetic <- data.frame(
  MunicipioAplicacion      = city,
  Sexo                     = sex,
  edad_Cal                 = age,
  Regimen                  = regimen,
  CAC_HTA                  = cac_hta,
  CAC_Diabetes             = cac_diabetes,
  CAC_Cancer               = cac_cancer,
  CAC_Artritis             = cac_artritis,
  CAC_PEH                  = cac_peh,
  CAC_VIH                  = cac_vih,
  nom_vac_1d               = nom_vac_1d,
  fec_vac_1d               = fec_vac_1d,
  nom_vac_2d               = nom_vac_2d,
  fec_vac_2d               = fec_vac_2d,
  nom_vac_rd               = nom_vac_rd,
  fec_vac_rd               = fec_vac_rd,
  Confirmado               = confirmado,
  FechaInicioSintomas      = symptom_onset,
  ServicioMayorComplejidad = hospitalised,
  NDEstadoVital            = died,
  stringsAsFactors         = FALSE
)

# Records lost to linkage: whole rows missing the surveillance fields
lost <- sample(seq_len(N), round(0.015 * N))
synthetic$Confirmado[lost] <- NA
synthetic$FechaInicioSintomas[lost] <- as.Date(NA)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(synthetic, OUT_FILE, row.names = FALSE, na = "")

message(sprintf(
"Wrote %s
  rows                              %7d
  complete primary schedule         %7d (%.1f%%)
  confirmed infections              %7d (%.1f%%)
  true vaccine effectiveness        %7.0f%%   (time-varying, from 2nd dose + %d days)",
  OUT_FILE, nrow(synthetic),
  sum(complete_schedule), 100 * mean(complete_schedule),
  sum(synthetic$Confirmado, na.rm = TRUE),
  100 * mean(synthetic$Confirmado, na.rm = TRUE),
  100 * TRUE_VE, IMMUNITY_LAG
))
