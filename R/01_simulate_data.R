# ---------------------------------------------------------------------------
# 01_simulate_data.R
#
# Generates a synthetic dataset that mirrors the structure of the linked
# national records used in the original analysis (immunisation registry,
# surveillance system and vital statistics) WITHOUT reproducing any real
# record. Variable names, categories, date ranges and missingness patterns
# match the original extract, so the rest of the pipeline runs unchanged.
#
# A true vaccine effect is built in (see TRUE_VE below), together with
# confounding by age, comorbidity and insurance scheme, so that the analysis
# scripts have something real to recover.
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
TRUE_VE      <- 0.55                       # true vaccine effectiveness
OUT_DIR      <- "data/synthetic"
OUT_FILE     <- file.path(OUT_DIR, "synthetic_cohort.csv")

CITIES   <- c("76001 - Cali", "23001 - Montería")
VACCINES <- c("PFIZER", "ASTRAZENECA", "JANSSEN", "MODERNA", "SINOVAC")
REGIMENS <- c("Contributivo", "Subsidiado", "Especial", "Excepcion")

# --- Demographics ----------------------------------------------------------

city <- sample(CITIES, N, replace = TRUE, prob = c(0.72, 0.28))

# Age: mixture giving a realistic adult-heavy distribution with a paediatric tail
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

cac_hta      <- rbinom(N, 1, 0.02 + 0.42 * p_age)
cac_diabetes <- rbinom(N, 1, 0.01 + 0.20 * p_age)
cac_cancer   <- rbinom(N, 1, 0.002 + 0.05 * p_age)
cac_artritis <- rbinom(N, 1, 0.003 + 0.08 * p_age)
cac_peh      <- rbinom(N, 1, 0.004 + 0.03 * p_age)   # chronic kidney disease
cac_vih      <- rbinom(N, 1, 0.004)

any_comorbidity <- as.integer(
  cac_hta | cac_diabetes | cac_cancer | cac_artritis | cac_peh | cac_vih
)

# --- Eligibility date under the National Vaccination Plan ------------------
# Reproduces the priority schedule: the date on which each person's priority
# group was opened. This is the time zero used for unvaccinated individuals.

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

# First dose: uniformly distributed between eligibility and the end of the
# roll-out period, with the delay shorter for older priority groups.
max_delay <- as.numeric(STUDY_END - eligibility_date) - 40
max_delay[max_delay < 1] <- 1
delay_1d <- round(runif(N) * max_delay * (0.35 + 0.65 * (1 - p_age)))

fec_vac_1d <- as.Date(rep(NA, N))
fec_vac_1d[vaccinated == 1] <- eligibility_date[vaccinated == 1] + delay_1d[vaccinated == 1]

# Second dose: most within the 14-28 day protocol window, some outside it
# (those become the "otro" category and are excluded from both strategies).
gap <- sample(c(sample(14:28, 1), 21), N, replace = TRUE)
gap <- round(rnorm(N, mean = 21, sd = 3))
gap[sample(seq_len(N), round(0.08 * N))] <- round(runif(round(0.08 * N), 30, 90))
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
nom_vac_2d[nom_vac_1d == "JANSSEN" & !is.na(nom_vac_1d)] <- NA_character_  # single dose

# Booster from late 2021 onwards, for a subset of those who completed the schedule
booster <- !is.na(fec_vac_2d) & rbinom(N, 1, 0.28) == 1
fec_vac_rd <- as.Date(rep(NA, N))
fec_vac_rd[booster] <- fec_vac_2d[booster] + round(runif(sum(booster), 120, 240))
fec_vac_rd[!is.na(fec_vac_rd) & fec_vac_rd > STUDY_END] <- as.Date(NA)
nom_vac_rd <- rep(NA_character_, N)
nom_vac_rd[!is.na(fec_vac_rd)] <- sample(VACCINES, sum(!is.na(fec_vac_rd)), replace = TRUE)

# --- Time zero --------------------------------------------------------------
# Vaccinated: second dose + 14 days. Unvaccinated: eligibility date.

complete_schedule <- !is.na(fec_vac_2d) & (fec_vac_2d - fec_vac_1d) >= 14 &
  (fec_vac_2d - fec_vac_1d) <= 28

time_zero <- as.Date(ifelse(
  complete_schedule, as.character(fec_vac_2d + 14), as.character(eligibility_date)
))

# --- Outcomes ---------------------------------------------------------------
# Time to laboratory-confirmed infection, simulated from an exponential hazard
# that depends on vaccination status, age and city.

baseline_hazard <- 0.0016
hr_vax <- 1 - TRUE_VE

lambda <- baseline_hazard *
  ifelse(complete_schedule, hr_vax, 1) *
  exp(0.15 * (city == CITIES[2])) *
  exp(-0.20 * p_age)                       # older people shielded more

time_to_event <- rexp(N, rate = lambda)
time_to_admin_end <- as.numeric(STUDY_END - time_zero)

confirmado <- as.integer(time_to_event <= pmin(time_to_admin_end, 360))
symptom_onset <- as.Date(rep(NA, N))
symptom_onset[confirmado == 1] <- time_zero[confirmado == 1] +
  round(time_to_event[confirmado == 1])

# Hospitalisation and death, conditional on infection and modified by age,
# comorbidity and vaccination.
p_hosp <- plogis(-4.4 + 3.1 * p_age + 0.55 * any_comorbidity - 0.85 * complete_schedule)
hospitalised <- ifelse(confirmado == 1, rbinom(N, 1, p_hosp), 0)

p_death <- plogis(-2.0 + 2.2 * p_age + 0.45 * any_comorbidity - 0.70 * complete_schedule)
died <- ifelse(hospitalised == 1, rbinom(N, 1, p_death), 0)

# --- Assemble ---------------------------------------------------------------
# Column names deliberately match the original registry extract.

synthetic <- data.frame(
  MunicipioAplicacion     = city,
  Sexo                    = sex,
  edad_Cal                = age,
  Regimen                 = regimen,
  CAC_HTA                 = cac_hta,
  CAC_Diabetes            = cac_diabetes,
  CAC_Cancer              = cac_cancer,
  CAC_Artritis            = cac_artritis,
  CAC_PEH                 = cac_peh,
  CAC_VIH                 = cac_vih,
  nom_vac_1d              = nom_vac_1d,
  fec_vac_1d              = fec_vac_1d,
  nom_vac_2d              = nom_vac_2d,
  fec_vac_2d              = fec_vac_2d,
  nom_vac_rd              = nom_vac_rd,
  fec_vac_rd              = fec_vac_rd,
  Confirmado              = confirmado,
  FechaInicioSintomas     = symptom_onset,
  ServicioMayorComplejidad = hospitalised,
  NDEstadoVital           = died,
  stringsAsFactors        = FALSE
)

# Records lost to linkage: whole rows missing the surveillance fields
lost <- sample(seq_len(N), round(0.015 * N))
synthetic$Confirmado[lost] <- NA
synthetic$FechaInicioSintomas[lost] <- as.Date(NA)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(synthetic, OUT_FILE, row.names = FALSE, na = "")

message(sprintf(
  "Wrote %s: %d rows, %d columns.\n  Complete primary schedule: %d (%.1f%%)\n  Confirmed infections: %d (%.1f%%)\n  True vaccine effectiveness built in: %.0f%%",
  OUT_FILE, nrow(synthetic), ncol(synthetic),
  sum(complete_schedule), 100 * mean(complete_schedule),
  sum(synthetic$Confirmado, na.rm = TRUE),
  100 * mean(synthetic$Confirmado, na.rm = TRUE),
  100 * TRUE_VE
))
