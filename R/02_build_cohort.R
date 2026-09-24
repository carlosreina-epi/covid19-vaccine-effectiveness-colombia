# ---------------------------------------------------------------------------
# 02_build_cohort.R
#
# Builds the analysis cohort from the registry extract.
#
# Design note, and the main methodological change relative to the original
# thesis code:
#
#   In the thesis, unvaccinated individuals were defined as those with no dose
#   recorded at any point during follow-up. That definition conditions group
#   membership on future behaviour, which is a known source of bias: a person
#   who was unvaccinated for ten months and then vaccinated contributes no
#   unexposed person-time, and the unexposed group is restricted to people who
#   never took up the vaccine, who differ systematically from those who did.
#
#   Here exposure is handled as a time-varying variable instead. Every person
#   contributes unexposed person-time from the date their priority group became
#   eligible until they are vaccinated (at which point unexposed follow-up is
#   censored), and those who complete the primary schedule then contribute
#   exposed person-time from 14 days after the second dose. The result is a
#   counting-process dataset with (tstart, tstop] intervals, analysed with Cox
#   models in 05_models.R.
#
#   The original "never vaccinated" definition is retained as a sensitivity
#   analysis so that the two can be compared directly; see the `analysis_set`
#   column.
#
#   Second design decision: the time scale is CALENDAR TIME, expressed as days
#   since the start of the vaccination programme, rather than time since each
#   individual's own origin. Vaccinated person-time necessarily falls later in
#   the calendar than unvaccinated person-time, and the force of infection
#   varied enormously between epidemic waves; using time since individual
#   origin would compare person-time from different phases of the epidemic and
#   inflate the estimated effectiveness. With calendar time as the time scale,
#   each risk set contains only person-time from the same day of the epidemic.
#
# Input:  data/synthetic/synthetic_cohort.csv  (or the real extract, same names)
# Output: data/derived/cohort_person_time.rds
#         data/derived/cohort_baseline.rds
# ---------------------------------------------------------------------------

source("R/00_functions.R")

# --- Parameters ------------------------------------------------------------

INPUT_FILE   <- "data/synthetic/synthetic_cohort.csv"
OUT_DIR      <- "data/derived"
STUDY_START  <- as.Date("2021-02-17")
STUDY_END    <- as.Date("2022-06-30")
IMMUNITY_LAG <- 14    # days after the second dose before exposure starts
MIN_GAP      <- 14    # protocol window between doses, lower bound
MAX_GAP      <- 28    # protocol window between doses, upper bound
MAX_FOLLOWUP <- 360   # days

CITIES <- c("76001 - Cali", "23001 - Montería")

# --- Read and clean --------------------------------------------------------

raw <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, na.strings = c("", "NA"))

date_cols <- c("fec_vac_1d", "fec_vac_2d", "fec_vac_rd", "FechaInicioSintomas")
raw[date_cols] <- lapply(raw[date_cols], as.Date)

d <- data.frame(
  city        = factor(raw$MunicipioAplicacion, levels = CITIES),
  sex         = clean_sex(raw$Sexo),
  age         = raw$edad_Cal,
  regimen     = factor(raw$Regimen),
  dose1_date  = raw$fec_vac_1d,
  dose2_date  = raw$fec_vac_2d,
  booster_date = raw$fec_vac_rd,
  vaccine     = raw$nom_vac_2d,
  infected    = raw$Confirmado,
  onset_date  = raw$FechaInicioSintomas,
  hospitalised = raw$ServicioMayorComplejidad,
  died        = raw$NDEstadoVital,
  stringsAsFactors = FALSE
)

d$any_comorbidity <- as.integer(
  raw$CAC_HTA == 1 | raw$CAC_Diabetes == 1 | raw$CAC_VIH == 1 |
  raw$CAC_PEH == 1 | raw$CAC_Cancer == 1 | raw$CAC_Artritis == 1
)
d$comorbidity <- comorbidity_group(
  raw$CAC_HTA, raw$CAC_Diabetes, raw$CAC_VIH,
  raw$CAC_PEH, raw$CAC_Cancer, raw$CAC_Artritis
)
d$age_group <- age_group(d$age)
d$id <- seq_len(nrow(d))

# --- Eligibility -----------------------------------------------------------

d$eligibility_date <- pnv_eligibility_date(d$age, d$any_comorbidity)

n_start <- nrow(d)
d <- d[!is.na(d$city), ]
n_city <- nrow(d)
d <- d[d$eligibility_date >= STUDY_START & d$eligibility_date <= STUDY_END, ]
n_elig <- nrow(d)

# Records with no surveillance linkage cannot contribute outcome information
d <- d[!is.na(d$infected), ]
n_linked <- nrow(d)

# --- Treatment strategies --------------------------------------------------
# Strategy 1: complete the primary schedule, two doses MIN_GAP to MAX_GAP days
#             apart. Strategy 2: remain unvaccinated.
# Individuals whose two doses fall outside the protocol window follow neither
# strategy: their exposed person-time is excluded, but their unexposed
# person-time before the first dose still counts.

d$dose_gap <- as.numeric(d$dose2_date - d$dose1_date)
d$complete_schedule <- !is.na(d$dose_gap) &
  d$dose_gap >= MIN_GAP & d$dose_gap <= MAX_GAP
d$exposure_start <- as.Date(ifelse(
  d$complete_schedule, as.character(d$dose2_date + IMMUNITY_LAG), NA
))

# Event date: symptom onset for confirmed infections
d$event_date <- d$onset_date
d$event_date[d$infected != 1] <- NA

# --- Unexposed person-time -------------------------------------------------
# From eligibility until the earliest of: first dose, event, administrative
# end, or MAX_FOLLOWUP days.

unexposed_end <- pmin(
  ifelse(is.na(d$dose1_date), Inf, as.numeric(d$dose1_date)),
  ifelse(is.na(d$event_date), Inf, as.numeric(d$event_date)),
  as.numeric(STUDY_END),
  as.numeric(d$eligibility_date) + MAX_FOLLOWUP,
  na.rm = TRUE
)

unexposed <- data.frame(
  id      = d$id,
  exposed = 0L,
  tstart  = as.numeric(d$eligibility_date - STUDY_START),
  tstop   = unexposed_end - as.numeric(STUDY_START),
  event   = as.integer(
    !is.na(d$event_date) & as.numeric(d$event_date) == unexposed_end
  )
)

# --- Exposed person-time ---------------------------------------------------
# From 14 days after the second dose until the earliest of: booster, event,
# administrative end, or MAX_FOLLOWUP days.

has_exposure <- d$complete_schedule & !is.na(d$exposure_start) &
  d$exposure_start <= STUDY_END

de <- d[has_exposure, ]

exposed_end <- pmin(
  ifelse(is.na(de$booster_date), Inf, as.numeric(de$booster_date)),
  ifelse(is.na(de$event_date), Inf, as.numeric(de$event_date)),
  as.numeric(STUDY_END),
  as.numeric(de$exposure_start) + MAX_FOLLOWUP,
  na.rm = TRUE
)

exposed <- data.frame(
  id      = de$id,
  exposed = 1L,
  tstart  = as.numeric(de$exposure_start - STUDY_START),
  tstop   = exposed_end - as.numeric(STUDY_START),
  event   = as.integer(
    !is.na(de$event_date) & as.numeric(de$event_date) == exposed_end
  )
)

person_time <- rbind(unexposed, exposed)

# Keep only intervals of positive length. An interval can be empty when an
# individual was infected before their exposure was due to start, in which
# case the event already belongs to their unexposed interval.
person_time <- person_time[person_time$tstop > person_time$tstart, ]
person_time$followup_days <- person_time$tstop - person_time$tstart

# Attach baseline covariates
baseline_cols <- c("id", "city", "sex", "age", "age_group", "regimen",
                   "comorbidity", "any_comorbidity", "vaccine",
                   "complete_schedule", "eligibility_date",
                   "hospitalised", "died")
person_time <- merge(person_time, d[, baseline_cols], by = "id", all.x = TRUE)

# --- Sensitivity set: the original "never vaccinated" definition -----------

never_vaccinated <- is.na(d$dose1_date) & is.na(d$dose2_date)
d$analysis_set <- ifelse(
  d$complete_schedule, "exposed",
  ifelse(never_vaccinated, "never vaccinated", "neither")
)

# --- Save and report -------------------------------------------------------

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(person_time, file.path(OUT_DIR, "cohort_person_time.rds"))
saveRDS(d,           file.path(OUT_DIR, "cohort_baseline.rds"))

cat(sprintf(
"Cohort construction
-------------------
Records read                              %7d
  with a study city                       %7d
  eligible within the study period        %7d
  linked to surveillance data             %7d

Person-time intervals                     %7d
  unexposed                               %7d
  exposed (complete primary schedule)     %7d

Events
  during unexposed time                   %7d
  during exposed time                     %7d

Median follow-up (days)
  unexposed                               %7.0f
  exposed                                 %7.0f

Sensitivity set (original definition)
  exposed                                 %7d
  never vaccinated                        %7d
  neither strategy                        %7d
",
n_start, n_city, n_elig, n_linked,
nrow(person_time),
sum(person_time$exposed == 0), sum(person_time$exposed == 1),
sum(person_time$event[person_time$exposed == 0]),
sum(person_time$event[person_time$exposed == 1]),
median(person_time$followup_days[person_time$exposed == 0]),
median(person_time$followup_days[person_time$exposed == 1]),
sum(d$analysis_set == "exposed"),
sum(d$analysis_set == "never vaccinated"),
sum(d$analysis_set == "neither")
))
