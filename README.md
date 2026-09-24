# COVID-19 vaccine effectiveness in two Colombian cities: a target trial emulation

Reproducible analysis pipeline for a population-based cohort study estimating the
effectiveness of Colombia's National Vaccination Plan against SARS-CoV-2 infection,
hospitalisation and death in Cali and Montería (February 2021 – June 2022).

This repository contains the full analysis code and a synthetic data generator, so
that every step can be run end to end without access to the original records.

## The causal question

Among residents of Cali and Montería who became eligible for vaccination under the
National Vaccination Plan, what would the risk of SARS-CoV-2 infection,
hospitalisation and death have been had everyone completed the primary vaccination
schedule, compared with had no one been vaccinated?

## Target trial protocol

| Component | Specification |
|---|---|
| **Eligibility** | Residents of Cali or Montería registered in the national vaccination and surveillance registries, alive and eligible for vaccination under the National Vaccination Plan between 17 Feb 2021 and 30 Jun 2022 |
| **Treatment strategies** | (1) Complete the primary schedule: two doses 14–28 days apart. (2) Remain unvaccinated |
| **Assignment** | Observational; confounding addressed by propensity score matching and covariate adjustment |
| **Time zero** | Vaccinated: date of second dose + 14 days. Unvaccinated: date on which the individual's priority group became eligible under the National Vaccination Plan, derived from age and documented comorbidities |
| **Outcomes** | Laboratory-confirmed SARS-CoV-2 infection; hospitalisation; death |
| **Follow-up** | From time zero to the outcome, 360 days, or 30 June 2022, whichever came first |
| **Causal contrast** | Per-protocol effect of completing the primary schedule |
| **Analysis** | Cox proportional hazards; vaccine effectiveness estimated as 1 − HR, overall, by follow-up window, and by age, sex, comorbidity, city and vaccine product |

Defining time zero for unvaccinated individuals from the national priority schedule
is the design feature that makes this emulation possible: it gives every person an
eligibility date that does not depend on their own later behaviour.

## Data

The analysis uses individual-level records from the national immunisation registry
(PAIweb), the national surveillance system (SIVIGILA) and vital statistics, linked
by national identifier under a data use agreement. **These data are not, and will not
be, included in this repository.**

`R/01_simulate_data.R` generates a synthetic dataset with the same variables,
distributions and missingness patterns, so that the pipeline is fully runnable.
Results obtained from synthetic data are illustrative and are not study findings.

## Repository structure

R/00_functions.R Shared helper functions
R/01_simulate_data.R Synthetic data generator
R/02_build_cohort.R Eligibility, time zero, treatment strategies
R/03_match.R Propensity score matching and balance diagnostics
R/04_table1.R Baseline characteristics
R/05_models.R Cox models, effectiveness by window and subgroup
R/06_diagnostics.R Proportional hazards, sensitivity analyses
analysis/report.qmd Reproducible report


## How to run

```r
install.packages("renv")
renv::restore()
source("R/01_simulate_data.R")
quarto::quarto_render("analysis/report.qmd")
```

## Known limitations

Stated explicitly rather than buried: unvaccinated individuals were identified from
the absence of any dose record, which conditions group membership on later behaviour;
a time-varying exposure definition is the preferred alternative and is implemented as
a sensitivity analysis. Outcome ascertainment depends on testing, which varied over
time and between cities. The propensity score model is deliberately parsimonious and
balance is reported for all baseline covariates.

## Citation

Reina Bolaños CA. *Evaluación de efectividad del Plan Nacional de Vacunación contra
el COVID-19 en dos ciudades colombianas durante el periodo de emergencia sanitaria.*
PhD thesis, Facultad Nacional de Salud Pública, Universidad de Antioquia, 2024.

## Author

Carlos Alberto Reina Bolaños, PhD — [ORCID 0000-0002-6367-9239](https://orcid.org/0000-0002-6367-9239)

## License

MIT (see `LICENSE`).