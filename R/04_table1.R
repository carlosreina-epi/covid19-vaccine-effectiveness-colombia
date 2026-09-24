# ---------------------------------------------------------------------------
# 04_table1.R
#
# Baseline characteristics of the matched cohort, by treatment strategy.
#
# Written with base R so that the pipeline has no dependency beyond `survival`.
# Output is both a data frame and a Markdown table, so the report in
# analysis/report.qmd and the README can reuse it without rerunning anything.
#
# Input:  data/derived/matched_pairs.rds
# Output: data/derived/table1.rds
#         analysis/table1.md
# ---------------------------------------------------------------------------

source("R/00_functions.R")

OUT_DIR <- "data/derived"
ANA_DIR <- "analysis"

m <- readRDS(file.path(OUT_DIR, "matched_pairs.rds"))

vaccinated   <- m[m$exposed == 1, ]
unvaccinated <- m[m$exposed == 0, ]

# --- Small helpers ----------------------------------------------------------

fmt_median <- function(x) {
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.0f (%.0f-%.0f)", q[2], q[1], q[3])
}

fmt_n_pct <- function(n, total) sprintf("%d (%.1f%%)", n, 100 * n / total)

row_continuous <- function(label, v1, v0) {
  data.frame(characteristic = label,
             vaccinated = fmt_median(v1),
             unvaccinated = fmt_median(v0),
             stringsAsFactors = FALSE)
}

row_categorical <- function(label, f1, f0) {
  lv <- union(levels(addNA(factor(f1), ifany = TRUE)),
              levels(addNA(factor(f0), ifany = TRUE)))
  lv <- lv[!is.na(lv)]
  header <- data.frame(characteristic = label, vaccinated = "",
                       unvaccinated = "", stringsAsFactors = FALSE)
  body <- do.call(rbind, lapply(lv, function(l) {
    data.frame(
      characteristic = paste0("  ", l),
      vaccinated   = fmt_n_pct(sum(f1 == l, na.rm = TRUE), length(f1)),
      unvaccinated = fmt_n_pct(sum(f0 == l, na.rm = TRUE), length(f0)),
      stringsAsFactors = FALSE
    )
  }))
  missing1 <- sum(is.na(f1)); missing0 <- sum(is.na(f0))
  if (missing1 + missing0 > 0) {
    body <- rbind(body, data.frame(
      characteristic = "  Missing",
      vaccinated   = fmt_n_pct(missing1, length(f1)),
      unvaccinated = fmt_n_pct(missing0, length(f0)),
      stringsAsFactors = FALSE))
  }
  rbind(header, body)
}

# --- Build the table --------------------------------------------------------

table1 <- rbind(
  data.frame(characteristic = "Individuals",
             vaccinated   = format(nrow(vaccinated), big.mark = ","),
             unvaccinated = format(nrow(unvaccinated), big.mark = ","),
             stringsAsFactors = FALSE),
  row_continuous("Age, years, median (IQR)", vaccinated$age, unvaccinated$age),
  row_categorical("Age group", vaccinated$age_group, unvaccinated$age_group),
  row_categorical("Sex", vaccinated$sex, unvaccinated$sex),
  row_categorical("City", vaccinated$city, unvaccinated$city),
  row_categorical("Comorbidity group", vaccinated$comorbidity, unvaccinated$comorbidity),
  row_categorical("Insurance scheme", vaccinated$regimen, unvaccinated$regimen),
  row_continuous("Follow-up, days, median (IQR)",
                 vaccinated$followup_days, unvaccinated$followup_days),
  data.frame(characteristic = "Confirmed infections",
             vaccinated   = fmt_n_pct(sum(vaccinated$event), nrow(vaccinated)),
             unvaccinated = fmt_n_pct(sum(unvaccinated$event), nrow(unvaccinated)),
             stringsAsFactors = FALSE)
)

names(table1) <- c("Characteristic", "Complete primary schedule", "Unvaccinated")

# --- Write ------------------------------------------------------------------

dir.create(ANA_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(table1, file.path(OUT_DIR, "table1.rds"))

md <- c(
  "# Table 1. Baseline characteristics of the matched cohort",
  "",
  paste0("| ", paste(names(table1), collapse = " | "), " |"),
  paste0("|", paste(rep("---", ncol(table1)), collapse = "|"), "|"),
  apply(table1, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")),
  "",
  "Matching was exact on city, sex, age group, comorbidity group and insurance",
  "scheme within risk sets, so the two columns are identical by construction on",
  "those variables; they are shown to document the cohort, not to test balance.",
  "Balance is reported as standardised mean differences in `figures/balance_smd.png`."
)
writeLines(md, file.path(ANA_DIR, "table1.md"))

# Print with cat rather than print(), which escapes non-ASCII characters and
# would show "Monter\303\255a" instead of "Montería" in the console.
w1 <- max(nchar(table1[[1]]))
w2 <- max(nchar(names(table1)[2]), max(nchar(table1[[2]])))
cat(sprintf("%-*s  %-*s  %s\n", w1, names(table1)[1], w2, names(table1)[2], names(table1)[3]))
cat(strrep("-", w1 + w2 + 20), "\n")
for (i in seq_len(nrow(table1))) {
  cat(sprintf("%-*s  %-*s  %s\n", w1, table1[i, 1], w2, table1[i, 2], table1[i, 3]))
}

cat(sprintf("\nWritten to %s and %s\n",
            file.path(OUT_DIR, "table1.rds"), file.path(ANA_DIR, "table1.md")))
