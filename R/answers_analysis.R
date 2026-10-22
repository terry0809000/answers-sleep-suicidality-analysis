# ANSWERS preregistration: prepare data and run the registered analyses
#
# This script follows the analysis set-ups in Question 1 of the
# preregistration. N1 is the primary network for H1a-H1e and H2; N2-N5 are
# sensitivity networks; P1 is a separate classifier for H3; and D1 describes
# the rarer suicide-related outcomes. A sensitivity network can therefore help
# us understand N1 without being allowed to replace its decision.
#


# 01 - Required packages and frozen analysis settings ----------------------

# Question 15 names the software used for the analysis.
registered_packages <- c(
  "readr", "dplyr", "tidyr", "janitor", "stringr", "ggplot2",
  "psych", "lavaan", "qgraph", "bootnet", "networktools", "glmnet",
  "rsample", "recipes", "pROC", "PRROC", "yardstick", "Matrix",
  "writexl", "renv"
)

# 02 - Frozen analysis constants and data-preparation functions --------------

# Question 15 says to use set.seed() and record the seed. We fix 20260817 before seeing any results
# and record that this number was chosen when the script was written.
ANALYSIS_SEED <- 20260817L

# These are the exact resampling totals specified in Questions 12, 14, 15 and
# 23: 1,000 edge bootstraps, 1,000 ordinary case-drop bootstraps, 1,000 custom
# repetitions at each removal fraction, and 1,000 P1 participant bootstraps.
N_EDGE_BOOT <- 1000L
N_CASE_BOOT <- 1000L
N_CUSTOM_BOOT <- 1000L
CUSTOM_DROP_FRACTIONS <- c(0.10, 0.25, 0.50, 0.75)
N_P1_BOOT <- 1000L

# One core keeps the random draws in the same order across machines. This can
# take longer, but it does not change the quantity being estimated.
N_CORES <- 1L
stopifnot(is.finite(ANALYSIS_SEED), ANALYSIS_SEED > 0L)

# Question 14 requests percentile edge intervals but does not state their
# level. The script uses the conventional 95% level and discloses this choice
# in the results workbook instead of treating it as prespecified.
EDGE_CI_LEVEL <- 0.95


# Each resampling stage receives a reproducible offset from the same seed. This
# keeps, for example, P1 folds unchanged if a network needs a different number
# of retry attempts on another computer.
stage_seed <- function(offset, base_seed = ANALYSIS_SEED) {
  value <- as.integer((as.double(base_seed) + offset) %% 2147483647)
  if (value == 0L) value <- 1L
  value
}

# Question 8 specifies the full answers-dataset-0.2.0.csv file.
# The nine-column harmonized file cannot be used because it omits the sleep,
# mental-health, interpersonal, impulsivity and suicidality measures needed by
# the registered analyses.
read_and_check_raw_data <- function(input_file) {
  raw <- read.csv(
    input_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # These checks reproduce the audited file structure stated in Question 8:
  # 971 rows, 294 columns, unique non-missing IDs, no repeated ID-visit pair
  # and visit 1 throughout. A mismatch stops the analysis before results.
  stopifnot(
    nrow(raw) == 971L,
    ncol(raw) == 294L,
    !anyNA(raw$id),
    anyDuplicated(raw$id) == 0L,
    anyDuplicated(raw[c("id", "visit")]) == 0L,
    all(raw$visit == 1)
  )
  raw
}

# Every constituent item is required for each reconstructed score. Keeping
# na.rm = FALSE means that one missing item makes the total missing instead of
# silently treating it as zero or prorating the scale. This follows the
# complete-item and no-imputation rules in Questions 8 and 9.
sum_items <- function(data, names) {
  rowSums(data[names], na.rm = FALSE)
}

# Question 8 defines cesd_no_sleep from items 1-10 and 12-20. Item 11 ("my
# sleep was restless") is deliberately excluded so that the depression node
# does not partly repeat sleep content already represented in the network.
# Items 4, 8, 12 and 16 are reverse-scored as 3 minus the response. Item 19 is
# retained, all 19 items are required and no clinical cutoff is applied.
score_registered_measures <- function(raw) {
  cesd_names <- paste0("cesd_", c(1:10, 12:20))
  cesd_scored <- raw[cesd_names]
  for (j in c(4, 8, 12, 16)) {
    cesd_scored[[paste0("cesd_", j)]] <- 3 - cesd_scored[[paste0("cesd_", j)]]
  }
  cesd_no_sleep <- rowSums(cesd_scored, na.rm = FALSE)

  # The remaining totals follow the exact scoring identities in Question 8.
  # ISI uses items 1-7, DDNSI uses items 1-5, GAD-7 uses items 1-7 and perceived
  # burdensomeness uses INQ items 1-6.
  isi_total <- sum_items(raw, paste0("isi_", 1:7))
  ddnsi_total <- sum_items(raw, paste0("ddnsi_", 1:5))
  gad_total <- sum_items(raw, paste0("gad_", 1:7))
  inq_perceivedburden <- sum_items(raw, paste0("inq_", 1:6))
  # Thwarted belongingness keeps INQ items 9, 11 and 12 in their original
  # direction and reverse-scores items 7, 8, 10, 13, 14 and 15 as 8 - response.
  inq_thwartedbelong <-
    sum_items(raw, paste0("inq_", c(9, 11, 12))) +
    rowSums(8 - raw[paste0("inq_", c(7, 8, 10, 13, 14, 15))], na.rm = FALSE)
  # The released UPPS scoring supports 5 - response for items 6, 8, 13 and 15.
  # This resolves the dictionary inconsistency documented in Question 8.
  upps_negativeurgency <- rowSums(
    5 - raw[paste0("upps_", c(6, 8, 13, 15))],
    na.rm = FALSE
  )
  # N2 uses the registered PSQI subjective-quality component, scored 0-3.
  psqi_c1 <- raw$psqi_9 - 1

  # The expected ranges stated in Question 8 are CES-D 0-57, ISI 0-28, DDNSI
  # 0-37, GAD-7 0-21, perceived burdensomeness 6-42, thwarted belongingness
  # 9-63, negative urgency 4-16 and PSQI component 1 0-3. The node audit below
  # records observed minima and maxima before estimation so these ranges and
  # coding directions can be checked in the results workbook.

  # Social jet lag is converted from minutes to hours without taking its
  # absolute value. The form deliberately keeps the sign because earlier and
  # later free-day midsleep represent different circadian shifts. A negative
  # edge is therefore directional rather than evidence of protection
  # (Question 8).
  socialjetlag_h <- raw$socialjetlag / 60

  # Question 8 distinguishes the number of endorsed C-SSRS ideation items from
  # the highest endorsed category. si_3mo_count is the count, whereas the
  # network endpoint si_3mo_maxseverity is the highest numbered endorsed item
  # from 1-5, or 0 if none is endorsed. This choice also addresses the
  # 61 non-cumulative response patterns found among the 885 primary participants,
  # as documented in Question 29.
  cssrs_ideation <- paste0("cssrs_3mo_", 1:5)
  si_3mo_count <- sum_items(raw, cssrs_ideation)
  si_3mo_maxseverity <- do.call(
    pmax,
    c(
      list(rep(0L, nrow(raw))),
      lapply(1:5, function(j) j * raw[[paste0("cssrs_3mo_", j)]])
    )
  )
  si_3mo_yn <- ifelse(
    is.na(si_3mo_maxseverity),
    NA_integer_,
    as.integer(si_3mo_maxseverity > 0)
  )

  list(
    cesd_no_sleep = cesd_no_sleep,
    isi_total = isi_total,
    ddnsi_total = ddnsi_total,
    gad_total = gad_total,
    inq_perceivedburden = inq_perceivedburden,
    inq_thwartedbelong = inq_thwartedbelong,
    upps_negativeurgency = upps_negativeurgency,
    psqi_c1 = psqi_c1,
    socialjetlag_h = socialjetlag_h,
    si_3mo_count = si_3mo_count,
    si_3mo_maxseverity = si_3mo_maxseverity,
    si_3mo_yn = si_3mo_yn
  )
}

# si_3mo_yn records whether any item 1-5 was endorsed. P1 and N4 use it only
# where the form specifically asks for binary any-ideation status.

# Each reconstructed score is checked against the released summary field
# before modelling. Matching both values and missingness helps detect an item,
# reversal or coding mistake without using network results to revise a score.
same_with_na <- function(x, y) {
  same_missingness <- is.na(x) == is.na(y)
  same_observed <- x[!is.na(x) & !is.na(y)] == y[!is.na(x) & !is.na(y)]
  all(same_missingness) && all(same_observed)
}
validate_registered_scores <- function(raw, scores) {
  stopifnot(
    same_with_na(
      scores$cesd_no_sleep + raw$cesd_11,
      raw$cesd_total
    ),
    same_with_na(scores$isi_total, raw$isi_total),
    same_with_na(scores$ddnsi_total, raw$ddnsi_total),
    same_with_na(scores$gad_total, raw$gad_total),
    same_with_na(
      scores$inq_perceivedburden,
      raw$inq_perceivedburden
    ),
    same_with_na(
      scores$inq_thwartedbelong,
      raw$inq_thwartedbelong
    ),
    same_with_na(
      scores$upps_negativeurgency,
      raw$upps_negativeurgency
    ),
    same_with_na(scores$psqi_c1, raw$psqi_c1),
    same_with_na(scores$si_3mo_count, raw$si_3mo_severity),
    same_with_na(
      scores$si_3mo_yn,
      as.integer(raw$si_3mo_yn == "Yes")
    )
  )
  invisible(TRUE)
}

# The released socialjetlag field remains the analysis variable. The msf-msw
# calculation is only a quality-control comparison because those components
# are rounded; its difference is reported but never used to exclude a person
# or replace the registered node.
audit_social_jetlag <- function(raw) {
  socialjetlag_qc_min <- 1440 * (raw$msf - raw$msw)
  maximum_difference <- max(
    abs(raw$socialjetlag - socialjetlag_qc_min),
    na.rm = TRUE
  )
  message(sprintf(
    "Maximum social-jet-lag QC difference: %.1f minutes (reported only; no exclusions).",
    maximum_difference
  ))
  invisible(maximum_difference)
}

# This smaller table keeps only the identifiers and variables needed for
# N1-N5, P1 and D1. Diagnoses, medication, substance use and other descriptive
# fields are not added as unregistered nodes or predictors. Removing unused columns does not select participants by
# symptoms or outcomes.
# No categories are collapsed, no values are winsorised, and no participant or
# node is removed after looking at outliers, correlations or centrality. Those
# restrictions are stated in Questions 6 and 8.
assemble_prepared_data <- function(raw, scores) {
  data.frame(
    id = as.character(raw$id),
    visit = raw$visit,
    age = raw$age,
    isi_1 = raw$isi_1,
    isi_2 = raw$isi_2,
    isi_3 = raw$isi_3,
    isi_7 = raw$isi_7,
    ddnsi_total = scores$ddnsi_total,
    socialjetlag_h = scores$socialjetlag_h,
    cesd_no_sleep = scores$cesd_no_sleep,
    gad_total = scores$gad_total,
    inq_perceivedburden = scores$inq_perceivedburden,
    inq_thwartedbelong = scores$inq_thwartedbelong,
    upps_negativeurgency = scores$upps_negativeurgency,
    si_3mo_maxseverity = scores$si_3mo_maxseverity,
    psqi_c1 = scores$psqi_c1,
    isi_total = scores$isi_total,
    si_3mo_yn = scores$si_3mo_yn,
    cssrs_3mo_5 = raw$cssrs_3mo_5,
    cssrs_3mo_6 = raw$cssrs_3mo_6,
    cssrs_3mo_7 = raw$cssrs_3mo_7,
    cssrs_3mo_8 = raw$cssrs_3mo_8,
    check.names = FALSE
  )
}

# Questions 2 and 4 define the primary population from unrounded decimal age
# as 18 <= age < 26. This keeps 885 participants, including the 14 recorded
# between ages 25.05 and 25.96; an age <= 25 rule would wrongly remove them.
define_analysis_cohorts <- function(prepared) {
  primary_18_25 <- prepared[
    !is.na(prepared$age) & prepared$age >= 18 & prepared$age < 26,
    ,
    drop = FALSE
  ]

  # The age-range sensitivity changes only the upper-age restriction. It keeps
  # all 945 participants with observed age, including the 60 aged 26 or older,
  # but still excludes the 26 with missing age. Age is not added as a node or
  # covariate because that would change the network estimand (Question 2,
  # PDF pp. 5-6).
  n1_names <- c(
    "id", "visit", "age", "isi_1", "isi_2", "isi_3", "isi_7",
    "ddnsi_total", "socialjetlag_h", "cesd_no_sleep", "gad_total",
    "inq_perceivedburden", "inq_thwartedbelong", "upps_negativeurgency",
    "si_3mo_maxseverity"
  )
  age_sensitivity <- prepared[!is.na(prepared$age), n1_names, drop = FALSE]

  list(
    primary_18_25 = primary_18_25,
    age_sensitivity = age_sensitivity
  )
}

# These checks reproduce the cohort sizes and D1 event counts printed in
# Questions 1 and 2 before any model is fitted. A mismatch pauses the analysis
# instead of allowing a different data release to pass unnoticed.
validate_analysis_cohorts <- function(cohorts) {
  primary_18_25 <- cohorts$primary_18_25
  age_sensitivity <- cohorts$age_sensitivity
  stopifnot(
    nrow(primary_18_25) == 885L,
    nrow(age_sensitivity) == 945L,
    sum(primary_18_25$cssrs_3mo_5) == 13L,
    sum(primary_18_25$cssrs_3mo_6) == 20L,
    sum(primary_18_25$cssrs_3mo_7) == 93L,
    sum(primary_18_25$cssrs_3mo_8) == 16L
  )
  invisible(TRUE)
}

# This workbook contains the two analysis-ready datasets: Primary_18_25 is the
# main cohort and Age_Sensitivity is the all-observed-age N1 cohort.
write_prepared_datasets <- function(cohorts, output_file) {
  writexl::write_xlsx(
    list(
      Primary_18_25 = cohorts$primary_18_25,
      Age_Sensitivity = cohorts$age_sensitivity
    ),
    path = output_file,
    format_headers = TRUE
  )
  message("Saved: ", output_file)
  invisible(output_file)
}

prepare_answers_data <- function(input_file, output_file = NULL, write_output = FALSE) {
  raw <- read_and_check_raw_data(input_file)
  registered_scores <- score_registered_measures(raw)
  validate_registered_scores(raw, registered_scores)
  social_jetlag_audit <- audit_social_jetlag(raw)
  prepared <- assemble_prepared_data(raw, registered_scores)
  cohorts <- define_analysis_cohorts(prepared)
  validate_analysis_cohorts(cohorts)
  written_path <- NA_character_
  if (isTRUE(write_output)) {
    if (is.null(output_file) || length(output_file) != 1L || !nzchar(output_file)) stop("A valid output_file is required when write_output is TRUE.", call. = FALSE)
    written_path <- write_prepared_datasets(cohorts, output_file)
  }
  list(raw = raw, scores = registered_scores, social_jetlag_audit = social_jetlag_audit, prepared = prepared, cohorts = cohorts, output_file = written_path)
}

# 03 - Registered network, classifier, and reporting definitions -------------

# Question 1 describes a cross-sectional, undirected network of conditional
# associations. Its edges cannot establish temporal order, mediation,
# causation or treatment effects, a limitation repeated in Question 13. The
# separate functions below follow that scientific hierarchy. The active run at
# the end calls them in the same order you would explain the study.

record_analysis_setup <- function(cohorts) {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(ANALYSIS_SEED)
  dirs <- list(root = analysis_dir, objects = file.path(analysis_dir, "objects"))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  # These choices are reported because the preregistration does not fix them numerically.
  implementation_choices <- data.frame(
    item = c(
      "analysis_seed", "edge_interval_level", "ordinary_case_drop_grid",
      "resampled_rank_ties", "p1_exact_tuning_ties", "p1_pr_auc_definition",
      "percentile_quantile_algorithm", "custom_case_drop_rounding",
      "p1_secondary_contrast", "p1_inner_fold_pairing",
      "bootnet_retry_policy", "estimation_warning_policy", "parallel_core_count"
    ),
    value = c(
      ANALYSIS_SEED, EDGE_CI_LEVEL, "bootnet defaults: 0.05-0.75 at 10 levels",
      paste("minimum rank for tables/first-place proportions;", "base-R average ranks for Spearman stability correlations"),
      "first in alpha order and descending lambda; tie count reported",
      "PRROC integral PR-AUC",
      "R quantile type 7; edge plot is drawn from the exported table",
      "drop round(n * fraction) participants using R round-to-even",
      "extended minus base", "same inner splits for base and extended within each outer split",
      paste("bootnet 1.8 internal retry limit: up to 10 failed candidate attempts",
            "per requested draw; failures are counted and redrawn; returned draws",
            "must equal the requested total"),
      paste("failed/nonconverged/non-PD/singular/nonfinite warnings invalidate fits;",
            "informational warnings retained where exposed; bootnet may suppress",
            "informational custom-estimator warnings inside resamples"),
      N_CORES
    ),
    pdf_status = c(
      "numeric value not printed", "level not printed", "grid not printed",
      "rule not printed", "rule not printed", "estimator not printed",
      "algorithm not printed", "integer rounding not printed", "sign not printed",
      "not required", "retry/redraw rule not printed", "warning rule not printed", "not prescribed"
    ),
    stringsAsFactors = FALSE
  )
  package_versions <- data.frame(
    package = registered_packages,
    version = vapply(registered_packages, function(x) as.character(utils::packageVersion(x)), character(1)),
    stringsAsFactors = FALSE
  )
  analysis_record <- data.frame(
    item = c("run_started_at", "run_mode", "input_file", "analysis_seed",
             "n_primary", "n_age_sensitivity", "edge_bootstraps",
             "case_bootstraps", "custom_bootstraps_per_fraction", "p1_bootstraps", "cores"),
    value = c(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "full fit from raw CSV", input_file,
              ANALYSIS_SEED, nrow(cohorts$primary_18_25), nrow(cohorts$age_sensitivity),
              N_EDGE_BOOT, N_CASE_BOOT, N_CUSTOM_BOOT, N_P1_BOOT, N_CORES),
    stringsAsFactors = FALSE
  )
  deviation_log <- data.frame(
    date = "", reason = "No deviations recorded", affected_analysis = "",
    relevant_results_already_inspected = NA, classification = "",
    stringsAsFactors = FALSE
  )
  list(
    dirs = dirs,
    private_result_path = file.path(dirs$objects, "ANSWERS_registered_results_private.rds"),
    registered_packages = registered_packages,
    analysis_record = analysis_record,
    package_versions = package_versions,
    implementation_choices = implementation_choices,
    deviation_log = deviation_log
  )
}

# Define the exact node set and scale treatment for every network.
# Questions 1 and 7 fix these choices before analysis.
# N1 alone decides H1a-H1e and H2; all other networks are sensitivity checks.
define_registered_networks <- function() {

  vulnerability_nodes <- c(
    "cesd_no_sleep", "gad_total", "inq_perceivedburden",
    "inq_thwartedbelong", "upps_negativeurgency"
  )
  n1_sleep <- c(
    "isi_1", "isi_2", "isi_3", "isi_7", "ddnsi_total",
    "socialjetlag_h"
  )

  # These five communities come from theory, not from a clustering algorithm.
  # Question 12 uses them to distinguish within-domain from between-domain
  # connections and Question 21 uses them for bridge expected influence.
  make_domains <- function(sleep, ideation) {
    c(
      stats::setNames(rep("sleep_nightmare_circadian", length(sleep)), sleep),
      cesd_no_sleep = "emotional_distress",
      gad_total = "emotional_distress",
      inq_perceivedburden = "interpersonal_needs",
      inq_thwartedbelong = "interpersonal_needs",
      upps_negativeurgency = "impulsivity",
      stats::setNames("suicidality", ideation)
    )
  }

  # Each specification stores the cohort, nodes, scale treatment and domains
  # together. Passing this object through every later function helps prevent a
  # sensitivity network from accidentally inheriting the wrong settings.
  make_spec <- function(
    name, sleep, ideation, ordered_levels, correlation, cohort,
    hypothesis_bearing = FALSE
  ) {
    nodes <- c(sleep, vulnerability_nodes, ideation)
    domains <- make_domains(sleep, ideation)[nodes]
    list(
      name = name,
      sleep = sleep,
      vulnerability = vulnerability_nodes,
      ideation = ideation,
      nodes = nodes,
      ordered_levels = ordered_levels,
      correlation = correlation,
      cohort = cohort,
      domains = domains,
      hypothesis_bearing = hypothesis_bearing
    )
  }

  # Question 10 treats the four ISI responses and the highest
  # endorsed ideation category as ordered categories. N2 additionally treats
  # psqi_c1 as ordered, and N4 orders No below Yes. The remaining totals and
  # social jet lag are treated as approximately continuous.
  isi_ordered <- list(
    isi_1 = 0:4, isi_2 = 0:4, isi_3 = 0:4, isi_7 = 0:4,
    si_3mo_maxseverity = 0:5
  )

  specs <- list(
    # N1 is the primary 12-node mixed-correlation network described in
    # Question 1. Its six sleep-domain nodes, five vulnerability nodes and
    # ordered ideation endpoint are the only network inputs for H1 and H2.
    N1 = make_spec(
      "N1", n1_sleep, "si_3mo_maxseverity", isi_ordered,
      "mixed", "primary", TRUE
    ),
    # N2 replaces DDNSI nightmare burden with PSQI subjective sleep quality.
    # Because the construct and recall period both change, it is not a direct
    # DDNSI-versus-PSQI importance comparison (Question 1, PDF p. 3).
    N2 = make_spec(
      "N2",
      c("isi_1", "isi_2", "isi_3", "isi_7", "psqi_c1", "socialjetlag_h"),
      "si_3mo_maxseverity",
      c(isi_ordered, list(psqi_c1 = 0:3)),
      "mixed", "primary"
    ),
    # N3 replaces four selected ISI items with the full ISI total, which also
    # contains items 4-6. The total is never placed beside its constituent
    # items, avoiding deterministic overlap. The resulting nine-node network
    # examines sensitivity to representing insomnia by four terms rather than
    # one (Questions 1, 6 and 12)
    N3 = make_spec(
      "N3", c("isi_total", "ddnsi_total", "socialjetlag_h"),
      "si_3mo_maxseverity", list(si_3mo_maxseverity = 0:5),
      "mixed", "primary"
    ),
    # N4 changes only the ideation endpoint: ordered maximum severity becomes
    # the ordered binary any-ideation indicator, No < Yes.
    N4 = make_spec(
      "N4", n1_sleep, "si_3mo_yn",
      list(
        isi_1 = 0:4, isi_2 = 0:4, isi_3 = 0:4, isi_7 = 0:4,
        si_3mo_yn = 0:1
      ),
      "mixed", "primary"
    ),
    # N5 keeps the complete N1 node set and changes only the correlation matrix
    # to Spearman. It runs routinely—even if N1 is estimable—and remains a
    # sensitivity analysis rather than a replacement for N1 (Questions 1, 11
    # and 29)
    N5 = make_spec(
      "N5", n1_sleep, "si_3mo_maxseverity", list(),
      "spearman", "primary"
    ),
    # This additional N1 rerun uses all 945 records with observed age while
    # retaining the same nodes, scoring, mixed correlation and EBICglasso
    # settings described in Question 2.
    N1_age_sensitivity = make_spec(
      "N1_age_sensitivity", n1_sleep, "si_3mo_maxseverity", isi_ordered,
      "mixed", "age_sensitivity"
    )
  )

  stopifnot(
    length(specs$N1$nodes) == 12L,
    length(specs$N2$nodes) == 12L,
    length(specs$N3$nodes) == 9L,
    length(specs$N4$nodes) == 12L,
    length(specs$N5$nodes) == 12L
  )
  specs
}


# The next small functions each handle one part of the network analysis:
# missingness, scale coding, estimation, summaries or resampling.

# Check missingness, coding and possible redundancy ----------------------
# Questions 8 and 9 report complete required data and therefore specify no
# imputation. If a frozen file differs, complete-case estimation is allowed
# only when every required variable has at most 5% missingness and the total
# loss of eligible cases is also at most 5%.

network_missingness <- function(data, nodes, analysis_name) {
  missing_prop <- colMeans(is.na(data[nodes]))
  complete <- stats::complete.cases(data[nodes])
  case_loss <- 1 - mean(complete)
  table <- data.frame(
    analysis = analysis_name,
    variable = nodes,
    missing_n = colSums(is.na(data[nodes])),
    missing_proportion = as.numeric(missing_prop),
    stringsAsFactors = FALSE
  )
  if (max(missing_prop) > 0.05 || case_loss > 0.05) {
    stop(
      analysis_name,
      " exceeds the registered 5% missingness/case-loss limit. Pause and log a revision."
    )
  }
  list(table = table, complete = complete, case_loss = case_loss)
}

# Only variables named as ordinal in Question 10 are turned into ordered
# factors. An empty factor level may disappear in a resample, but populated
# categories are never combined. A node with fewer than two observed values
# cannot contribute a correlation and makes that resample non-estimable.
format_network_data <- function(data, spec) {
  out <- data[spec$nodes]
  if (length(spec$ordered_levels) > 0L) {
    for (name in names(spec$ordered_levels)) {
      allowed <- spec$ordered_levels[[name]]
      observed <- sort(unique(stats::na.omit(out[[name]])))
      if (!all(observed %in% allowed)) {
        stop(spec$name, ": invalid ordered values in ", name)
      }
      out[[name]] <- ordered(out[[name]], levels = allowed)
    }
  }
  out <- out[stats::complete.cases(out), , drop = FALSE]
  # Dropping an unused label is not category collapsing: every category that
  # is actually present keeps its registered value and order.
  out[] <- lapply(out, function(x) if (is.factor(x)) droplevels(x) else x)
  if (any(vapply(out, function(x) length(unique(x)) < 2L, logical(1)))) {
    stop(spec$name, ": at least one resampled node has insufficient variation")
  }
  out
}

# N5 needs numeric scores for Spearman ranks. Converting an ordered factor
# through its displayed value preserves 0-4 and 0-5 rather than replacing
# them with R's internal 1-based factor codes.
numeric_network_data <- function(data) {
  as.data.frame(lapply(data, function(x) {
    if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  }), check.names = FALSE)
}

# Question 12 asks for category frequencies, variance and possible
# redundancy to be examined if a matrix is not positive definite. These
# tables make that audit available before fitting; they are not used for
# data-driven node deletion. The recorded minimum and maximum values also
# make the range and coding-direction checks from Question 8 inspectable.
node_audit <- function(data, spec) {
  d <- data[spec$nodes]
  node_rows <- lapply(names(d), function(name) {
    x <- d[[name]]
    data.frame(
      network = spec$name,
      node = name,
      registered_ordered = name %in% names(spec$ordered_levels),
      missing_n = sum(is.na(x)),
      unique_n = length(unique(stats::na.omit(x))),
      minimum = if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
      maximum = if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
      variance = stats::var(as.numeric(x), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  frequency_rows <- lapply(names(d), function(name) {
    tab <- table(d[[name]], useNA = "ifany")
    data.frame(
      network = spec$name,
      node = name,
      value = names(tab),
      n = as.integer(tab),
      stringsAsFactors = FALSE
    )
  })
  pairs <- utils::combn(names(d), 2, simplify = FALSE)
  redundancy_rows <- lapply(pairs, function(pair) {
    x <- as.numeric(d[[pair[1]]])
    y <- as.numeric(d[[pair[2]]])
    ok <- stats::complete.cases(x, y)
    rho <- if (sum(ok) > 2L && stats::sd(x[ok]) > 0 && stats::sd(y[ok]) > 0) {
      suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
    } else {
      NA_real_
    }
    data.frame(
      network = spec$name,
      node1 = pair[1],
      node2 = pair[2],
      exactly_identical = sum(ok) > 0L && all(x[ok] == y[ok]),
      spearman = rho,
      stringsAsFactors = FALSE
    )
  })
  list(
    nodes = do.call(rbind, node_rows),
    frequencies = do.call(rbind, frequency_rows),
    redundancy = do.call(rbind, redundancy_rows)
  )
}


# Estimate the registered correlation and EBICglasso network --------------
# Questions 10-12 prescribe two correlation approaches. N1-N4 use cor_auto
# so explicitly ordered measures receive ordinal correlations and the scale
# totals and social jet lag contribute continuous correlations. N5 uses
# Spearman ranks on the same N1 scores as a robustness check for skewed and
# floor-concentrated variables.
#
# EBICglasso then estimates a sparse partial-correlation network. An edge is
# a signed conditional association after accounting for the other included
# nodes; it is not a temporal pathway, mediator or causal effect.

warning_requires_failure <- function(
  warnings, matrix_warnings_are_diagnostic = FALSE
) {
  # Informational warnings are retained for review. Warnings that say an
  # estimate failed, did not converge or became non-finite invalidate the fit.
  # A non-positive-definite or singular correlation matrix is handled slightly
  # differently: Question 12 asks us to calculate its diagnostics and nearPD
  # distance before reporting it as non-estimable.
  failed_estimation <- grepl(
    paste(
      "failed|non[- ]?conver|not finite|nan|infinite|could not|",
      "estimation error",
      sep = ""
    ),
    warnings,
    ignore.case = TRUE,
    perl = TRUE
  )
  matrix_problem <- grepl(
    "non[- ]?positive|not positive definite|singular",
    warnings,
    ignore.case = TRUE,
    perl = TRUE
  )
  if (matrix_warnings_are_diagnostic) {
    failed_estimation & !matrix_problem
  } else {
    failed_estimation | matrix_problem
  }
}

estimate_prespecified_network <- function(data, spec) {
  # Formatting can fail if complete cases disappear or a node has fewer than
  # two populated values in a resample. That draw is recorded as
  # non-estimable instead of terminating the remaining analyses.
  formatted <- tryCatch(
    format_network_data(data, spec),
    error = function(e) e
  )
  if (inherits(formatted, "error")) {
    return(list(
      estimable = FALSE,
      reason = paste(
        "network data formatting failed:",
        conditionMessage(formatted)
      ),
      n_cc = NA_integer_,
      warnings = character()
    ))
  }
  d <- formatted
  n_cc <- nrow(d)
  if (n_cc <= ncol(d)) {
    return(list(estimable = FALSE, reason = "n_cc <= number of nodes", n_cc = n_cc))
  }

  prefix_warnings <- function(prefix, warnings) {
    if (length(warnings) == 0L) character() else paste0(prefix, warnings)
  }
  correlation_warnings <- character()
  R <- tryCatch(
    withCallingHandlers(
      {
        if (identical(spec$correlation, "mixed")) {
          # This is the exact N1-N4 call printed in Question 12. Automatic
          # detection is off so only the explicitly declared factors are
          # treated as ordinal.
          qgraph::cor_auto(
            d,
            detectOrdinal = FALSE,
            forcePD = FALSE,
            missing = "listwise"
          )
        } else {
          # N5 retains the N1 scores and changes only the correlation method,
          # exactly matching the same-node sensitivity analysis in the form.
          stats::cor(
            numeric_network_data(d),
            method = "spearman",
            use = "complete.obs"
          )
        }
      },
      warning = function(w) {
        correlation_warnings <<- c(correlation_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(R, "error")) {
    return(list(
      estimable = FALSE,
      reason = paste("correlation estimation failed:", conditionMessage(R)),
      n_cc = n_cc,
      warnings = prefix_warnings("correlation: ", correlation_warnings)
    ))
  }
  correlation_warning_failed <- warning_requires_failure(
    correlation_warnings,
    matrix_warnings_are_diagnostic = TRUE
  )
  if (any(correlation_warning_failed)) {
    return(list(
      estimable = FALSE,
      reason = paste(
        "correlation estimation issued a material warning:",
        paste(
          unique(correlation_warnings[correlation_warning_failed]),
          collapse = " | "
        )
      ),
      n_cc = n_cc,
        warnings = prefix_warnings("correlation: ", correlation_warnings)
    ))
  }

  # EBICglasso needs a finite, full-rank, positive-definite correlation
  # matrix. These checks record whether that condition is met. Failure leads
  # to a transparent non-estimable result rather than an unregistered repair.
  R <- as.matrix(R)
  dimnames(R) <- list(spec$nodes, spec$nodes)
  finite <- all(is.finite(R))
  symmetric <- finite && isTRUE(all.equal(
    R, t(R), tolerance = 1e-10, check.attributes = FALSE
  ))
  eigenvalues <- if (symmetric) {
    eigen(R, symmetric = TRUE, only.values = TRUE)$values
  } else {
    rep(NA_real_, ncol(R))
  }
  matrix_rank <- if (finite) qr(R)$rank else NA_integer_
  min_eigen <- if (symmetric) min(eigenvalues) else NA_real_
  positive_definite <-
    finite && symmetric && matrix_rank == ncol(R) && min_eigen > 0

  diagnostic <- data.frame(
    network = spec$name,
    n_cc = n_cc,
    finite = finite,
    symmetric = symmetric,
    matrix_rank = matrix_rank,
    dimension = ncol(R),
    minimum_eigenvalue = min_eigen,
    positive_definite = positive_definite,
    stringsAsFactors = FALSE
  )

  near_pd <- NULL
  if (!positive_definite) {
    # Question 12 allows nearPD only to show how much mathematical adjustment
    # would be needed. Its adjusted matrix is never passed to EBICglasso, so
    # a failed registered matrix remains reported as non-estimable.
    near_pd <- tryCatch({
      npd <- Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)
      adjusted <- as.matrix(npd$mat)
      data.frame(
        network = spec$name,
        frobenius_adjustment = npd$normF,
        maximum_absolute_adjustment = max(abs(adjusted - R)),
        adjusted_minimum_eigenvalue = min(eigen(
          adjusted, symmetric = TRUE, only.values = TRUE
        )$values),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        network = spec$name,
        frobenius_adjustment = NA_real_,
        maximum_absolute_adjustment = NA_real_,
        adjusted_minimum_eigenvalue = NA_real_,
        stringsAsFactors = FALSE
      )
    })
    return(list(
      estimable = FALSE,
      reason = "correlation matrix is non-positive-definite or singular",
      n_cc = n_cc,
      data = d,
      correlation = R,
      diagnostic = diagnostic,
      near_pd = near_pd,
      warnings = prefix_warnings("correlation: ", correlation_warnings)
    ))
  }

  # The arguments below reproduce the exact Question 12 call: gamma 0.50,
  # no diagonal penalty, 100 lambdas, lambda ratio 0.01, a positive-
  # definiteness check, no refit and no extra edge thresholding.
  ebic_warnings <- character()
  graph <- tryCatch(
    withCallingHandlers(
      qgraph::EBICglasso(
        S = R,
        n = n_cc,
        gamma = 0.50,
        penalize.diagonal = FALSE,
        nlambda = 100,
        lambda.min.ratio = 0.01,
        checkPD = TRUE,
        refit = FALSE,
        threshold = FALSE
      ),
      warning = function(w) {
        ebic_warnings <<- c(ebic_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  all_estimation_warnings <- c(
    prefix_warnings("correlation: ", correlation_warnings),
    prefix_warnings("EBICglasso: ", ebic_warnings)
  )
  if (any(warning_requires_failure(ebic_warnings))) {
    return(list(
      estimable = FALSE,
      reason = paste(
        "EBICglasso issued a material warning:",
        paste(unique(ebic_warnings[warning_requires_failure(
          ebic_warnings
        )]), collapse = " | ")
      ),
      n_cc = n_cc,
      data = d,
      correlation = R,
      diagnostic = diagnostic,
      warnings = all_estimation_warnings
    ))
  }
  if (inherits(graph, "error")) {
    return(list(
      estimable = FALSE,
      reason = paste("EBICglasso failed:", conditionMessage(graph)),
      n_cc = n_cc,
      data = d,
      correlation = R,
      diagnostic = diagnostic,
      warnings = all_estimation_warnings
    ))
  }

  graph <- as.matrix(graph)
  dimnames(graph) <- list(spec$nodes, spec$nodes)
  if (!all(is.finite(graph)) || !isTRUE(all.equal(graph, t(graph), tolerance = 1e-10))) {
    return(list(
      estimable = FALSE,
      reason = "EBICglasso returned a non-finite or asymmetric matrix",
      n_cc = n_cc,
      data = d,
      correlation = R,
      diagnostic = diagnostic,
      warnings = all_estimation_warnings
    ))
  }

  list(
    estimable = TRUE,
    reason = "estimable",
    n_cc = n_cc,
    data = d,
    correlation = R,
    graph = graph,
    diagnostic = diagnostic,
    near_pd = near_pd,
    warnings = all_estimation_warnings
  )
}


# Calculate the full-sample network summaries ------------------------------
# These quantities describe cross-sectional conditional associations.
# In line with Questions 1, 13 and 21-22, none is interpreted as a causal
# pathway, mediator, treatment target or individual-risk prediction. A zero
# EBICglasso edge is also not proof that no association exists, because
# regularisation can shrink weak estimates to zero.

# An undirected network matrix stores the same edge above and below its
# diagonal. Keeping only the upper triangle reports each unique edge once.
edge_table <- function(graph) {
  index <- which(upper.tri(graph), arr.ind = TRUE)
  data.frame(
    node1 = rownames(graph)[index[, 1]],
    node2 = colnames(graph)[index[, 2]],
    weight = graph[index],
    stringsAsFactors = FALSE
  )
}

# Questions 12, 21 and 22 define the custom H2 summary. For each of the five
# vulnerability nodes, A is the mean signed edge to that network's sleep
# nodes, B is the signed edge to its ideation node and T = (A + B) / 2.
# Averaging the sleep edges before combining the components prevents the
# larger sleep domain from receiving extra weight simply because it contains
# more nodes. A and B keep their signs, so they can offset each other; a
# negative T means net negative connectivity, not a protective effect.
targeted_connectivity <- function(graph, spec) {
  rows <- lapply(spec$vulnerability, function(v) {
    A <- mean(graph[v, spec$sleep])
    B <- graph[v, spec$ideation]
    data.frame(
      node = v,
      A_sleep_mean_signed = A,
      B_ideation_signed = B,
      T_targeted = (A + B) / 2,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  # Rank 1 is the largest T. Using the minimum rank for a tie makes every
  # tied leader visibly rank first, while unique_first separately records
  # whether only one candidate occupies that position.
  out$rank <- rank(-out$T_targeted, ties.method = "min")
  out$shared_first <- out$rank == 1L
  out$unique_first <- out$shared_first & sum(out$shared_first) == 1L
  out
}

# Strength adds the absolute size of all edges attached to a node. Expected
# influence keeps the signs, so positive and negative connections can cancel.
# Question 21 labels both as secondary summaries rather than H2 rules.
ordinary_centrality <- function(graph) {
  data.frame(
    node = rownames(graph),
    strength = rowSums(abs(graph)),
    expected_influence = rowSums(graph),
    stringsAsFactors = FALSE
  )
}

# One-step bridge expected influence summarizes a node's signed connections
# to nodes in the other theory-defined domains. It is normalized because the
# five domains have unequal sizes. This secondary measure is not mediation
# and does not show what would happen after an intervention (Question 21)
bridge_expected_influence <- function(graph, spec) {
  bridge <- networktools::bridge(
    graph,
    communities = unname(spec$domains[colnames(graph)]),
    nodes = colnames(graph),
    normalize = TRUE
  )
  values <- bridge[["Bridge Expected Influence (1-step)"]]
  data.frame(
    node = names(values),
    normalized_bridge_expected_influence_1step = as.numeric(values),
    community = unname(spec$domains[names(values)]),
    stringsAsFactors = FALSE
  )
}

# Questions 17-20 and 24-28 are marked N/A in the form. Accordingly, the
# script does not add group comparisons, formal tests between individual
# edge weights, data-driven clustering or cluster-stability analyses.

plot_network <- function(graph, spec) {
  # The graph is a visual summary only. The spring layout and colours help a
  # reader see nodes and domains, but their positions have no decision-rule
  # meaning; Question 14 asks conclusions to rely on numerical estimates and
  # resampling uncertainty rather than appearance alone.
  groups <- split(seq_along(spec$nodes), unname(spec$domains[spec$nodes]))
  qgraph::qgraph(
    graph,
    layout = "spring",
    groups = groups,
    labels = spec$nodes,
    theme = "colorblind",
    cut = 0,
    minimum = 0,
    legend = TRUE,
    legend.cex = 0.45,
    label.cex = 0.55,
    title = spec$name,
    details = TRUE
  )
}


# Assess edge accuracy and ordinary-centrality stability -------------------
# Question 14 uses 1,000 nonparametric participant bootstraps to show how
# much each edge estimate changes across resampled datasets. Question 23
# uses a different, case-dropping bootstrap for strength and expected
# influence: this asks whether reduced-sample centrality values remain
# strongly correlated with their full-sample values when part of the sample
# is removed. Edge uncertainty and case-drop stability therefore answer
# related but different questions.

bootnet_estimator <- function(data, spec) {
  # bootnet calls this function on the participant rows selected in each
  # resample. The ordered factors, correlation matrix and EBICglasso network
  # are rebuilt every time, as required by Question 12. Bootstrapping one
  # fixed correlation matrix would not repeat the registered analysis.
  tryCatch({
    fit <- estimate_prespecified_network(data, spec)
    if (!isTRUE(fit$estimable)) {
      stop(fit$reason, call. = FALSE)
    }
    # Material warnings have already made the fit non-estimable above.
    # bootnet may suppress informational estimator warnings, so the script
    # does not claim that its resample-warning table is exhaustive.
    list(graph = fit$graph)
  }, error = function(e) {
    bootnet_retry_audit$failed_attempts <-
      bootnet_retry_audit$failed_attempts + 1L
    bootnet_retry_audit$failure_reasons <- c(
      bootnet_retry_audit$failure_reasons,
      conditionMessage(e)
    )
    stop(e)
  })
}

# bootnet 1.8 does not emit a warning containing its retry count. Because this
# script deliberately uses one core, this small in-process audit environment
# records every rejected candidate draw at the estimator boundary instead.
bootnet_retry_audit <- new.env(parent = emptyenv())
reset_bootnet_retry_audit <- function() {
  bootnet_retry_audit$failed_attempts <- 0L
  bootnet_retry_audit$failure_reasons <- character()
  invisible(NULL)
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

format_bootnet_retry_reasons <- function(reasons) {
  if (length(reasons) == 0L) return(character())
  counts <- sort(table(reasons), decreasing = TRUE)
  sprintf(
    "rejected candidate draw retried (%d): %s",
    as.integer(counts),
    names(counts)
  )
}

canonical_edge_key <- function(node1, node2) {
  paste(pmin(node1, node2), pmax(node1, node2), sep = "\r")
}

edge_bootstrap_summary <- function(boot, graph, network_name) {
  # Every unique edge receives its full-sample estimate, percentile interval
  # and proportion of successful bootstrap estimates above zero. These are
  # uncertainty descriptions, not ordinary significance-test intervals.
  bt <- as.data.frame(boot$bootTable)
  bt <- bt[bt$type == "edge", , drop = FALSE]
  bt$key <- canonical_edge_key(bt$node1, bt$node2)

  original <- edge_table(graph)
  original$key <- canonical_edge_key(original$node1, original$node2)
  alpha_tail <- (1 - EDGE_CI_LEVEL) / 2

  rows <- lapply(seq_len(nrow(original)), function(i) {
    values <- bt$value[bt$key == original$key[i]]
    if (length(values) == 0L) {
      stop(network_name, ": an edge is absent from bootnet$bootTable")
    }
    interval <- stats::quantile(
      values,
      probs = c(alpha_tail, 1 - alpha_tail),
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
    data.frame(
      network = network_name,
      node1 = original$node1[i],
      node2 = original$node2[i],
      estimate = original$weight[i],
      percentile_lower = interval[1],
      percentile_upper = interval[2],
      proportion_above_zero = mean(values > 0, na.rm = TRUE),
      successful_bootstraps = sum(is.finite(values)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

draw_edge_accuracy_plot <- function(intervals, spec) {
  # Plotting from the exported table keeps its visual and numeric intervals
  # identical. N1's five H1 edges are highlighted for the reader, but the
  # highlight and interval do not change the point-estimate-only H1 rule.
  h1_nodes <- c("ddnsi_total", "isi_1", "isi_2", "isi_3", "isi_7")
  plot_data <- intervals
  plot_data$edge <- paste(plot_data$node1, plot_data$node2, sep = " -- ")
  plot_data$is_H1_edge <- spec$name == "N1" & (
    (plot_data$node1 == "si_3mo_maxseverity" & plot_data$node2 %in% h1_nodes) |
      (plot_data$node2 == "si_3mo_maxseverity" & plot_data$node1 %in% h1_nodes)
  )
  plot_data$edge <- factor(
    plot_data$edge,
    levels = plot_data$edge[order(plot_data$estimate)]
  )
  print(
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = edge, y = estimate, colour = is_H1_edge)
    ) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
      ggplot2::geom_errorbar(ggplot2::aes(
        ymin = percentile_lower,
        ymax = percentile_upper
      ), width = 0.25) +
      ggplot2::geom_point(size = 1.4) +
      ggplot2::coord_flip() +
      ggplot2::scale_colour_manual(
        values = c(`FALSE` = "#4D4D4D", `TRUE` = "#D55E00"),
        labels = c(`FALSE` = "Other edge", `TRUE` = "Registered H1 edge")
      ) +
      ggplot2::labs(
        title = paste0(spec$name, ": bootstrap edge accuracy"),
        x = NULL,
        y = paste0(round(100 * EDGE_CI_LEVEL), "% percentile interval"),
        colour = NULL
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(legend.position = "bottom")
  )
}

draw_case_bootnet_plot <- function(boot) {
  print(plot(boot, statistics = c("strength", "expectedInfluence")))
}

run_bootnet_resamples <- function(data, spec, graph) {
  # Both bootstrap types repeat the same correlation and EBICglasso estimate
  # used in the full sample. bootnet retries an estimator internally when a
  # particular resample cannot be estimated. Its `maxErrors` argument is not
  # used for an already estimated bootnetResult in bootnet 1.8, so it is not
  # passed here. Failed attempts are counted at the estimator boundary and
  # reported in the final workbook.
  estimated <- bootnet::estimateNetwork(
    data,
    default = "none",
    fun = bootnet_estimator,
    spec = spec,
    labels = spec$nodes,
    weighted = TRUE,
    signed = TRUE,
    directed = FALSE,
    verbose = FALSE
  )

  # The nonparametric bootstrap samples participants with replacement and
  # provides the 1,000 edge-accuracy replicates required by Question 14.
  reset_bootnet_retry_audit()
  edge_run <- capture_warnings(bootnet::bootnet(
    estimated,
    nBoots = N_EDGE_BOOT,
    type = "nonparametric",
    statistics = "edge",
    nCores = N_CORES,
    memorysaver = TRUE,
    verbose = FALSE
  ))
  edge_failed_attempts <- bootnet_retry_audit$failed_attempts
  edge_failure_reasons <- bootnet_retry_audit$failure_reasons
  edge_boot <- edge_run$value
  edge_intervals <- edge_bootstrap_summary(edge_boot, graph, spec$name)

  # The ordinary case-drop analysis removes increasing proportions of cases
  # and recalculates strength and expected influence. The CS coefficient asks
  # how much data can be removed while retaining correlation 0.70 with the
  # full-sample values. Question 23 treats values above 0.25 as minimally
  # interpretable and above 0.50 as preferable.
  reset_bootnet_retry_audit()
  case_run <- capture_warnings(bootnet::bootnet(
    estimated,
    nBoots = N_CASE_BOOT,
    type = "case",
    statistics = c("strength", "expectedInfluence"),
    caseMin = 0.05,
    caseMax = 0.75,
    caseN = 10,
    nCores = N_CORES,
    memorysaver = TRUE,
    verbose = FALSE
  ))
  case_failed_attempts <- bootnet_retry_audit$failed_attempts
  case_failure_reasons <- bootnet_retry_audit$failure_reasons
  case_boot <- case_run$value
  cs <- bootnet::corStability(
    case_boot,
    cor = 0.70,
    statistics = c("strength", "expectedInfluence"),
    verbose = FALSE
  )
  cs_table <- data.frame(
    network = spec$name,
    statistic = names(cs),
    cs_cor_0_70 = as.numeric(cs),
    minimally_interpretable_above_0_25 = as.numeric(cs) > 0.25,
    preferable_above_0_50 = as.numeric(cs) > 0.50,
    stringsAsFactors = FALSE
  )

  retry_table <- data.frame(
    network = spec$name,
    bootstrap = c("edge_nonparametric", "ordinary_case_drop"),
    requested_successes = c(N_EDGE_BOOT, N_CASE_BOOT),
    returned_successes = c(
      length(unique(edge_boot$bootTable$name)),
      length(unique(case_boot$bootTable$name))
    ),
    failed_attempts_resampled = c(
      edge_failed_attempts,
      case_failed_attempts
    ),
    stringsAsFactors = FALSE
  )
  if (any(retry_table$returned_successes != retry_table$requested_successes)) {
    stop(
      spec$name,
      ": bootnet returned fewer successful draws than requested",
      call. = FALSE
    )
  }

  edge_warnings <- c(
    edge_run$warnings,
    format_bootnet_retry_reasons(edge_failure_reasons)
  )
  case_warnings <- c(
    case_run$warnings,
    format_bootnet_retry_reasons(case_failure_reasons)
  )
  all_warnings <- c(edge_warnings, case_warnings)
  warning_table <- if (length(all_warnings) > 0L) {
    data.frame(
      network = rep(spec$name, length(all_warnings)),
      bootstrap = c(
        rep("edge_nonparametric", length(edge_warnings)),
        rep("ordinary_case_drop", length(case_warnings))
      ),
      warning = all_warnings,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      network = character(), bootstrap = character(), warning = character(),
      stringsAsFactors = FALSE
    )
  }

  list(
    estimated = estimated,
    edge = edge_boot,
    edge_intervals = edge_intervals,
    case = case_boot,
    cs = cs_table,
    retries = retry_table,
    warnings = warning_table
  )
}


# Assess the custom index and bridge stability ------------------------------
# Standard bootnet summaries do not include A, B, T or normalized bridge
# expected influence. Questions 12 and 23 therefore specify a separate loop
# that removes 10%, 25%, 50% and 75% of participants and re-estimates the
# complete network 1,000 times at each level. The 75% condition is a stress
# test, not the sole basis for judging stability, and failed fits are kept in
# the audit rather than silently discarded.

quantile_row <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(c(mean = NA, sd = NA, p05 = NA, median = NA, p95 = NA))
  }
  qs <- stats::quantile(
    values,
    probs = c(0.05, 0.50, 0.95),
    names = FALSE,
    type = 7
  )
  c(
    mean = mean(values),
    sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
    p05 = qs[1],
    median = qs[2],
    p95 = qs[3]
  )
}

run_custom_case_drop <- function(data, spec, full_graph) {
  # The full-sample T values provide the reference ordering against which
  # each reduced-sample ordering is compared.
  full_target <- targeted_connectivity(full_graph, spec)
  full_T <- stats::setNames(full_target$T_targeted, full_target$node)

  component_rows <- list()
  bridge_rows <- list()
  rank_rows <- list()
  failure_rows <- list()
  warning_rows <- list()
  component_i <- bridge_i <- rank_i <- failure_i <- warning_i <- 0L

  for (drop_fraction in CUSTOM_DROP_FRACTIONS) {
    # The form specifies percentages but not how to turn a percentage into a
    # whole number of people. Rounding n*fraction is therefore disclosed as
    # an implementation choice, and both kept_n and dropped_n are reported.
    drop_n <- round(nrow(data) * drop_fraction)
    keep_n <- nrow(data) - drop_n

    for (replicate in seq_len(N_CUSTOM_BOOT)) {
      keep <- sample.int(nrow(data), size = keep_n, replace = FALSE)
      fit <- tryCatch(
        estimate_prespecified_network(data[keep, , drop = FALSE], spec),
        error = function(e) list(
          estimable = FALSE,
          reason = paste("unexpected error:", conditionMessage(e))
        )
      )

      if (!isTRUE(fit$estimable)) {
        failure_i <- failure_i + 1L
        failure_rows[[failure_i]] <- data.frame(
          network = spec$name,
          drop_fraction = drop_fraction,
          replicate = replicate,
          kept_n = keep_n,
          reason = fit$reason,
          stringsAsFactors = FALSE
        )
        next
      }

      if (length(fit$warnings) > 0L) {
        for (warning_text in unique(fit$warnings)) {
          warning_i <- warning_i + 1L
          warning_rows[[warning_i]] <- data.frame(
            network = spec$name,
            drop_fraction = drop_fraction,
            replicate = replicate,
            warning = warning_text,
            stringsAsFactors = FALSE
          )
        }
      }

      metrics <- tryCatch(
        {
          target_value <- targeted_connectivity(fit$graph, spec)
          bridge_value <- bridge_expected_influence(fit$graph, spec)
          if (
            !setequal(target_value$node, spec$vulnerability) ||
              any(!is.finite(as.matrix(target_value[c(
                "A_sleep_mean_signed", "B_ideation_signed", "T_targeted"
              )]))) ||
              !setequal(bridge_value$node, spec$nodes) ||
              any(!is.finite(
                bridge_value$normalized_bridge_expected_influence_1step
              ))
          ) {
            stop("targeted or bridge output is incomplete/non-finite")
          }
          list(target = target_value, bridge = bridge_value)
        },
        error = function(e) e
      )
      if (inherits(metrics, "error")) {
        failure_i <- failure_i + 1L
        failure_rows[[failure_i]] <- data.frame(
          network = spec$name,
          drop_fraction = drop_fraction,
          replicate = replicate,
          kept_n = keep_n,
          reason = paste("metric extraction failed:", conditionMessage(metrics)),
          stringsAsFactors = FALSE
        )
        next
      }

      target <- metrics$target
      target_T <- stats::setNames(target$T_targeted, target$node)
      # This Spearman correlation compares the five vulnerability-node ranks
      # in the reduced sample with their full-sample ordering. It may be
      # undefined when regularisation creates too many ties; the script
      # reports how many correlations are defined instead of replacing an
      # undefined value with zero (Question 23, PDF p. 21).
      rank_correlation <- suppressWarnings(stats::cor(
        full_T[names(target_T)],
        target_T,
        method = "spearman"
      ))

      target_long <- rbind(
        data.frame(node = target$node, component = "A_sleep_mean_signed",
                   value = target$A_sleep_mean_signed),
        data.frame(node = target$node, component = "B_ideation_signed",
                   value = target$B_ideation_signed),
        data.frame(node = target$node, component = "T_targeted",
                   value = target$T_targeted)
      )
      component_i <- component_i + 1L
      component_rows[[component_i]] <- transform(
        target_long,
        network = spec$name,
        drop_fraction = drop_fraction,
        replicate = replicate,
        kept_n = keep_n,
        rank = rep(target$rank, 3L),
        shared_first = rep(target$shared_first, 3L),
        unique_first = rep(target$unique_first, 3L)
      )

      bridge <- metrics$bridge
      bridge_i <- bridge_i + 1L
      bridge_rows[[bridge_i]] <- transform(
        bridge,
        network = spec$name,
        drop_fraction = drop_fraction,
        replicate = replicate,
        kept_n = keep_n
      )

      rank_i <- rank_i + 1L
      rank_rows[[rank_i]] <- data.frame(
        network = spec$name,
        drop_fraction = drop_fraction,
        replicate = replicate,
        kept_n = keep_n,
        spearman_T_vs_full = rank_correlation,
        stringsAsFactors = FALSE
      )
    }
  }

  components <- if (length(component_rows)) do.call(rbind, component_rows) else NULL
  bridges <- if (length(bridge_rows)) do.call(rbind, bridge_rows) else NULL
  ranks <- if (length(rank_rows)) do.call(rbind, rank_rows) else NULL
  failures <- if (length(failure_rows)) do.call(rbind, failure_rows) else {
    data.frame(
      network = character(), drop_fraction = numeric(), replicate = integer(),
      kept_n = integer(), reason = character(), stringsAsFactors = FALSE
    )
  }
  warnings <- if (length(warning_rows)) do.call(rbind, warning_rows) else {
    data.frame(
      network = character(), drop_fraction = numeric(), replicate = integer(),
      warning = character(), stringsAsFactors = FALSE
    )
  }

  if (is.null(components) || is.null(bridges) || is.null(ranks)) {
    stop(spec$name, ": no custom case-drop re-estimation succeeded")
  }

  # The detailed draws are condensed into readable distributions of A, B,
  # T and bridge expected influence for each removal fraction.
  component_groups <- split(
    components,
    interaction(
      components$drop_fraction, components$node, components$component,
      drop = TRUE, lex.order = TRUE
    )
  )
  component_summary <- do.call(rbind, lapply(component_groups, function(x) {
    data.frame(
      network = spec$name,
      drop_fraction = x$drop_fraction[1],
      node = x$node[1],
      component = x$component[1],
      successful_repetitions = length(unique(x$replicate)),
      as.list(quantile_row(x$value)),
      stringsAsFactors = FALSE
    )
  }))

  bridge_groups <- split(
    bridges,
    interaction(
      bridges$drop_fraction, bridges$node,
      drop = TRUE, lex.order = TRUE
    )
  )
  bridge_summary <- do.call(rbind, lapply(bridge_groups, function(x) {
    data.frame(
      network = spec$name,
      drop_fraction = x$drop_fraction[1],
      node = x$node[1],
      community = x$community[1],
      successful_repetitions = length(unique(x$replicate)),
      as.list(quantile_row(x$normalized_bridge_expected_influence_1step)),
      stringsAsFactors = FALSE
    )
  }))

  # Rank-first proportions and the later rank-distribution table make close
  # or tied leaders visible. This is important because H2 requires perceived
  # burdensomeness to be uniquely first, not merely close to first.
  target_only <- components[components$component == "T_targeted", , drop = FALSE]
  rank_groups <- split(
    target_only,
    interaction(target_only$drop_fraction, target_only$node, drop = TRUE)
  )
  rank_first_summary <- do.call(rbind, lapply(rank_groups, function(x) {
    data.frame(
      network = spec$name,
      drop_fraction = x$drop_fraction[1],
      node = x$node[1],
      successful_network_fits = nrow(x),
      mean_rank = mean(x$rank),
      shared_first_proportion = mean(x$shared_first),
      unique_first_proportion = mean(x$unique_first),
      tied_first_proportion = mean(x$shared_first & !x$unique_first),
      stringsAsFactors = FALSE
    )
  }))

  rank_correlation_summary <- do.call(rbind, lapply(
    split(ranks, ranks$drop_fraction),
    function(x) {
      values <- x$spearman_T_vs_full
      finite_values <- values[is.finite(values)]
      summaries <- if (length(finite_values) > 0L) {
        c(
          mean = mean(finite_values),
          median = stats::median(finite_values),
          p05 = stats::quantile(
            finite_values, 0.05, names = FALSE, type = 7
          )
        )
      } else {
        c(mean = NA_real_, median = NA_real_, p05 = NA_real_)
      }
      data.frame(
        network = spec$name,
        drop_fraction = x$drop_fraction[1],
        successful_network_fits = nrow(x),
        defined_rank_correlations = length(finite_values),
        mean_spearman = summaries["mean"],
        median_spearman = summaries["median"],
        fifth_percentile_spearman = summaries["p05"],
        stringsAsFactors = FALSE
      )
    }
  ))

  success_counts <- as.data.frame(table(
    factor(ranks$drop_fraction, levels = CUSTOM_DROP_FRACTIONS)
  ), stringsAsFactors = FALSE)
  names(success_counts) <- c("drop_fraction", "successful_repetitions")
  success_counts$drop_fraction <- as.numeric(as.character(success_counts$drop_fraction))
  success_counts$network <- spec$name
  success_counts$original_n <- nrow(data)
  success_counts$dropped_n <- round(
    success_counts$original_n * success_counts$drop_fraction
  )
  success_counts$kept_n <- success_counts$original_n - success_counts$dropped_n
  success_counts$requested_repetitions <- N_CUSTOM_BOOT
  success_counts$failed_repetitions <- N_CUSTOM_BOOT - success_counts$successful_repetitions
  success_counts <- success_counts[c(
    "network", "drop_fraction", "original_n", "dropped_n", "kept_n",
    "requested_repetitions",
    "successful_repetitions", "failed_repetitions"
  )]

  list(
    components = components,
    bridges = bridges,
    rank_correlations = ranks,
    failures = failures,
    warnings = warnings,
    component_summary = component_summary,
    bridge_summary = bridge_summary,
    rank_first_summary = rank_first_summary,
    rank_correlation_summary = rank_correlation_summary,
    success_counts = success_counts
  )
}


# Estimate and assess one registered network -------------------------------
# Each specification follows the same sequence: choose its registered
# cohort, check missingness and coding, estimate its correlation matrix and
# EBICglasso network, calculate full-sample summaries, and then run both
# forms of stability analysis. A non-estimable network is recorded honestly
# rather than being replaced by a different method.

run_network <- function(spec, source_data, spec_index) {
  # The audit is completed before estimation so a coding or missing-data
  # problem is not diagnosed after looking at the network result.
  missingness <- network_missingness(source_data, spec$nodes, spec$name)
  audit <- node_audit(source_data, spec)
  fit <- estimate_prespecified_network(source_data, spec)

  if (!isTRUE(fit$estimable)) {
    return(list(
      spec = spec,
      estimable = FALSE,
      reason = fit$reason,
      fit = fit,
      missingness = missingness,
      audit = audit
    ))
  }

  # These four summaries have distinct roles: edges include the H1 targets;
  # targeted connectivity supplies H2; ordinary centrality and bridge
  # expected influence are secondary descriptions.
  edges <- transform(edge_table(fit$graph), network = spec$name)
  targeted <- transform(targeted_connectivity(fit$graph, spec), network = spec$name)
  centrality <- transform(ordinary_centrality(fit$graph), network = spec$name)
  bridge <- transform(bridge_expected_influence(fit$graph, spec), network = spec$name)

  set.seed(stage_seed(1000 * spec_index))
  message(spec$name, ": running registered bootnet resamples")
  boot <- run_bootnet_resamples(fit$data, spec, fit$graph)
  set.seed(stage_seed(50000 + 1000 * spec_index))
  message(spec$name, ": running targeted/bridge case-drop resamples")
  custom <- run_custom_case_drop(fit$data, spec, fit$graph)

  list(
    spec = spec,
    estimable = TRUE,
    reason = "estimable",
    fit = fit,
    missingness = missingness,
    audit = audit,
    edges = edges,
    targeted = targeted,
    centrality = centrality,
    bridge = bridge,
    boot = boot,
    custom = custom
  )
}


# Estimate all registered networks -----------------------------------------
# N1 is estimated first because it alone answers H1 and H2. N2-N5 then test
# specific sensitivities, and the observed-age N1 rerun checks the age rule.
estimate_registered_networks <- function(specs, cohorts, analysis_setup) {
  dirs <- analysis_setup$dirs
  expected_networks <- c("N1", "N2", "N3", "N4", "N5", "N1_age_sensitivity")
  if (!identical(names(specs), expected_networks)) {
    stop("The six registered network specifications are not present in their required order.", call. = FALSE)
  }

  network_results <- vector("list", length(specs))
  names(network_results) <- names(specs)
  for (i in seq_along(specs)) {
    network_name <- names(specs)[i]
    source_data <- if (identical(specs[[i]]$cohort, "primary")) {
      cohorts$primary_18_25
    } else {
      cohorts$age_sensitivity
    }
    network_results[[i]] <- run_network(specs[[i]], source_data, i)
    if (!isTRUE(network_results[[i]]$estimable)) {
      stop(
        network_name, " is not estimable: ", network_results[[i]]$reason,
        call. = FALSE
      )
    }

    message(network_name, ": complete")
  }

  # Each of the six networks contributes its graph, edge-accuracy plot, and case-drop plot.
  figure_path <- file.path(dirs$root, "ANSWERS_registered_figures.pdf")
  grDevices::pdf(figure_path, width = 11, height = 8.5, useDingbats = FALSE)
  device_open <- TRUE
  on.exit(if (device_open) grDevices::dev.off(), add = TRUE)
  for (result in network_results) {
    plot_network(result$fit$graph, result$spec)
    draw_edge_accuracy_plot(result$boot$edge_intervals, result$spec)
    draw_case_bootnet_plot(result$boot$case)
  }
  grDevices::dev.off()
  device_open <- FALSE
  network_status <- do.call(rbind, lapply(network_results, function(x) {
    data.frame(
      network = x$spec$name,
      role = if (isTRUE(x$spec$hypothesis_bearing)) {
        "primary hypothesis-bearing"
      } else if (identical(x$spec$name, "N1_age_sensitivity")) {
        "age-range sensitivity"
      } else {
        "prespecified sensitivity"
      },
      cohort = x$spec$cohort,
      estimable = x$estimable,
      reason = x$reason,
      n_cc = if (!is.null(x$fit$n_cc)) x$fit$n_cc else NA_integer_,
      node_count = length(x$spec$nodes),
      possible_edges = choose(length(x$spec$nodes), 2),
      full_sample_warning_count = length(unique(x$fit$warnings)),
      stringsAsFactors = FALSE
    )
  }))

  list(
    network_results = network_results,
    network_status = network_status,
    figure_path = figure_path
  )
}


# Turn the N1 estimates into the prespecified H1 and H2 answers -------------
answer_network_questions <- function(network_results) {

  # Question 16 maps H1a-H1e to five separate N1 edges: ideation with DDNSI,
  # ISI 1, ISI 2, ISI 3 and ISI 7. A positive full-sample N1 estimate is the
  # sole directional-consistency rule for each component. The interval and
  # bootstrap proportion above zero describe uncertainty but are not extra
  # significance tests. No all-five, majority or other omnibus H1 result is
  # created, and no multiplicity test is added to these directional estimates.
  h1_other_nodes <- c(
    H1a = "ddnsi_total",
    H1b = "isi_1",
    H1c = "isi_2",
    H1d = "isi_3",
    H1e = "isi_7"
  )

  make_h1_table <- function(result, decision_network = FALSE) {
    if (!isTRUE(result$estimable)) {
      return(data.frame(
        network = result$spec$name,
        hypothesis = names(h1_other_nodes),
        sleep_endpoint = unname(h1_other_nodes),
        ideation_endpoint = "si_3mo_maxseverity",
        estimate = NA_real_,
        percentile_lower = NA_real_,
        percentile_upper = NA_real_,
        proportion_above_zero = NA_real_,
        directionally_consistent = NA,
        decision_role = if (decision_network) "decision" else "descriptive only",
        stringsAsFactors = FALSE
      ))
    }

    intervals <- result$boot$edge_intervals
    rows <- lapply(seq_along(h1_other_nodes), function(i) {
      node <- unname(h1_other_nodes[i])
      key <- canonical_edge_key(node, "si_3mo_maxseverity")
      hit <- intervals[
        canonical_edge_key(intervals$node1, intervals$node2) == key,
        , drop = FALSE
      ]
      if (nrow(hit) != 1L) stop(result$spec$name, ": H1 edge extraction failed")
      data.frame(
        network = result$spec$name,
        hypothesis = names(h1_other_nodes)[i],
        sleep_endpoint = node,
        ideation_endpoint = "si_3mo_maxseverity",
        estimate = hit$estimate,
        percentile_lower = hit$percentile_lower,
        percentile_upper = hit$percentile_upper,
        proportion_above_zero = hit$proportion_above_zero,
        directionally_consistent = if (decision_network) hit$estimate > 0 else NA,
        decision_role = if (decision_network) "decision" else "descriptive only",
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }

  h1_decisions <- make_h1_table(network_results$N1, decision_network = TRUE)
  h1_age_descriptive <- make_h1_table(
    network_results$N1_age_sensitivity,
    decision_network = FALSE
  )
  # Question 2 makes the all-observed-age rerun descriptive. It compares the
  # directions and magnitudes of the five H1 edges but cannot change the N1
  # decisions or establish equivalence between age groups.
  h1_age_comparison <- merge(
    h1_decisions,
    h1_age_descriptive,
    by = c("hypothesis", "sleep_endpoint", "ideation_endpoint"),
    suffixes = c("_primary_N1", "_age_sensitivity"),
    sort = FALSE
  )
  h1_age_comparison$direction_same <- with(
    h1_age_comparison,
    sign(estimate_primary_N1) == sign(estimate_age_sensitivity)
  )
  h1_age_comparison$magnitude_difference_age_minus_primary <- with(
    h1_age_comparison,
    estimate_age_sensitivity - estimate_primary_N1
  )

  # Questions 16 and 22 support H2 only when perceived burdensomeness has a
  # strictly larger full-sample N1 T value than depression, anxiety, thwarted
  # belongingness and negative urgency. A tied first place is non-support;
  # resampled ranks remain available to show whether the ordering is fragile.
  if (isTRUE(network_results$N1$estimable)) {
    h2_values <- network_results$N1$targeted
    pb <- h2_values$T_targeted[h2_values$node == "inq_perceivedburden"]
    comparator <- h2_values$T_targeted[
      h2_values$node != "inq_perceivedburden"
    ]
    h2_supported <- length(pb) == 1L && pb > max(comparator)
    h2_decision <- data.frame(
      hypothesis = "H2",
      network = "N1",
      perceived_burdensomeness_T = pb,
      largest_comparator_T = max(comparator),
      strictly_largest = h2_supported,
      decision = if (h2_supported) "supported" else "not supported",
      stringsAsFactors = FALSE
    )
  } else {
    h2_decision <- data.frame(
      hypothesis = "H2", network = "N1",
      perceived_burdensomeness_T = NA_real_, largest_comparator_T = NA_real_,
      strictly_largest = NA, decision = "N1 non-estimable",
      stringsAsFactors = FALSE
    )
  }

  # The same H2 ordering is noted for the age-range rerun only as sensitivity
  # information. It never supplies the primary H2 decision.
  age_h2_description <- data.frame(
    network = c("N1", "N1_age_sensitivity"),
    perceived_burdensomeness_uniquely_first = vapply(
      network_results[c("N1", "N1_age_sensitivity")],
      function(x) {
        if (!isTRUE(x$estimable)) return(NA)
        value <- x$targeted$T_targeted[x$targeted$node == "inq_perceivedburden"]
        others <- x$targeted$T_targeted[x$targeted$node != "inq_perceivedburden"]
        value > max(others)
      },
      logical(1)
    ),
    decision_role = c("H2 decision", "descriptive sensitivity only"),
    stringsAsFactors = FALSE
  )

  # This ordered table answers the second research question directly by
  # showing which vulnerability node has the largest N1 targeted-connectivity
  # value. The separate H2 row then checks whether that node is specifically
  # perceived burdensomeness and whether it is uniquely first.
  n1_targeted_ranking <- if (isTRUE(network_results$N1$estimable)) {
    out <- network_results$N1$targeted
    out[order(out$rank, -out$T_targeted), , drop = FALSE]
  } else {
    data.frame(note = "N1 was non-estimable", stringsAsFactors = FALSE)
  }

  list(
    h1_decisions = h1_decisions,
    h1_age_descriptive = h1_age_descriptive,
    h1_age_comparison = h1_age_comparison,
    h2_decision = h2_decision,
    age_h2_description = age_h2_description,
    n1_targeted_ranking = n1_targeted_ranking
  )
}


# Run the separate P1 classifier and evaluate H3 ---------------------------
classify_current_ideation <- function(primary_18_25) {
  # Questions 1 and 15 define P1 as concurrent classification of any
  # past-three-month ideation, not prediction of a future suicide attempt or
  # creation of an individual risk label. The base model contains the five
  # affective, interpersonal and impulsivity variables. The extended model
  # adds only the six registered sleep, nightmare and circadian variables, so
  # their comparison isolates the added information from that block.

  p1_base_predictors <- c(
    "cesd_no_sleep", "gad_total", "inq_perceivedburden",
    "inq_thwartedbelong", "upps_negativeurgency"
  )
  p1_extended_predictors <- c(
    p1_base_predictors,
    "isi_1", "isi_2", "isi_3", "isi_7", "ddnsi_total", "socialjetlag_h"
  )
  p1_models <- list(base = p1_base_predictors, extended = p1_extended_predictors)
  # No interactions, nonlinear terms, class weighting, over-sampling,
  # under-sampling, synthetic data or extra predictors are introduced because
  # none is part of the P1 plan in Question 15.
  alpha_grid <- c(0, 0.25, 0.50, 0.75, 1)

  # P1 repeats the same missing-data safeguard used by the networks. No
  # outcome or predictor is imputed; exceeding either registered 5% limit
  # stops the analysis before classifier results are examined.
  p1_data_all <- primary_18_25[c(
    "id", "si_3mo_yn", unique(unlist(p1_models, use.names = FALSE))
  )]
  p1_required_variables <- setdiff(names(p1_data_all), "id")
  p1_missing <- colMeans(is.na(p1_data_all[p1_required_variables]))
  p1_missingness <- data.frame(
    variable = p1_required_variables,
    missing_n = colSums(is.na(p1_data_all[p1_required_variables])),
    missing_proportion = as.numeric(p1_missing),
    stringsAsFactors = FALSE
  )
  p1_complete <- stats::complete.cases(p1_data_all)
  p1_case_loss <- 1 - mean(p1_complete)
  if (max(p1_missing) > 0.05 || p1_case_loss > 0.05) {
    stop("P1 exceeds the registered 5% missingness/case-loss limit. Pause and log a revision.")
  }
  p1_data <- p1_data_all[p1_complete, , drop = FALSE]
  p1_data$stratum <- factor(
    p1_data$si_3mo_yn,
    levels = c(0, 1),
    labels = c("No", "Yes")
  )
  stopifnot(
    !anyNA(p1_data$id),
    anyDuplicated(p1_data$id) == 0L,
    all(p1_data$si_3mo_yn %in% c(0L, 1L)),
    length(unique(p1_data$si_3mo_yn)) == 2L
  )

  p1_sample_audit <- data.frame(
    preregistered_expected_cohort_n = 885L,
    preregistered_expected_events = 287L,
    frozen_cohort_n = nrow(p1_data_all),
    analytic_complete_case_n = nrow(p1_data),
    frozen_cohort_observed_events = sum(p1_data_all$si_3mo_yn == 1L, na.rm = TRUE),
    analytic_events = sum(p1_data$si_3mo_yn == 1L),
    analytic_non_events = sum(p1_data$si_3mo_yn == 0L),
    event_prevalence = mean(p1_data$si_3mo_yn == 1L),
    maximum_variable_missing_proportion = max(p1_missing),
    complete_case_loss_proportion = p1_case_loss,
    complete_release_matches_expected_n = nrow(p1_data_all) == 885L,
    complete_release_matches_expected_events = if (
      anyNA(p1_data_all$si_3mo_yn)
    ) NA else sum(p1_data_all$si_3mo_yn) == 287L,
    stringsAsFactors = FALSE
  )

  # Predictors are centred and scaled from training participants only. The
  # learned means and standard deviations are then applied unchanged to the
  # corresponding validation or test participants. This prevents held-out
  # information from leaking into preprocessing (Questions 8 and 15).
  fit_scaler <- function(data, predictors) {
    center <- vapply(data[predictors], mean, numeric(1))
    scale <- vapply(data[predictors], stats::sd, numeric(1))
    if (any(!is.finite(scale)) || any(scale == 0)) {
      stop("P1 training split contains a zero-variance predictor")
    }
    list(center = center, scale = scale)
  }

  apply_scaler <- function(data, predictors, scaler) {
    x <- as.matrix(data[predictors])
    x <- sweep(x, 2, scaler$center, FUN = "-")
    x <- sweep(x, 2, scaler$scale, FUN = "/")
    storage.mode(x) <- "double"
    x
  }

  binomial_deviance_from_link <- function(y, eta) {
    # Inner tuning uses mean binomial deviance. This numerically stable form
    # calculates the same quantity directly from glmnet's linear predictor
    # without clipping probabilities; Question 15 reserves clipping for the
    # calibration logit calculation later in the script.
    softplus <- pmax(eta, 0) + log1p(exp(-abs(eta)))
    colMeans(2 * (softplus - y * eta))
  }

  # Within an outer-training set, five-fold inner cross-validation chooses
  # alpha and lambda by the lowest mean validation deviance. The untouched
  # outer-test participants do not take part in this choice. Every alpha uses
  # its registered 100-value descending path with lambda.min.ratio = 0.0001.
  tune_one_p1_model <- function(outer_train, inner_folds, predictors) {
    outer_scaler <- fit_scaler(outer_train, predictors)
    outer_x <- apply_scaler(outer_train, predictors, outer_scaler)
    outer_y <- outer_train$si_3mo_yn
    candidates <- list()
    candidate_i <- 0L

    for (alpha in alpha_grid) {
      # glmnet derives the alpha-specific path from outer-training data. The
      # explicit check prevents silent use of fewer than the registered 100
      # lambda values.
      path_fit <- glmnet::glmnet(
        x = outer_x,
        y = outer_y,
        family = "binomial",
        alpha = alpha,
        nlambda = 100,
        lambda.min.ratio = 0.0001,
        standardize = FALSE,
        intercept = TRUE
      )
      lambda_path <- path_fit$lambda
      if (length(lambda_path) != 100L) {
        stop(
          "glmnet returned ", length(lambda_path),
          " lambdas; the registered path requires exactly 100. Review before results."
        )
      }

      fold_deviance <- matrix(
        NA_real_, nrow = nrow(inner_folds), ncol = length(lambda_path)
      )
      for (inner_i in seq_len(nrow(inner_folds))) {
        inner_train <- rsample::analysis(inner_folds$splits[[inner_i]])
        inner_validation <- rsample::assessment(inner_folds$splits[[inner_i]])

        # Each inner-training fold learns its own scaling values, which are
        # then applied unchanged to that fold's validation participants.
        inner_scaler <- fit_scaler(inner_train, predictors)
        inner_x <- apply_scaler(inner_train, predictors, inner_scaler)
        validation_x <- apply_scaler(inner_validation, predictors, inner_scaler)

        fit <- glmnet::glmnet(
          x = inner_x,
          y = inner_train$si_3mo_yn,
          family = "binomial",
          alpha = alpha,
          lambda = lambda_path,
          standardize = FALSE,
          intercept = TRUE
        )
        eta <- predict(
          fit,
          newx = validation_x,
          s = lambda_path,
          type = "link",
          exact = FALSE
        )
        fold_deviance[inner_i, ] <- binomial_deviance_from_link(
          inner_validation$si_3mo_yn,
          eta
        )
      }

      mean_deviance <- colMeans(fold_deviance)
      candidate_i <- candidate_i + 1L
      candidates[[candidate_i]] <- data.frame(
        alpha = alpha,
        lambda = lambda_path,
        lambda_index_descending = seq_along(lambda_path),
        mean_inner_binomial_deviance = mean_deviance,
        stringsAsFactors = FALSE
      )
    }

    candidate_table <- do.call(rbind, candidates)
    minimum <- min(candidate_table$mean_inner_binomial_deviance)
    tied <- which(candidate_table$mean_inner_binomial_deviance == minimum)
    # Exact ties are resolved by the frozen alpha order, then descending lambda.
    # The count is retained so this unregistered mechanical rule is auditable.
    selected <- candidate_table[tied[1], , drop = FALSE]
    selected$exact_minimum_tie_count <- length(tied)
    selected_lambda_path <- candidate_table$lambda[
      candidate_table$alpha == selected$alpha
    ]

    list(
      selected = selected,
      selected_lambda_path = selected_lambda_path,
      scaler = outer_scaler,
      x = outer_x,
      y = outer_y
    )
  }

  # glmnet can normally stop a path early after a small deviance gain. Turning
  # that shortcut off preserves all 100 registered values; the original global
  # controls are restored automatically when this function finishes.
  previous_glmnet_control <- glmnet::glmnet.control()
  glmnet::glmnet.control(fdev = 0, devmax = 1)
  on.exit(
    do.call(glmnet::glmnet.control, previous_glmnet_control),
    add = TRUE
  )

  # A stage-specific seed keeps P1 folds stable even if an earlier network has
  # a different number of failed bootstrap attempts on another system.
  p1_seed <- as.integer((as.double(ANALYSIS_SEED) + 100000) %% 2147483647)
  if (p1_seed == 0L) p1_seed <- 1L
  set.seed(p1_seed)
  # Question 15 specifies five repetitions of stratified 10-fold outer
  # cross-validation, giving 50 held-out evaluations. Stratification keeps the
  # No/Yes balance as similar as possible across folds, while both models use
  # the same splits for a fair paired comparison.
  outer_folds <- rsample::vfold_cv(
    p1_data,
    v = 10,
    repeats = 5,
    strata = stratum
  )
  stopifnot(nrow(outer_folds) == 50L)

  prediction_rows <- list()
  tuning_rows <- list()
  fold_rows <- list()
  prediction_i <- tuning_i <- fold_i <- 0L

  for (outer_i in seq_len(nrow(outer_folds))) {
    outer_train <- rsample::analysis(outer_folds$splits[[outer_i]])
    outer_test <- rsample::assessment(outer_folds$splits[[outer_i]])
    repeat_label <- as.character(outer_folds$id[outer_i])
    fold_label <- as.character(outer_folds$id2[outer_i])

    # The two models also share one inner split object inside each outer split.
    # The form explicitly requires identical outer folds; pairing the inner
    # folds additionally avoids needless random imbalance during tuning.
    set.seed(stage_seed(100000 + outer_i))
    inner_folds <- rsample::vfold_cv(
      outer_train,
      v = 5,
      strata = stratum
    )

    fold_i <- fold_i + 1L
    fold_rows[[fold_i]] <- data.frame(
      id = as.character(outer_test$id),
      repetition = repeat_label,
      outer_fold = fold_label,
      stringsAsFactors = FALSE
    )

    for (model_name in names(p1_models)) {
      predictors <- p1_models[[model_name]]
      tuned <- tune_one_p1_model(outer_train, inner_folds, predictors)
      selected <- tuned$selected

      # After tuning, the chosen model is fitted from all outer-training
      # participants. Their scaling values are applied to the untouched test
      # fold, which receives one genuinely held-out probability.
      test_x <- apply_scaler(outer_test, predictors, tuned$scaler)
      final_fit <- glmnet::glmnet(
        x = tuned$x,
        y = tuned$y,
        family = "binomial",
        alpha = selected$alpha,
        lambda = tuned$selected_lambda_path,
        standardize = FALSE,
        intercept = TRUE
      )
      probability <- as.numeric(predict(
        final_fit,
        newx = test_x,
        s = selected$lambda,
        type = "response",
        exact = FALSE
      ))
      if (any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
        stop("P1 produced an invalid outer-test probability")
      }

      prediction_i <- prediction_i + 1L
      prediction_rows[[prediction_i]] <- data.frame(
        id = as.character(outer_test$id),
        outcome = outer_test$si_3mo_yn,
        repetition = repeat_label,
        outer_fold = fold_label,
        model = model_name,
        probability = probability,
        stringsAsFactors = FALSE
      )

      tuning_i <- tuning_i + 1L
      tuning_rows[[tuning_i]] <- data.frame(
        repetition = repeat_label,
        outer_fold = fold_label,
        model = model_name,
        outer_training_n = nrow(outer_train),
        outer_test_n = nrow(outer_test),
        alpha = selected$alpha,
        lambda = selected$lambda,
        minimum_inner_binomial_deviance = selected$mean_inner_binomial_deviance,
        exact_minimum_tie_count = selected$exact_minimum_tie_count,
        stringsAsFactors = FALSE
      )
    }
  }

  p1_predictions <- do.call(rbind, prediction_rows)
  p1_tuning <- do.call(rbind, tuning_rows)
  p1_fold_manifest <- do.call(rbind, fold_rows)

  # Every participant should receive one held-out prediction in each of the
  # five repetitions from each model. These checks confirm that none is
  # duplicated or omitted and that base and extended predictions refer to the
  # same participant-repetition pairs.
  prediction_counts <- as.data.frame(table(
    p1_predictions$id,
    p1_predictions$model
  ), stringsAsFactors = FALSE)
  names(prediction_counts) <- c("id", "model", "n_predictions")
  stopifnot(
    nrow(p1_predictions) == nrow(p1_data) * 5L * 2L,
    all(prediction_counts$n_predictions == 5L),
    anyDuplicated(p1_predictions[c("id", "repetition", "model")]) == 0L
  )
  base_keys <- p1_predictions[
    p1_predictions$model == "base", c("id", "repetition")
  ]
  extended_keys <- p1_predictions[
    p1_predictions$model == "extended", c("id", "repetition")
  ]
  base_key_text <- sort(paste(
    base_keys$id, base_keys$repetition, sep = "\r"
  ))
  extended_key_text <- sort(paste(
    extended_keys$id, extended_keys$repetition, sep = "\r"
  ))
  stopifnot(identical(base_key_text, extended_key_text))

  p1_metric_names <- c(
    "brier", "roc_auc", "pr_auc",
    "calibration_intercept", "calibration_slope"
  )
  p1_calibration_warning_log <- character()

  # Within each outer repetition, predictions from its ten held-out folds are
  # pooled before scoring. Brier score is primary and is better when lower
  # because it averages squared probability error. ROC-AUC, PR-AUC,
  # calibration intercept and calibration slope are secondary, as stated in
  # Question 15.
  p1_metric_bundle <- function(data) {
    y <- data$outcome
    p <- data$probability
    if (length(unique(y)) != 2L) stop("P1 metric sample does not contain both outcomes")

    # Probabilities are clipped only here to avoid infinite logits at exactly
    # zero or one. Their original values remain unchanged for Brier score and
    # both AUC measures. The offset model estimates calibration-in-the-large;
    # the second model estimates the calibration slope described in the form.
    p_calibration <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    z <- stats::qlogis(p_calibration)
    calibration_data <- data.frame(y = y, z = z)
    calibration_warnings <- character()
    calibration_fits <- withCallingHandlers(
      list(
        intercept = stats::glm(
          y ~ 1 + offset(z),
          data = calibration_data,
          family = stats::binomial()
        ),
        slope = stats::glm(
          y ~ z,
          data = calibration_data,
          family = stats::binomial()
        )
      ),
      warning = function(w) {
        calibration_warnings <<- c(
          calibration_warnings,
          conditionMessage(w)
        )
        invokeRestart("muffleWarning")
      }
    )
    calibration_coefficients <- c(
      intercept = unname(stats::coef(calibration_fits$intercept)[1]),
      slope = unname(stats::coef(calibration_fits$slope)["z"])
    )
    if (length(calibration_warnings) > 0L) {
      p1_calibration_warning_log <<- c(
        p1_calibration_warning_log,
        calibration_warnings
      )
    }
    # A fitted-probability warning can occur in a bootstrap sample even when
    # both calibration models converge and return finite coefficients. Because
    # calibration is secondary and H3 is decided by Brier score, that warning
    # is recorded rather than allowed to erase an otherwise valid bootstrap.
    # A non-converged or non-finite calibration fit still stops P1.
    if (
      !isTRUE(calibration_fits$intercept$converged) ||
        !isTRUE(calibration_fits$slope$converged) ||
        any(!is.finite(calibration_coefficients))
    ) {
      stop(
        "P1 calibration failed to converge or returned a non-finite value. ",
        paste(unique(calibration_warnings), collapse = " | ")
      )
    }

    values <- c(
      brier = mean((y - p)^2),
      roc_auc = as.numeric(pROC::auc(pROC::roc(
        response = y,
        predictor = p,
        levels = c(0, 1),
        direction = "<",
        quiet = TRUE
      ))),
      pr_auc = PRROC::pr.curve(
        scores.class0 = p[y == 1L],
        scores.class1 = p[y == 0L],
        curve = FALSE
      )$auc.integral,
      calibration_intercept = calibration_coefficients["intercept"],
      calibration_slope = calibration_coefficients["slope"]
    )
    names(values)[names(values) == "calibration_intercept.intercept"] <-
      "calibration_intercept"
    names(values)[names(values) == "calibration_slope.slope"] <-
      "calibration_slope"
    if (any(!is.finite(values))) stop("P1 produced a non-finite registered metric")
    if (!identical(names(values), p1_metric_names)) {
      stop("P1 metric names do not match the registered metric set")
    }
    values
  }

  # Metrics are first calculated separately in each of the five repetitions,
  # after pooling the ten test folds, and are then averaged across repetitions.
  p1_group_keys <- unique(p1_predictions[c("model", "repetition")])
  p1_repetition_metrics <- do.call(rbind, lapply(
    seq_len(nrow(p1_group_keys)),
    function(i) {
      key <- p1_group_keys[i, ]
      subset <- p1_predictions[
        p1_predictions$model == key$model &
          p1_predictions$repetition == key$repetition,
        , drop = FALSE
      ]
      values <- p1_metric_bundle(subset)
      data.frame(
        model = key$model,
        repetition = key$repetition,
        as.list(values),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  ))

  p1_mean_metrics <- do.call(rbind, lapply(names(p1_models), function(model_name) {
    rows <- p1_repetition_metrics[p1_repetition_metrics$model == model_name, ]
    data.frame(
      model = model_name,
      as.list(colMeans(rows[p1_metric_names])),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  stopifnot(nrow(p1_mean_metrics) == 2L)

  base_mean <- p1_mean_metrics[p1_mean_metrics$model == "base", p1_metric_names]
  extended_mean <- p1_mean_metrics[
    p1_mean_metrics$model == "extended", p1_metric_names
  ]
  p1_paired_differences <- data.frame(
    metric = p1_metric_names,
    contrast = "extended_minus_base",
    estimate = as.numeric(unlist(extended_mean - base_mean, use.names = FALSE)),
    stringsAsFactors = FALSE
  )

  # Question 15 requires 1,000 outcome-stratified participant bootstraps for
  # paired 95% intervals. Sampling a participant carries all five repeated
  # predictions and both models together, preserving the repeated and paired
  # structure instead of treating prediction rows as independent observations.
  p1_boot_seed <- as.integer((as.double(ANALYSIS_SEED) + 200000) %% 2147483647)
  if (p1_boot_seed == 0L) p1_boot_seed <- 1L
  set.seed(p1_boot_seed)

  id_outcome <- unique(p1_predictions[c("id", "outcome")])
  id_outcome <- id_outcome[order(id_outcome$id), , drop = FALSE]
  ids_by_outcome <- split(id_outcome$id, id_outcome$outcome)
  prediction_indices <- split(seq_len(nrow(p1_predictions)), p1_predictions$id)
  p1_boot_rows <- vector("list", N_P1_BOOT)

  for (bootstrap_i in seq_len(N_P1_BOOT)) {
    sampled_ids <- unlist(lapply(ids_by_outcome, function(ids) {
      sample(ids, size = length(ids), replace = TRUE)
    }), use.names = FALSE)
    sampled_rows <- unlist(prediction_indices[sampled_ids], use.names = FALSE)
    draw <- p1_predictions[sampled_rows, , drop = FALSE]

    averaged <- lapply(names(p1_models), function(model_name) {
      model_draw <- draw[draw$model == model_name, , drop = FALSE]
      repeat_values <- lapply(
        split(model_draw, model_draw$repetition),
        p1_metric_bundle
      )
      Reduce("+", repeat_values) / length(repeat_values)
    })
    names(averaged) <- names(p1_models)
    difference <- averaged$extended - averaged$base

    p1_boot_rows[[bootstrap_i]] <- rbind(
      data.frame(
        bootstrap = bootstrap_i, scope = "base",
        metric = names(averaged$base), value = as.numeric(averaged$base)
      ),
      data.frame(
        bootstrap = bootstrap_i, scope = "extended",
        metric = names(averaged$extended), value = as.numeric(averaged$extended)
      ),
      data.frame(
        bootstrap = bootstrap_i, scope = "extended_minus_base",
        metric = names(difference), value = as.numeric(difference)
      )
    )
  }
  p1_bootstrap_values <- do.call(rbind, p1_boot_rows)

  p1_bootstrap_ci <- do.call(rbind, lapply(
    split(
      p1_bootstrap_values,
      interaction(
        p1_bootstrap_values$scope,
        p1_bootstrap_values$metric,
        drop = TRUE,
        lex.order = TRUE
      )
    ),
    function(x) {
      interval <- stats::quantile(
        x$value,
        probs = c(0.025, 0.975),
        names = FALSE,
        type = 7
      )
      data.frame(
        scope = x$scope[1],
        metric = x$metric[1],
        estimate = if (x$scope[1] == "extended_minus_base") {
          p1_paired_differences$estimate[
            p1_paired_differences$metric == x$metric[1]
          ]
        } else {
          p1_mean_metrics[
            p1_mean_metrics$model == x$scope[1], x$metric[1]
          ]
        },
        percentile_lower = interval[1],
        percentile_upper = interval[2],
        bootstrap_repetitions = nrow(x),
        stringsAsFactors = FALSE
      )
    }
  ))

  # H3 defines the contrast as extended Brier minus base Brier. A negative
  # point estimate is directionally consistent with H3; an interval wholly
  # below zero is stronger internal evidence. Neither result is external
  # validation or prospective individual-risk prediction (Question 16).
  h3_row <- p1_bootstrap_ci[
    p1_bootstrap_ci$scope == "extended_minus_base" &
      p1_bootstrap_ci$metric == "brier",
    , drop = FALSE
  ]
  h3_decision <- data.frame(
    hypothesis = "H3",
    contrast = "mean Brier: extended minus base",
    estimate = h3_row$estimate,
    percentile_lower = h3_row$percentile_lower,
    percentile_upper = h3_row$percentile_upper,
    directionally_consistent_with_H3 = h3_row$estimate < 0,
    interval_wholly_below_zero_stronger_evidence = h3_row$percentile_upper < 0,
    qualification = "internal concurrent classification; not external validation or prospective risk prediction",
    stringsAsFactors = FALSE
  )

  # Keep any benign calibration warnings visible in the ordinary P1 sample
  # summary instead of creating another output file just for warning text.
  p1_sample_audit$calibration_warning_count <-
    length(p1_calibration_warning_log)
  p1_sample_audit$calibration_warning_types <- if (
    length(p1_calibration_warning_log) == 0L
  ) {
    "none"
  } else {
    paste(unique(p1_calibration_warning_log), collapse = " | ")
  }

  list(
    sample_audit = p1_sample_audit,
    missingness = p1_missingness,
    fold_manifest = p1_fold_manifest,
    tuning = p1_tuning,
    predictions = p1_predictions,
    repetition_metrics = p1_repetition_metrics,
    mean_metrics = p1_mean_metrics,
    paired_differences = p1_paired_differences,
    bootstrap_values = p1_bootstrap_values,
    bootstrap_intervals = p1_bootstrap_ci,
    h3_decision = h3_decision,
    p1_seed = p1_seed,
    p1_boot_seed = p1_boot_seed
  )
}


# Describe the rare D1 outcomes --------------------------------------------
describe_rare_outcomes <- function(primary_18_25) {
  # Question 1 asks only for separate frequencies and percentages of plan with
  # intent, preparatory acts, non-suicidal self-injury and suicide attempt.
  # These indicators can overlap within one participant and are not combined
  # into a composite. Their small event counts are why they are not used as
  # network nodes or multivariable-classifier outcomes.
  d1_variables <- c(
    cssrs_3mo_5 = "Suicide plan with intent in past 3 months",
    cssrs_3mo_6 = "Preparatory acts in past 3 months",
    cssrs_3mo_7 = "Non-suicidal self-injury in past 3 months",
    cssrs_3mo_8 = "Suicide attempt in past 3 months"
  )
  d1_counts <- do.call(rbind, lapply(names(d1_variables), function(variable) {
    endorsed <- sum(primary_18_25[[variable]] == 1L)
    data.frame(
      variable = variable,
      clinical_label = unname(d1_variables[variable]),
      endorsed_n = endorsed,
      denominator = nrow(primary_18_25),
      percentage = 100 * endorsed / nrow(primary_18_25),
      non_mutually_exclusive = TRUE,
      stringsAsFactors = FALSE
    )
  }))
  stopifnot(identical(d1_counts$endorsed_n, c(13L, 20L, 93L, 16L)))
  d1_counts
}


# Organize the answers and reproducibility record --------------------------
save_analysis_results <- function(
  analysis_setup, network_analysis, network_answers,
  classification_results, rare_outcome_counts
) {
  dirs <- analysis_setup$dirs
  private_result_path <- analysis_setup$private_result_path
  registered_packages <- analysis_setup$registered_packages
  network_results <- network_analysis$network_results
  network_status <- network_analysis$network_status

  # The workbook keeps the analysis record, package versions, disclosed choices, and session information together.
  session_info_lines <- capture.output(utils::sessionInfo())
  run_information <- rbind(
    data.frame(section = "analysis record", item = analysis_setup$analysis_record$item,
               value = as.character(analysis_setup$analysis_record$value), note = "", stringsAsFactors = FALSE),
    data.frame(section = "package version", item = analysis_setup$package_versions$package,
               value = analysis_setup$package_versions$version, note = "", stringsAsFactors = FALSE),
    data.frame(section = "choice not fixed numerically in the form",
               item = analysis_setup$implementation_choices$item,
               value = as.character(analysis_setup$implementation_choices$value),
               note = analysis_setup$implementation_choices$pdf_status, stringsAsFactors = FALSE),
    data.frame(section = "session information", item = paste0("line_", seq_along(session_info_lines)),
               value = session_info_lines, note = "", stringsAsFactors = FALSE)
  )
  deviation_log <- analysis_setup$deviation_log

  h1_decisions <- network_answers$h1_decisions
  h1_age_descriptive <- network_answers$h1_age_descriptive
  h1_age_comparison <- network_answers$h1_age_comparison
  h2_decision <- network_answers$h2_decision
  age_h2_description <- network_answers$age_h2_description

  p1_sample_audit <- classification_results$sample_audit
  p1_missingness <- classification_results$missingness
  p1_fold_manifest <- classification_results$fold_manifest
  p1_tuning <- classification_results$tuning
  p1_predictions <- classification_results$predictions
  p1_repetition_metrics <- classification_results$repetition_metrics
  p1_mean_metrics <- classification_results$mean_metrics
  p1_paired_differences <- classification_results$paired_differences
  p1_bootstrap_values <- classification_results$bootstrap_values
  p1_bootstrap_ci <- classification_results$bootstrap_intervals
  h3_decision <- classification_results$h3_decision
  p1_seed <- classification_results$p1_seed
  p1_boot_seed <- classification_results$p1_boot_seed

  # The workbook collects aggregate, non-sensitive outputs in named sheets so
  # that primary decisions, sensitivity results, uncertainty and diagnostics
  # remain distinguishable. Questions 15 and 29 also require the seed, package
  # versions, session information, renv lockfile and deviation log. Together,
  # these records make the run reproducible and expose any later departure
  # from the preregistration.

  bind_network_rows <- function(extractor, include_empty = FALSE) {
    pieces <- lapply(network_results, extractor)
    pieces <- pieces[!vapply(pieces, is.null, logical(1))]
    if (!include_empty) {
      pieces <- pieces[vapply(pieces, nrow, integer(1)) > 0L]
    }
    if (length(pieces) == 0L) {
      return(data.frame(note = "No rows available", stringsAsFactors = FALSE))
    }
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    out
  }

  network_missingness_table <- bind_network_rows(function(x) {
    out <- x$missingness$table
    out$complete_case_loss_proportion <- x$missingness$case_loss
    out
  })
  node_audit_table <- bind_network_rows(function(x) x$audit$nodes)
  frequency_audit_table <- bind_network_rows(function(x) x$audit$frequencies)
  redundancy_audit_table <- bind_network_rows(function(x) x$audit$redundancy)
  matrix_diagnostic_table <- bind_network_rows(function(x) x$fit$diagnostic)
  near_pd_table <- bind_network_rows(function(x) x$fit$near_pd)
  all_edges <- bind_network_rows(function(x) x$edges)
  all_edge_intervals <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$boot$edge_intervals else NULL
  })
  all_targeted <- bind_network_rows(function(x) x$targeted)
  all_centrality <- bind_network_rows(function(x) x$centrality)
  all_bridge <- bind_network_rows(function(x) x$bridge)
  all_retry_counts <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$boot$retries else NULL
  })
  bootstrap_warning_pieces <- lapply(network_results, function(x) {
    if (isTRUE(x$estimable) && nrow(x$boot$warnings) > 0L) x$boot$warnings else NULL
  })
  bootstrap_warning_pieces <- bootstrap_warning_pieces[
    !vapply(bootstrap_warning_pieces, is.null, logical(1))
  ]
  all_bootstrap_warnings <- if (length(bootstrap_warning_pieces)) {
    do.call(rbind, bootstrap_warning_pieces)
  } else {
    data.frame(
      network = character(), bootstrap = character(), warning = character(),
      stringsAsFactors = FALSE
    )
  }
  all_cs <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$boot$cs else NULL
  })
  all_custom_component_summary <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$custom$component_summary else NULL
  })
  all_custom_bridge_summary <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$custom$bridge_summary else NULL
  })
  all_custom_rank_first <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$custom$rank_first_summary else NULL
  })
  all_custom_rank_correlations <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$custom$rank_correlation_summary else NULL
  })
  all_custom_rank_distribution <- bind_network_rows(function(x) {
    if (!isTRUE(x$estimable)) return(NULL)
    draws <- x$custom$components[
      x$custom$components$component == "T_targeted",
      ,
      drop = FALSE
    ]
    if (nrow(draws) == 0L) return(NULL)

    # Reporting the full rank distribution as well as first-place frequency
    # keeps near-ties visible, as requested in Questions 12, 16 and 23. Minimum
    # ranks implement the disclosed table tie rule.
    counts <- stats::aggregate(
      rep(1L, nrow(draws)),
      by = list(
        network = draws$network,
        drop_fraction = draws$drop_fraction,
        node = draws$node,
        resampled_rank = draws$rank
      ),
      FUN = sum
    )
    names(counts)[ncol(counts)] <- "repetitions_at_rank"
    denominators <- stats::aggregate(
      replicate ~ network + drop_fraction + node,
      data = draws,
      FUN = function(z) length(unique(z))
    )
    names(denominators)[ncol(denominators)] <- "successful_network_fits"
    out <- merge(
      counts,
      denominators,
      by = c("network", "drop_fraction", "node"),
      sort = FALSE
    )
    out$rank_proportion <-
      out$repetitions_at_rank / out$successful_network_fits
    full_rank <- x$targeted[c("node", "rank")]
    names(full_rank)[2] <- "full_sample_rank"
    merge(out, full_rank, by = "node", all.x = TRUE, sort = FALSE)
  })
  all_custom_success <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$custom$success_counts else NULL
  })
  all_custom_failures <- bind_network_rows(function(x) {
    if (isTRUE(x$estimable)) x$custom$failures else NULL
  })
  custom_warning_pieces <- lapply(network_results, function(x) {
    if (isTRUE(x$estimable) && nrow(x$custom$warnings) > 0L) x$custom$warnings else NULL
  })
  custom_warning_pieces <- custom_warning_pieces[
    !vapply(custom_warning_pieces, is.null, logical(1))
  ]
  all_custom_warnings <- if (length(custom_warning_pieces)) {
    do.call(rbind, custom_warning_pieces)
  } else {
    data.frame(
      network = character(), drop_fraction = numeric(), replicate = integer(),
      warning = character(), stringsAsFactors = FALSE
    )
  }

  full_warning_pieces <- lapply(network_results, function(x) {
    if (length(x$fit$warnings) == 0L) return(NULL)
    data.frame(
      network = x$spec$name,
      warning = unique(x$fit$warnings),
      stringsAsFactors = FALSE
    )
  })
  full_warning_pieces <- full_warning_pieces[
    !vapply(full_warning_pieces, is.null, logical(1))
  ]
  full_network_warnings <- if (length(full_warning_pieces)) {
    do.call(rbind, full_warning_pieces)
  } else {
    data.frame(network = character(), warning = character(), stringsAsFactors = FALSE)
  }

  warning_events <- rbind(
    data.frame(
      source = rep("full_sample", nrow(full_network_warnings)),
      network = full_network_warnings$network,
      context = rep("full_sample", nrow(full_network_warnings)),
      drop_fraction = rep("", nrow(full_network_warnings)),
      warning = full_network_warnings$warning, stringsAsFactors = FALSE
    ),
    data.frame(
      source = rep("bootnet", nrow(all_bootstrap_warnings)),
      network = all_bootstrap_warnings$network,
      context = all_bootstrap_warnings$bootstrap,
      drop_fraction = rep("", nrow(all_bootstrap_warnings)),
      warning = all_bootstrap_warnings$warning, stringsAsFactors = FALSE
    ),
    data.frame(
      source = rep("custom_case_drop", nrow(all_custom_warnings)),
      network = all_custom_warnings$network,
      context = rep("custom_case_drop", nrow(all_custom_warnings)),
      drop_fraction = as.character(all_custom_warnings$drop_fraction),
      warning = all_custom_warnings$warning, stringsAsFactors = FALSE
    )
  )
  warning_counts <- if (nrow(warning_events) > 0L) {
    out <- stats::aggregate(
      rep(1L, nrow(warning_events)),
      by = warning_events,
      FUN = sum
    )
    names(out)[ncol(out)] <- "occurrences"
    out
  } else {
    data.frame(note = "No estimation warnings", stringsAsFactors = FALSE)
  }

  # This compact table is the first place to read after the run. It applies
  # only the decision rules written in the preregistration: positive N1 point
  # estimates for H1a-H1e, a uniquely largest perceived-burdensomeness T for
  # H2, and an extended-minus-base Brier difference below zero for H3.
  h1_answer_rows <- data.frame(
    research_question = rep("RQ1", nrow(h1_decisions)),
    hypothesis = h1_decisions$hypothesis,
    analysis = rep("N1", nrow(h1_decisions)),
    quantity = paste(
      h1_decisions$sleep_endpoint,
      "-- si_3mo_maxseverity edge"
    ),
    estimate = h1_decisions$estimate,
    comparison_value = rep(0, nrow(h1_decisions)),
    percentile_lower = h1_decisions$percentile_lower,
    percentile_upper = h1_decisions$percentile_upper,
    prespecified_rule = rep(
      "full-sample regularised N1 edge > 0",
      nrow(h1_decisions)
    ),
    answer = ifelse(
      is.na(h1_decisions$directionally_consistent),
      "N1 non-estimable",
      ifelse(
        h1_decisions$directionally_consistent,
        "directionally consistent with H1",
        "not directionally consistent with H1"
      )
    ),
    decision_role = rep("primary network decision", nrow(h1_decisions)),
    stringsAsFactors = FALSE
  )
  h2_answer_row <- data.frame(
    research_question = "RQ2",
    hypothesis = "H2",
    analysis = "N1",
    quantity = "perceived-burdensomeness targeted connectivity",
    estimate = h2_decision$perceived_burdensomeness_T,
    comparison_value = h2_decision$largest_comparator_T,
    percentile_lower = NA_real_,
    percentile_upper = NA_real_,
    prespecified_rule = "PB T strictly greater than every comparator T",
    answer = h2_decision$decision,
    decision_role = "primary network decision",
    stringsAsFactors = FALSE
  )
  h3_answer_row <- data.frame(
    research_question = "Supplementary P1",
    hypothesis = "H3",
    analysis = "P1",
    quantity = "mean Brier difference: extended minus base",
    estimate = h3_decision$estimate,
    comparison_value = 0,
    percentile_lower = h3_decision$percentile_lower,
    percentile_upper = h3_decision$percentile_upper,
    prespecified_rule = "extended-minus-base mean Brier < 0",
    answer = ifelse(
      h3_decision$directionally_consistent_with_H3,
      "directionally consistent with H3",
      "not directionally consistent with H3"
    ),
    decision_role = "supplementary classifier decision",
    stringsAsFactors = FALSE
  )
  prespecified_answers <- rbind(
    h1_answer_rows,
    h2_answer_row,
    h3_answer_row
  )

  aggregate_workbook <- list(
    Prespecified_answers = prespecified_answers,
    Run_information = run_information,
    Deviation_log = deviation_log,
    RQ2_N1_ranking = network_answers$n1_targeted_ranking,
    Network_status = network_status,
    Missingness = network_missingness_table,
    Node_audit = node_audit_table,
    Category_frequencies = frequency_audit_table,
    Redundancy_audit = redundancy_audit_table,
    Matrix_diagnostics = matrix_diagnostic_table,
    nearPD_diagnostic_only = near_pd_table,
    Edges = all_edges,
    Edge_percentile_intervals = all_edge_intervals,
    Targeted_connectivity = all_targeted,
    Centrality = all_centrality,
    Bridge_EI_1step_normalized = all_bridge,
    Bootstrap_fit_counts = all_retry_counts,
    Recorded_warning_counts = warning_counts,
    CS_coefficients = all_cs,
    Custom_component_summary = all_custom_component_summary,
    Custom_bridge_summary = all_custom_bridge_summary,
    Custom_rank_first = all_custom_rank_first,
    Custom_rank_correlation = all_custom_rank_correlations,
    Custom_rank_distribution = all_custom_rank_distribution,
    Custom_fit_counts = all_custom_success,
    Custom_failure_reasons = all_custom_failures,
    H1_decisions = h1_decisions,
    H1_age_comparison = h1_age_comparison,
    H2_decision = h2_decision,
    H2_age_description = age_h2_description,
    P1_sample = p1_sample_audit,
    P1_missingness = p1_missingness,
    P1_tuning = p1_tuning,
    P1_repeat_metrics = p1_repetition_metrics,
    P1_mean_metrics = p1_mean_metrics,
    P1_differences = p1_paired_differences,
    P1_bootstrap_CI = p1_bootstrap_ci,
    H3_decision = h3_decision,
    D1_counts = rare_outcome_counts
  )
  # Participant IDs, fold assignments, individual probabilities and raw
  # resample draws remain in the private RDS. They are deliberately excluded
  # from the aggregate reporting workbook, consistent with the restricted-data
  # statement in Question 29.
  result_object <- list(
    configuration = list(
      seed = ANALYSIS_SEED,
      p1_seed = p1_seed,
      p1_bootstrap_seed = p1_boot_seed,
      edge_interval_level = EDGE_CI_LEVEL,
      n_edge_boot = N_EDGE_BOOT,
      n_case_boot = N_CASE_BOOT,
      n_custom_boot_per_fraction = N_CUSTOM_BOOT,
      custom_drop_fractions = CUSTOM_DROP_FRACTIONS,
      n_p1_boot = N_P1_BOOT,
      n_cores = N_CORES
    ),
    network_status = network_status,
    networks = network_results,
    hypotheses = list(
      prespecified_answers = prespecified_answers,
      H1 = h1_decisions,
      H1_age_descriptive = h1_age_descriptive,
      H2 = h2_decision,
      H2_age_descriptive = age_h2_description,
      RQ2_N1_ranking = network_answers$n1_targeted_ranking,
      H3 = h3_decision
    ),
    p1 = list(
      sample_audit = p1_sample_audit,
      missingness = p1_missingness,
      fold_manifest = p1_fold_manifest,
      tuning = p1_tuning,
      predictions = p1_predictions,
      repetition_metrics = p1_repetition_metrics,
      mean_metrics = p1_mean_metrics,
      paired_differences = p1_paired_differences,
      bootstrap_values = p1_bootstrap_values,
      bootstrap_intervals = p1_bootstrap_ci
    ),
    d1 = rare_outcome_counts,
    diagnostic_tables = list(
      full_sample_warnings = full_network_warnings,
      bootstrap_warnings = all_bootstrap_warnings,
      custom_warnings = all_custom_warnings,
      custom_failures = all_custom_failures
    ),
    output_files = list(
      prepared_datasets = output_file,
      results_workbook = file.path(dirs$root, "ANSWERS_registered_results.xlsx"),
      figures = network_analysis$figure_path,
      original_figures = network_analysis$original_figure_path,
      private_rds = private_result_path,
      session_info = file.path(dirs$root, "sessionInfo.txt"),
      renv_lock = file.path(dirs$root, "renv.lock")
    )
  )
  writexl::write_xlsx(
    aggregate_workbook,
    path = result_object$output_files$results_workbook,
    format_headers = TRUE
  )
  writeLines(session_info_lines, con = result_object$output_files$session_info)
  renv::snapshot(
    project = PROJECT_DIR,
    lockfile = result_object$output_files$renv_lock,
    packages = registered_packages,
    prompt = FALSE
  )
  saveRDS(result_object, private_result_path, compress = "xz")

  message("Registered analysis saved under: ", dirs$root)
  message("Prespecified answers (also in the first workbook sheet):")
  print(
    prespecified_answers[c(
      "research_question", "hypothesis", "quantity", "estimate", "answer"
    )],
    row.names = FALSE
  )
  invisible(result_object)
}

# =============================================================================
# 04 - RUN THE COMPLETE REGISTERED ANALYSIS FROM RAW DATA
# =============================================================================

# Separate centrality panels prevent nearly overlapping strength and expected-influence traces from hiding one another.
plot_centrality_page <- function(network) {
  centrality <- network$centrality
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(8, 4, 4, 1), oma = c(2, 0, 0, 0))
  graphics::barplot(centrality$strength, names.arg = centrality$node, las = 2, col = "#2A9D8F", border = NA, ylab = "Sum of absolute edge weights", main = paste(network$spec$name, "Strength"))
  graphics::barplot(centrality$expected_influence, names.arg = centrality$node, las = 2, col = "#E76F51", border = NA, ylab = "Sum of signed edge weights", main = paste(network$spec$name, "Expected influence"))
  graphics::mtext("Separate panels keep the nearly overlapping measures visible.", side = 1, outer = TRUE, line = 0.5, cex = 0.8)
}

plot_bridge_page <- function(network) {
  bridge <- network$bridge
  values <- bridge$normalized_bridge_expected_influence_1step
  ordering <- order(values)
  communities <- factor(bridge$community, levels = unique(bridge$community))
  community_colours <- grDevices::hcl.colors(nlevels(communities), palette = "Dark 3")
  colours <- community_colours[as.integer(communities)]
  # dotchart uses the numeric codes and labels of a factor to separate groups.
  graphics::dotchart(values[ordering], labels = bridge$node[ordering], groups = communities[ordering], color = colours[ordering], pch = 19, xlab = "Normalized one-step bridge expected influence", main = paste(network$spec$name, "Bridge expected influence"))
  graphics::abline(v = 0, lty = 2, col = "grey50")
}

write_detailed_figures <- function(network_results, path) {
  grDevices::pdf(path, width = 11, height = 8.5, useDingbats = FALSE)
  device_open <- TRUE
  on.exit(if (device_open) grDevices::dev.off(), add = TRUE)
  for (network in network_results) {
    plot_network(network$fit$graph, network$spec)
    draw_edge_accuracy_plot(network$boot$edge_intervals, network$spec)
    draw_case_bootnet_plot(network$boot$case)
    plot_centrality_page(network)
    plot_bridge_page(network)
  }
  grDevices::dev.off()
  device_open <- FALSE
  invisible(path)
}

run_complete_analysis <- function() {
  message("Preparing the registered datasets from the raw CSV...")
  prepared <- prepare_answers_data(INPUT_FILE, output_file, write_output = TRUE)
  analysis_setup <- record_analysis_setup(prepared$cohorts)
  specifications <- define_registered_networks()

  message("Fitting all six registered networks and resampling each from scratch...")
  network_analysis <- estimate_registered_networks(specifications, prepared$cohorts, analysis_setup)
  network_analysis$original_figure_path <- network_analysis$figure_path
  network_answers <- answer_network_questions(network_analysis$network_results)

  message("Running the registered P1 classifier and D1 descriptions...")
  classification <- classify_current_ideation(prepared$cohorts$primary_18_25)
  d1 <- describe_rare_outcomes(prepared$cohorts$primary_18_25)

  detailed_figures <- file.path(OUTPUT_DIR, "ANSWERS_registered_figures_detailed.pdf")
  write_detailed_figures(network_analysis$network_results, detailed_figures)
  network_analysis$figure_path <- detailed_figures

  result <- save_analysis_results(analysis_setup, network_analysis, network_answers, classification, d1)
  message("Complete analysis saved under: ", OUTPUT_DIR)
  result
}

