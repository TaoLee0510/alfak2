# Scoreboard aggregation.
#
# Datasets are averaged within landscape and scenario first, then landscapes
# are weighted equally. Paired comparisons resample landscapes, the independent
# unit; patients and scenarios of one landscape are not independent replicates.

stt_metric_columns <- c("rmse", "rmse_zero", "median_abs_error", "median_abs_error_zero",
                        "extreme_fraction", "bias", "error_sd", "truth_mean", "truth_sd",
                        "skill", "slope", "spearman", "auc_beneficial",
                        "edge_rmse", "edge_sign_accuracy", "coverage95", "interval_score95",
                        "prediction_fraction")

stt_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# One row per landscape x scenario x method x distance.
stt_landscape_scenario_means <- function(metrics) {
  key <- c("landscape_id", "wavelength_class", "scenario", "method", "stratum")
  stats::aggregate(metrics[stt_metric_columns], metrics[key], stt_mean)
}

stt_pooled_scoreboard <- function(cells, metrics) {
  per_landscape <- stats::aggregate(cells[stt_metric_columns],
                                    cells[c("landscape_id", "method", "stratum")], stt_mean)
  board <- stats::aggregate(per_landscape[stt_metric_columns],
                            per_landscape[c("method", "stratum")], stt_mean)
  per_fit <- metrics[!duplicated(metrics[c("dataset_id", "method")]), ]
  runtime <- data.frame(method = sort(unique(per_fit$method)), stringsAsFactors = FALSE)
  runtime$failure_fraction <- as.numeric(tapply(per_fit$status != "ok", per_fit$method, mean)[runtime$method])
  runtime$median_seconds <- as.numeric(tapply(per_fit$seconds, per_fit$method, stats::median)[runtime$method])
  merge(board, runtime, by = "method")
}

# Paired difference of a metric between a candidate and a baseline, by
# distance, with a landscape bootstrap interval.
stt_paired_comparison <- function(cells, candidate, baseline, metric = "rmse",
                                  n_boot = 2000L, seed = 20261110L) {
  per_landscape <- stats::aggregate(cells[metric], cells[c("landscape_id", "method", "stratum")], stt_mean)
  rows <- lapply(c("d0", "d1", "d2"), function(d) {
    pick <- function(method) {
      part <- per_landscape[per_landscape$method == method & per_landscape$stratum == d, ]
      stats::setNames(part[[metric]], part$landscape_id)
    }
    a <- pick(candidate)
    b <- pick(baseline)
    common <- intersect(names(a)[is.finite(a)], names(b)[is.finite(b)])
    if (length(common) < 3L) return(NULL)
    difference <- a[common] - b[common]
    set.seed(seed)
    boot <- replicate(n_boot, mean(sample(difference, replace = TRUE)))
    data.frame(
      candidate = candidate, baseline = baseline, metric = metric, stratum = d,
      n_landscapes = length(common),
      baseline_mean = mean(b[common]), candidate_mean = mean(a[common]),
      difference = mean(difference),
      ci_low = unname(stats::quantile(boot, 0.025)), ci_high = unname(stats::quantile(boot, 0.975)),
      relative_change = mean(difference) / mean(b[common]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Scenario x distance cells: candidate RMSE relative to the baseline.
stt_safety_cells <- function(cells, candidate, baseline, metric = "rmse") {
  cells <- cells[cells$stratum %in% c("d0", "d1", "d2"), ]
  by_scenario <- stats::aggregate(cells[metric], cells[c("scenario", "method", "stratum")], stt_mean)
  a <- by_scenario[by_scenario$method == candidate, c("scenario", "stratum", metric)]
  b <- by_scenario[by_scenario$method == baseline, c("scenario", "stratum", metric)]
  names(a)[3L] <- "candidate_value"
  names(b)[3L] <- "baseline_value"
  out <- merge(a, b, by = c("scenario", "stratum"))
  out$ratio <- out$candidate_value / out$baseline_value
  cbind(candidate = candidate, baseline = baseline, out, stringsAsFactors = FALSE)
}

stt_markdown_table <- function(x, digits = 3L) {
  x[] <- lapply(x, function(col) if (is.numeric(col)) formatC(col, digits = digits, format = "f") else as.character(col))
  header <- paste0("| ", paste(names(x), collapse = " | "), " |")
  rule <- paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|")
  body <- apply(x, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  c(header, rule, body)
}

# The acceptance rule for a candidate, evaluated on the three distances.
#
# Predicting no fitness differences is the bar for error, because the original
# ALFA-K is far above it (a few unconstrained fits dominate its RMSE) and the
# current alfak2 pipeline sits on it. Rank agreement is compared with both
# existing methods, since it is what they do carry and it does not depend on
# the reference.
stt_acceptance <- function(cells, candidate, acceptance) {
  methods <- unique(cells$method)
  criterion <- function(label, baseline, metric, lower_is_better) {
    if (!baseline %in% methods) return(NULL)
    paired <- stt_paired_comparison(cells, candidate, baseline, metric)
    if (is.null(paired)) return(NULL)
    data.frame(
      candidate = candidate, criterion = label, stratum = paired$stratum,
      n_landscapes = paired$n_landscapes, baseline_mean = paired$baseline_mean,
      candidate_mean = paired$candidate_mean, difference = paired$difference,
      ci_low = paired$ci_low, ci_high = paired$ci_high,
      passed = if (lower_is_better) paired$ci_high < 0 else paired$ci_low > 0,
      stringsAsFactors = FALSE
    )
  }
  rows <- list(
    criterion("RMSE below the no-information answer", acceptance$zero, "rmse", TRUE),
    criterion("RMSE below current alfak2", acceptance$current, "rmse", TRUE),
    criterion("Spearman above current alfak2", acceptance$current, "spearman", FALSE)
  )
  for (original in acceptance$original) {
    rows[[length(rows) + 1L]] <- criterion(paste("Spearman above", original), original, "spearman", FALSE)
  }
  if (acceptance$current %in% methods) {
    safety <- stt_safety_cells(cells, candidate, acceptance$current)
    worst <- stats::aggregate(ratio ~ stratum, safety, max)
    rows[[length(rows) + 1L]] <- data.frame(
      candidate = candidate, criterion = "Worst scenario RMSE ratio to current alfak2 at most 1.10",
      stratum = worst$stratum, n_landscapes = NA_integer_, baseline_mean = 1, candidate_mean = worst$ratio,
      difference = worst$ratio - 1, ci_low = NA_real_, ci_high = NA_real_, passed = worst$ratio <= 1.10,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

stt_stage_report <- function(config, split) {
  metrics <- utils::read.csv(file.path(config$output_dir, sprintf("metrics_%s.csv", split)),
                             stringsAsFactors = FALSE)
  cells <- stt_landscape_scenario_means(metrics)
  board <- stt_pooled_scoreboard(cells, metrics)
  # A method scored on fewer datasets is averaged over fewer landscapes, so its
  # pooled row is not comparable; the paired tables use common landscapes only.
  board$n_datasets <- as.integer(tapply(metrics$dataset_id, metrics$method,
                                        function(x) length(unique(x)))[board$method])
  board <- board[order(board$stratum, board$method), ]
  methods <- unique(metrics$method)
  acceptance <- config$acceptance
  baselines <- c(acceptance$zero, acceptance$current, acceptance$original)
  # Methods that are neither baselines nor named candidates (sensitivity
  # variants, for instance) are tabulated but not judged.
  candidates <- intersect(acceptance$candidates, methods)
  # Every method against the no-information answer, and the current alfak2 and
  # every candidate against the original ALFA-K settings.
  pairs <- list()
  for (method in setdiff(methods, acceptance$zero)) {
    if (acceptance$zero %in% methods) pairs[[length(pairs) + 1L]] <- c(method, acceptance$zero)
    if (!method %in% acceptance$original) {
      for (original in intersect(acceptance$original, methods)) pairs[[length(pairs) + 1L]] <- c(method, original)
    }
    if (method %in% candidates && acceptance$current %in% methods) {
      pairs[[length(pairs) + 1L]] <- c(method, acceptance$current)
    }
  }
  # RMSE and median absolute error: lower is better. Spearman and AUC: higher
  # is better, and unlike the errors they do not depend on the reference.
  paired <- do.call(rbind, lapply(pairs, function(pair) {
    do.call(rbind, lapply(c("rmse", "median_abs_error", "spearman", "auc_beneficial"), function(metric) {
      stt_paired_comparison(cells, pair[[1L]], pair[[2L]], metric)
    }))
  }))
  if (is.null(paired)) stop("No paired comparison is available yet.", call. = FALSE)
  verdict <- do.call(rbind, lapply(candidates, function(candidate) stt_acceptance(cells, candidate, acceptance)))
  report_dir <- file.path(config$output_dir, "report")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(cells, file.path(report_dir, sprintf("cells_%s.csv", split)), row.names = FALSE)
  utils::write.csv(board, file.path(report_dir, sprintf("scoreboard_%s.csv", split)), row.names = FALSE)
  utils::write.csv(paired, file.path(report_dir, sprintf("paired_%s.csv", split)), row.names = FALSE)
  if (!is.null(verdict)) {
    utils::write.csv(verdict, file.path(report_dir, sprintf("acceptance_%s.csv", split)), row.names = FALSE)
  }
  show <- c("stratum", "method", "n_datasets", "rmse", "rmse_zero", "truth_sd", "median_abs_error",
            "slope", "spearman", "auc_beneficial", "edge_sign_accuracy", "coverage95",
            "extreme_fraction", "failure_fraction", "median_seconds")
  lines <- c(
    sprintf("# Sparse two-timepoint scoreboard (%s landscapes)", split), "",
    sprintf("%d datasets, %d landscapes, %d scenarios. `rmse_zero` is the error of predicting no fitness differences; `truth_sd` is the error of knowing only the mean fitness of the stratum.",
            length(unique(metrics$dataset_id)), length(unique(metrics$landscape_id)), length(unique(metrics$scenario))),
    "", "## Pooled over scenarios", "", stt_markdown_table(board[show]),
    "", "## Acceptance", "",
    if (is.null(verdict)) "No candidate method has been scored yet." else {
      c(sprintf("Candidates passing every criterion: %s.",
                {
                  ok <- vapply(candidates, function(candidate) all(verdict$passed[verdict$candidate == candidate]), logical(1))
                  if (any(ok)) paste(candidates[ok], collapse = ", ") else "none"
                }), "",
        stt_markdown_table(verdict, 4L))
    },
    "", "## Paired differences (candidate minus baseline, landscape bootstrap 95% interval)", "",
    stt_markdown_table(paired[c("candidate", "baseline", "metric", "stratum", "n_landscapes", "baseline_mean",
                                "candidate_mean", "difference", "ci_low", "ci_high", "relative_change")], 4L)
  )
  path <- file.path(report_dir, sprintf("scoreboard_%s.md", split))
  writeLines(lines, path)
  cat(lines, sep = "\n")
  invisible(path)
}
