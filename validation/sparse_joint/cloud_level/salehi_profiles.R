#!/usr/bin/env Rscript
# Cloud-level study, step 0: what the stored Salehi fits of the paper
# benchmark (round 05a8855) say about the history of the balance prior. Read
# only: no fit is repeated; a two-time-point fit's posterior is rebuilt at its
# stored mode to read the residuals of the balance prior.
#
# For every alfak2 fit (designs "ends" and "all"; methods alfak2,
# alfak2_pairwise, alfak2_state_space) one row with the chosen history, the
# fitted scales, the domain's coverage, and:
#   cloud_cells_t0      cells at the first time point on domain karyotypes
#                       that are not major
#   major_cells_t0      cells at the first time point on majors
#   neutral_per_day     what the majors would feed their neighbours per day of
#                       history if the neighbours were as fit as the majors:
#                       sum over edges major -> non-major of 2 P_ij b n_i
#   implied_history     cloud_cells_t0 / neutral_per_day: the history under
#                       which a neutral cloud has the observed size
#   balance_residual_*  mean and sd over the fed states of log abundance at
#                       the first time point minus the balance prior's centre
#                       (two time points only)
# and the profile of the marginal likelihood over the candidate histories
# (one row per candidate, in a second table).
#
#   Rscript salehi_profiles.R OUT_DIR

suppressPackageStartupMessages(library(alfak2))
out_dir <- commandArgs(TRUE)[1]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
root <- "/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_05a8855/fits/salehi"
ids <- list.files(root)
fits <- list()
profiles <- list()
for (id in ids) for (design in c("ends", "all")) for (method in c("alfak2", "alfak2_pairwise", "alfak2_state_space")) {
  path <- file.path(root, id, design, method, "minobs_5", "alfak2_fit.Rds")
  if (!file.exists(path)) next
  fit <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(fit) || !inherits(fit, "alfak2_sparse_joint_fit")) next
  d <- fit$diagnostics
  s <- fit$summary
  counts <- fit$data$counts
  n_time <- ncol(counts)
  key <- data.frame(id = id, design = design, method = method, stringsAsFactors = FALSE)
  scales <- fit$trend$scales
  row <- cbind(key, data.frame(
    n_time = n_time, dt = fit$data$dt, cells_t0 = sum(counts[, 1L]), cells_t1 = sum(counts[, n_time]),
    history = d$history[[1L]], n_candidates = if (is.null(d$history_profile)) 1L else nrow(d$history_profile),
    objective = d$objective, convergence = d$convergence, n_states = d$n_states, n_major = d$n_major,
    dropped_karyotypes = d$dropped_karyotypes, dropped_cells = d$dropped_cells,
    coverage_t0 = d$domain_coverage$fraction[[1L]],
    coverage_t1 = d$domain_coverage$fraction[[nrow(d$domain_coverage)]],
    fitted_birth_rate = d$fitted_birth_rate,
    scale_dosage = scales[["dosage"]], scale_curvature = scales[["curvature"]],
    scale_residual = scales[["residual"]],
    shell_coef = fit$trend$coefficients[[length(fit$trend$coefficients)]],
    cloud_cells_t0 = NA_real_, major_cells_t0 = NA_real_, neutral_per_day = NA_real_,
    cloud_states = NA_real_, cloud_states_seen_t0 = NA_real_,
    var_major = NA_real_, sd_relative_minor_observed = NA_real_,
    balance_residual_mean = NA_real_, balance_residual_sd = NA_real_, balance_fed = NA_real_
  ))
  # The plain summary (two time points, anchored): states x time points.
  if (!is.null(s$major) && !is.null(s$count_t0) && nrow(s) == d$n_states) {
    major <- s$major
    y0 <- s$count_t0
    row$cloud_cells_t0 <- sum(y0[!major])
    row$major_cells_t0 <- sum(y0[major])
    row$cloud_states <- sum(!major)
    row$cloud_states_seen_t0 <- sum(!major & y0 > 0)
    kernel <- alfak2:::sparse_division_kernel(alfak2:::parse_karyotypes(s$karyotype), fit$data$beta)
    feeds <- major[kernel$from] & !major[kernel$to]
    row$neutral_per_day <- sum(2 * exp(kernel$log_prob[feeds]) * d$birth_rate * y0[kernel$from[feeds]])
    w <- y0[major] / sum(y0[major])
    if (sum(y0[major]) > 0) row$var_major <- sum(w * (s$fitness_mean[major] - sum(w * s$fitness_mean[major]))^2)
    seen_minor <- !major & (y0 + s[[paste0("count_t", n_time - 1L)]]) > 0
    if (sum(seen_minor) > 1L) row$sd_relative_minor_observed <- stats::sd(s$relative_fitness[seen_minor])
  }
  if (n_time == 2L && identical(d$initial_state, "balance")) {
    residual <- tryCatch({
      rebuilt <- alfak2:::sparse_joint_rebuild(fit, posterior = FALSE)
      problem <- alfak2:::sparse_joint_problem(fit$data, fit$settings)
      parameters <- fit$posterior[names(problem$parameters)]
      map <- problem$map
      map$log_scale <- factor(rep(NA, length(parameters$log_scale)))
      objective <- TMB::MakeADFun(problem$tmb_data, parameters, map = map, random = problem$random,
                                  DLL = "alfak2", silent = TRUE, inner.control = fit$settings$inner_control)
      objective$fn(numeric(0))
      report <- objective$report(objective$env$last.par)
      fed <- seq_along(report$balance_residual) > problem$domain$n_major & report$balance_residual != 0
      report$balance_residual[fed]
    }, error = function(e) NULL)
    if (length(residual)) {
      row$balance_residual_mean <- mean(residual)
      row$balance_residual_sd <- stats::sd(residual)
      row$balance_fed <- length(residual)
    }
  }
  fits[[length(fits) + 1L]] <- row
  if (!is.null(d$history_profile)) {
    profiles[[length(profiles) + 1L]] <- cbind(key, d$history_profile[c("history", "objective", "convergence")])
  }
  cat(id, design, method, "\n")
}
fits <- do.call(rbind, fits)
fits$implied_history <- fits$cloud_cells_t0 / fits$neutral_per_day
utils::write.csv(fits, file.path(out_dir, "salehi_fits.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, profiles), file.path(out_dir, "salehi_history_profiles.csv"), row.names = FALSE)
cat("Fits:", nrow(fits), " profile rows:", sum(vapply(profiles, nrow, integer(1))), "\n")
