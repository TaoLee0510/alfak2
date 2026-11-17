# Level of the balance prior: study setup (README.md in this directory).
#
# The data come from the development simulator of
# benchmark/sparse_two_timepoint, whose files are sourced unchanged. The study
# keeps its own copies of the inputs it fits, under results/inputs/<set>, so
# that one directory holds everything a fit reads:
#   dev      the development benchmark's scoreboard datasets (seed base 71000)
#   short    the short-interval diagnosis set (seed base 93000)
#   diverse  new: populations as diverse as the Salehi passages, by three
#            different mechanisms (seed base 97000)
# Held-out landscapes (82000, 94000, 96000, 98000) are not touched.

cl_repo <- local({
  here <- if (exists("cl_script_dir")) cl_script_dir else getwd()
  normalizePath(file.path(here, "..", "..", ".."))
})
for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R")) {
  source(file.path(cl_repo, "benchmark", "sparse_two_timepoint", part))
}

cl_output <- file.path(cl_repo, "validation", "sparse_joint", "cloud_level", "results")
cl_sets <- c("dev", "short", "diverse")
cl_input_dir <- function(set) file.path(cl_output, "inputs", set)
cl_index_path <- function(set) file.path(cl_input_dir(set), "dataset_index.csv")
cl_dataset_path <- function(set, dataset_id) file.path(cl_input_dir(set), "datasets", paste0(dataset_id, ".rds"))
cl_nodes_path <- function(set, dataset_id) file.path(cl_input_dir(set), "nodes", paste0(dataset_id, ".rds"))
cl_fit_path <- function(set, method, dataset_id) file.path(cl_output, "fits", set, method, paste0(dataset_id, ".rds"))
cl_tasks_path <- function(name = "tasks.tsv") file.path(cl_output, name)

# The alfak2 build: bf86c85 with cl1.patch, which adds `balance_level` to
# fit_sparse_joint() for two time points. At its defaults it gives the fits
# of bf86c85 (results/same_default.R). Installed outside the source tree.
cl_build <- "cl1"
cl_library <- function(build = cl_build) file.path(cl_output, "lib", paste0("alfak2_", build))

# The missegregation rate every fit assumes, as the paper benchmark does for
# the Salehi lines. The diverse set's "rate" levels simulate a higher one.
cl_assumed_beta <- 5e-5

# Arms. `history = "dt"` fixes the history at the interval.
cl_methods <- list(
  # The package default: balance prior, history chosen among the candidates.
  base = list(),
  # No balance prior.
  free = list(initial_state = "free"),
  # Balance prior with a free level, history chosen among the candidates.
  level = list(balance_level = "free"),
  # Balance prior with a free level, history fixed at the interval: one fit.
  level_dt = list(balance_level = "free", history = "dt"),
  # Corrected-rule study (rule/): the same with the history at the population's
  # age (the first time point's day in the simulation), so that the level is
  # measured against the balance a cloud can reach over the population's
  # whole history.
  level_age = list(balance_level = "free", history = "age")
)

# The diverse set: three ways in which a population comes to hold as many
# rare karyotypes as a Salehi passage (half its cells on karyotypes seen
# once, a largest clone of 13%, half to two thirds of its cells in the
# model's domain). The fit always assumes cl_assumed_beta and no calling
# errors. Levels chosen by results/pilot_levels.R (medians of the share of
# singleton cells, the largest clone and the domain's coverage in brackets):
#   rate8   wavelength 0.8, true missegregation rate 8 times the assumed one
#           (0.45, 0.15, 0.68)
#   rate2   wavelength 1.6, rate 2 times the assumed one (0.42, 0.10, 0.65)
#   smooth  wavelength 2.4 at the assumed rate (0.23, 0.08, 0.82): as diverse
#           as a population gets when the model's assumptions hold
#   call2   wavelength 0.8 at the assumed rate; each chromosome of each
#           sampled cell is called one copy off with probability 2%
#   call1   wavelength 1.6 at the assumed rate; calling errors of 1%
cl_levels <- data.frame(
  level = c("rate8", "rate2", "smooth", "call2", "call1"),
  mechanism = c("rate", "rate", "smooth", "call", "call"),
  wavelength = c(0.8, 1.6, 2.4, 0.8, 1.6),
  beta_multiple = c(8, 2, 1, 1, 1),
  call_error = c(0, 0, 0, 0.02, 0.01),
  stringsAsFactors = FALSE
)

# `set`, `seed_base` and `prefix` give the confirmation set `diverse2` (see
# cl_confirm_sets) the same design on new landscapes.
cl_diverse_config <- function(per_level = 12L, set = "diverse", seed_base = 97000L, prefix = "cl") {
  config <- stt_config(cl_input_dir(set))
  grid <- cl_levels[rep(seq_len(nrow(cl_levels)), each = per_level), ]
  index <- seq_len(nrow(grid))
  config$landscapes <- data.frame(
    landscape_id = sprintf("%s_%s_L%02d", prefix, grid$level, rep(seq_len(per_level), nrow(cl_levels))), split = "cl",
    wavelength_class = grid$level, wavelength = grid$wavelength,
    landscape_seed = seed_base + index, level = grid$level, mechanism = grid$mechanism,
    beta_multiple = grid$beta_multiple, call_error = grid$call_error, stringsAsFactors = FALSE
  )
  # Only the population at the first time point is kept: at 16 times the
  # rate a snapshot holds several hundred thousand karyotypes.
  config$abm$end_day <- 1440
  config$abm$record_days <- 1440
  # One interval of 15 days and 400 cells per time point, as a typical pair
  # of Salehi passages.
  config$scenarios <- data.frame(scenario = "t1440_15d", depth_t0 = 400, depth_t1 = 400, interval = 15,
                                 t0 = 1440, bottleneck = 0, shift = 0, role = "diverse", stringsAsFactors = FALSE)
  config
}

# Calling errors: every chromosome of every sampled cell is moved one copy up
# or down with probability `rate` (kept within the viable copy numbers).
cl_call_errors <- function(counts, rate, seed, min_cn = 1L, max_cn = 8L) {
  set.seed(seed)
  karyotypes <- stt_parse_karyotypes(rownames(counts))
  called <- lapply(seq_len(ncol(counts)), function(t) {
    cells <- karyotypes[rep(seq_len(nrow(karyotypes)), counts[, t]), , drop = FALSE]
    hit <- matrix(stats::runif(length(cells)) < rate, nrow(cells))
    step <- matrix(sample(c(-1, 1), length(cells), replace = TRUE), nrow(cells))
    moved <- cells + hit * step
    moved[moved < min_cn] <- min_cn + 1
    moved[moved > max_cn] <- max_cn - 1
    table(stt_format_karyotypes(moved))
  })
  labels <- unique(unlist(lapply(called, names)))
  out <- vapply(called, function(tab) {
    v <- stats::setNames(numeric(length(labels)), labels)
    v[names(tab)] <- as.numeric(tab)
    v
  }, numeric(length(labels)))
  dimnames(out) <- list(labels, colnames(counts))
  out[order(rowSums(out), decreasing = TRUE), , drop = FALSE]
}

# How diverse the counts are, by the statistics quoted for the Salehi
# passages: share of the cells of a time point on karyotypes seen once in it,
# share of its largest karyotype, and share of its cells in the model's domain.
cl_diversity <- function(counts) {
  counts <- as.matrix(counts)
  domain <- alfak2:::sparse_joint_domain(counts, 0.01, 3, 5000, 1L, 8L)
  coverage <- alfak2:::sparse_joint_domain_coverage(counts, domain, 0.01, 3)
  c(singleton_t0 = sum(counts[counts[, 1L] == 1, 1L]) / sum(counts[, 1L]),
    singleton_t1 = sum(counts[counts[, 2L] == 1, 2L]) / sum(counts[, 2L]),
    top_t0 = max(counts[, 1L]) / sum(counts[, 1L]), top_t1 = max(counts[, 2L]) / sum(counts[, 2L]),
    coverage_t0 = coverage$fraction[[1L]], coverage_t1 = coverage$fraction[[2L]],
    n_major = domain$n_major, n_karyotypes = nrow(counts))
}

# One fit with everything the scores read from it.
cl_fit <- function(dataset, nodes, arguments) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  if (identical(arguments$history, "dt")) arguments$history <- dataset$interval
  if (identical(arguments$history, "age")) arguments$history <- dataset$t0
  beta <- if (is.null(dataset$beta_assumed)) dataset$beta else dataset$beta_assumed
  warned <- character()
  started <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(
    do.call(alfak2::fit_sparse_joint, c(
      list(counts, dt = dataset$interval, beta = beta, birth_rate = 1, min_cn = 1L, max_cn = 8L), arguments)),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  seconds <- proc.time()[["elapsed"]] - started
  prediction <- stats::predict(fit, nodes$karyotype)
  compat <- tryCatch(
    alfak2::as_alfakR(fit, minobs = 10, nboot = 2, n0 = NULL, nb = NULL, correct_efflux = FALSE, xval = FALSE),
    error = function(e) NULL)
  d <- fit$diagnostics
  list(
    status = "ok",
    prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                            fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                            shell = prediction$shell, in_domain = prediction$in_domain, stringsAsFactors = FALSE),
    majors = fit$summary$karyotype[fit$summary$major],
    landscape = if (is.null(compat)) NULL else compat$landscape[c("k", "mean", "fq", "nn")],
    scales = fit$trend$scales,
    shell_coef = fit$trend$coefficients[[length(fit$trend$coefficients)]],
    diagnostics = d[intersect(c("convergence", "objective", "n_states", "n_major", "dropped_karyotypes",
                                "dropped_cells", "fitted_birth_rate", "initial_state", "history", "balance_level",
                                "balance_level_sd"), names(d))],
    history_profile = d$history_profile, coverage = d$domain_coverage, warnings = warned, seconds = seconds
  )
}

# Direction accuracy by distance (validation/sparse_joint/direction_accuracy.R):
# for neighbouring karyotypes one copy of one chromosome apart, whether the
# predicted difference in fitness has the sign of the true one; a tie counts
# one half.
cl_direction <- function(table) {
  cn <- stt_parse_karyotypes(table$karyotype)
  index_of <- stats::setNames(seq_len(nrow(table)), table$karyotype)
  edges <- do.call(rbind, lapply(seq_len(ncol(cn)), function(chr) {
    moved <- cn
    moved[, chr] <- moved[, chr] + 1
    child <- index_of[stt_format_karyotypes(moved)]
    keep <- which(!is.na(child))
    cbind(keep, unname(child[keep]))
  }))
  d <- pmax(table$distance[edges[, 1L]], table$distance[edges[, 2L]])
  dt <- table$truth_relative[edges[, 2L]] - table$truth_relative[edges[, 1L]]
  dp <- table$pred_relative[edges[, 2L]] - table$pred_relative[edges[, 1L]]
  dp[abs(dp) < 1e-12] <- 0
  score <- ifelse(dp == 0, 0.5, as.numeric(sign(dp) == sign(dt)))
  out <- vapply(0:2, function(k) {
    use <- dt != 0 & d == k
    if (any(use)) mean(score[use]) else NA_real_
  }, numeric(1))
  stats::setNames(out, paste0("direction_d", 0:2))
}

cl_wasserstein <- function(ka, ma, kb, mb) {
  a <- is.finite(ma) & ma > 0
  b <- is.finite(mb) & mb > 0
  if (!any(a) || !any(b)) return(NA_real_)
  tryCatch(transport::wasserstein(
    transport::wpp(stt_parse_karyotypes(ka[a]), mass = ma[a] / sum(ma[a])),
    transport::wpp(stt_parse_karyotypes(kb[b]), mass = mb[b] / sum(mb[b]))), error = function(e) NA_real_)
}

# The forecast that selection alone makes over the interval, from the cells
# sampled at the first time point: their frequencies times exp(fitness x
# interval), renormalised, with the predicted and with the true fitness. As
# the paper's Fig. 2d, the size of the predicted change (Wasserstein distance
# from the first time point to the forecast) is set against the observed
# change (from the first to the second time point's cells).
cl_forecast <- function(table, interval) {
  observed <- table$distance == 0L
  first <- observed & table$count_t0 > 0
  second <- observed & table$count_t1 > 0
  k0 <- table$karyotype[first]
  w <- table$count_t0[first] / sum(table$count_t0[first])
  variance <- function(f) sum(w * (f - sum(w * f))^2)
  forecast <- function(f) {
    m <- w * exp((f - max(f)) * interval)
    m / sum(m)
  }
  pred <- forecast(table$pred_relative[first])
  true <- forecast(table$truth_relative[first])
  k1 <- table$karyotype[second]
  m1 <- table$count_t1[second]
  c(var_pred = variance(table$pred_relative[first]), var_true = variance(table$truth_relative[first]),
    change_observed = cl_wasserstein(k0, w, k1, m1),
    change_pred = cl_wasserstein(k0, w, k0, pred), change_true = cl_wasserstein(k0, w, k0, true),
    error_pred = cl_wasserstein(k0, pred, k1, m1), error_true = cl_wasserstein(k0, true, k1, m1))
}

# Landscape of as_alfakR() (frequent karyotypes with at least 10 cells, their
# neighbours and the neighbours' neighbours) against the truth, as the paper
# benchmark scores a landscape.
cl_landscape <- function(landscape, truth) {
  if (is.null(landscape)) return(c(landscape_r = NA_real_, landscape_r_fq = NA_real_, landscape_n_fq = NA_real_))
  f <- stt_true_growth_rate(landscape$k, truth)
  c(landscape_r = stt_safe_cor(landscape$mean, f, "pearson"),
    landscape_r_fq = stt_safe_cor(landscape$mean[landscape$fq], f[landscape$fq], "pearson"),
    landscape_n_fq = sum(landscape$fq))
}

# All scores of one fit, as one row.
cl_score <- function(dataset, nodes, fit) {
  truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  score <- stt_score_predictions(nodes, truth_rate, fit$prediction, 5)
  metrics <- score$metrics
  pick <- function(column) {
    out <- metrics[[column]][match(c("d0", "d1", "d2"), metrics$stratum)]
    stats::setNames(out, paste0(column, "_d", 0:2))
  }
  d <- fit$diagnostics
  value <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)[[1L]]
  profile <- fit$history_profile
  c(pick("rmse"), pick("rmse_zero"), pick("bias"), pick("spearman"), pick("slope"), pick("coverage95"),
    cl_direction(score$nodes), cl_forecast(score$nodes, dataset$interval), cl_landscape(fit$landscape, dataset$truth),
    scale_dosage = fit$scales[["dosage"]], scale_curvature = fit$scales[["curvature"]],
    scale_residual = fit$scales[["residual"]], shell_coef = fit$shell_coef,
    history = value(d$history), balance_level = value(d$balance_level),
    balance_level_sd = value(d$balance_level_sd), objective = value(d$objective),
    convergence = value(d$convergence), n_states = value(d$n_states), n_major = value(d$n_major),
    candidates = if (is.null(profile)) 1 else nrow(profile),
    gap_second = if (is.null(profile) || nrow(profile) < 2L) NA_real_ else diff(sort(profile$objective)[1:2]),
    warnings = length(fit$warnings), seconds = fit$seconds)
}

# Confirmation of the switching rule (README, "Step 1 of the follow-up"):
# the designs of dev, short and diverse on new landscapes. Seed bases 101000,
# 102000 and 103000 are used nowhere else.
cl_confirm_sets <- list(
  dev2 = list(kind = "dev", seed_base = 101000L, prefix = "c2d"),
  short2 = list(kind = "short", seed_base = 102000L, prefix = "c2s"),
  diverse2 = list(kind = "diverse", seed_base = 103000L, prefix = "c2v"),
  # Corrected-rule study (rule/README.md), step 0 (2): populations of the
  # ALFA-K paper's size and age, at intervals of one, three and seven passages.
  paper3 = list(kind = "paper", seed_base = 104000L, prefix = "c3p")
)

# Populations of the ALFA-K paper's size (each passage grows from 2e5 to 2e7
# cells) and age (first time point at day 120), at intervals of 6, 18 and 42
# days (one, three and seven passages), 1000 cells per time point, as in
# lattice_residual/setup.R. Wavelengths 0.2-1.6, five landscapes each, under
# three conditions: the assumed missegregation rate is the true one (`rate1`),
# the true rate is 8 times the assumed one (`rate8`), or each chromosome of a
# sampled cell is called one copy off with probability 2% (`call2`).
cl_paper_config <- function(set, seed_base, prefix) {
  config <- stt_config(cl_input_dir(set))
  conditions <- data.frame(level = c("rate1", "rate8", "call2"), mechanism = c("none", "rate", "call"),
                           beta_multiple = c(1, 8, 1), call_error = c(0, 0, 0.02), stringsAsFactors = FALSE)
  classes <- data.frame(wavelength_class = c("w0.2", "w0.4", "w0.8", "w1.6"), wavelength = c(0.2, 0.4, 0.8, 1.6),
                        stringsAsFactors = FALSE)
  grid <- merge(conditions, classes[rep(seq_len(nrow(classes)), each = 5L), ], by = NULL)
  grid <- grid[order(grid$level, grid$wavelength), ]
  index <- seq_len(nrow(grid))
  config$landscapes <- data.frame(
    landscape_id = sprintf("%s_%s_%s_L%02d", prefix, grid$level, sub("w", "w", grid$wavelength_class),
                           stats::ave(index, grid$level, grid$wavelength, FUN = seq_along)),
    split = set, wavelength_class = grid$wavelength_class, wavelength = grid$wavelength,
    landscape_seed = seed_base + index, level = grid$level, mechanism = grid$mechanism,
    beta_multiple = grid$beta_multiple, call_error = grid$call_error, stringsAsFactors = FALSE)
  config$abm$pop_size <- 2e5
  config$abm$max_pop <- 2e7
  config$abm$end_day <- 120
  config$abm$record_days <- 60
  intervals <- c(6, 18, 42)
  config$scenarios <- data.frame(scenario = sprintf("t0120_%02dd", intervals), depth_t0 = 1000, depth_t1 = 1000,
                                 interval = intervals, t0 = 120, bottleneck = 0, shift = 0, role = "paper",
                                 stringsAsFactors = FALSE)
  config
}

# The development scoreboard's design: wavelengths 0.4-1.6, three landscapes
# each, one patient, its 14 scenarios from day 1440.
cl_dev_config <- function(set, seed_base, prefix) {
  config <- stt_config(cl_input_dir(set))
  classes <- data.frame(wavelength_class = c("w0.4", "w0.8", "w1.2", "w1.6"), wavelength = c(0.4, 0.8, 1.2, 1.6),
                        stringsAsFactors = FALSE)
  grid <- classes[rep(seq_len(nrow(classes)), each = 3L), ]
  index <- seq_len(nrow(grid))
  config$landscapes <- data.frame(landscape_id = sprintf("%s_L%02d", prefix, index), split = set,
                                  wavelength_class = grid$wavelength_class, wavelength = grid$wavelength,
                                  landscape_seed = seed_base + index, stringsAsFactors = FALSE)
  config$patients_per_landscape <- 1L
  config$scoreboard_patients <- 1L
  config
}

# The short-interval diagnosis set's design (short_interval/setup.R): wavelengths
# 0.2-1.6, five landscapes each, first time point at day 120 or 1440,
# intervals of 6, 18 and 42 days, 1000 cells at each time point.
cl_short_config <- function(set, seed_base, prefix) {
  config <- stt_config(cl_input_dir(set))
  classes <- data.frame(wavelength_class = c("w0.2", "w0.4", "w0.8", "w1.6"), wavelength = c(0.2, 0.4, 0.8, 1.6),
                        stringsAsFactors = FALSE)
  grid <- classes[rep(seq_len(nrow(classes)), each = 5L), ]
  index <- seq_len(nrow(grid))
  config$landscapes <- data.frame(landscape_id = sprintf("%s_L%02d", prefix, index), split = set,
                                  wavelength_class = grid$wavelength_class, wavelength = grid$wavelength,
                                  landscape_seed = seed_base + index, stringsAsFactors = FALSE)
  config$abm$end_day <- 1440
  design <- expand.grid(interval = c(6, 18, 42), t0 = c(120, 1440))
  config$scenarios <- data.frame(scenario = sprintf("t%04d_%02dd", design$t0, design$interval),
                                 depth_t0 = 1000, depth_t1 = 1000, interval = design$interval, t0 = design$t0,
                                 bottleneck = 0, shift = 0, role = "short", stringsAsFactors = FALSE)
  config
}

cl_set_config <- function(set) {
  spec <- cl_confirm_sets[[set]]
  if (is.null(spec)) {
    if (set == "diverse") return(cl_diverse_config())
    stop("No configuration for set ", set, call. = FALSE)
  }
  switch(spec$kind,
    dev = cl_dev_config(set, spec$seed_base, spec$prefix),
    short = cl_short_config(set, spec$seed_base, spec$prefix),
    paper = cl_paper_config(set, spec$seed_base, spec$prefix),
    diverse = cl_diverse_config(set = set, seed_base = spec$seed_base, prefix = spec$prefix))
}
