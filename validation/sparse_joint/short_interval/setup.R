# Short sampling intervals and young populations: diagnosis of two weaknesses
# of fit_sparse_joint() (README.md in this directory).
#
# The data come from the development simulator of
# benchmark/sparse_two_timepoint, whose files are sourced unchanged. This
# study has its own landscapes, seeds and output directory, so nothing of the
# development benchmark (its landscape numbering, caches or results) is
# touched.

si_repo <- local({
  here <- if (exists("si_script_dir")) si_script_dir else getwd()
  normalizePath(file.path(here, "..", "..", ".."))
})
for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R")) {
  source(file.path(si_repo, "benchmark", "sparse_two_timepoint", part))
}

si_default_output <- file.path(si_repo, "validation", "sparse_joint", "short_interval", "results")

# Library holding a given alfak2 build, installed outside the source tree.
si_library <- function(build, output_dir = si_default_output) {
  file.path(output_dir, "lib", paste0("alfak2_", build))
}

# `set = "heldout"`: the same design on landscapes no one has looked at
# (seed base 94000, ids "sh_L.."), in results/heldout_short, for the final
# confirmation only.
si_config <- function(output_dir = NULL, set = c("diagnosis", "heldout")) {
  set <- match.arg(set)
  if (is.null(output_dir)) {
    output_dir <- if (set == "diagnosis") si_default_output else file.path(si_default_output, "heldout_short")
  }
  config <- stt_config(output_dir)
  # Wavelength 0.2 is more rugged than any development class; it is the
  # regime of the second weakness. Five landscapes per class, as landscapes
  # are the independent unit.
  classes <- data.frame(wavelength_class = c("w0.2", "w0.4", "w0.8", "w1.6"),
                        wavelength = c(0.2, 0.4, 0.8, 1.6), stringsAsFactors = FALSE)
  per_class <- 5L
  grid <- classes[rep(seq_len(nrow(classes)), each = per_class), ]
  index <- seq_len(nrow(grid))
  # Seed base 93000: the development benchmark uses 71000 and 82000.
  config$landscapes <- data.frame(
    landscape_id = sprintf(if (set == "diagnosis") "si_L%02d" else "sh_L%02d", index), split = "si",
    wavelength_class = grid$wavelength_class, wavelength = grid$wavelength,
    landscape_seed = (if (set == "diagnosis") 93000L else 94000L) + index, stringsAsFactors = FALSE
  )
  config$set <- set
  # Libraries are shared by both sets.
  config$library_root <- si_default_output
  config$patients_per_landscape <- 1L
  config$scoreboard_patients <- 1L
  # Day 120 is about the age of the paper's ABM populations at their windows;
  # day 1440 is the first sample of the development benchmark. The two
  # separate a short interval from a young population. Intervals of 6, 18 and
  # 42 days match one, three and seven passages of the paper's ABM.
  config$abm$end_day <- 1440
  t0_values <- c(120, 1440)
  intervals <- c(6, 18, 42)
  design <- expand.grid(interval = intervals, t0 = t0_values)
  config$scenarios <- data.frame(
    scenario = sprintf("t%04d_%02dd", design$t0, design$interval),
    depth_t0 = 1000, depth_t1 = 1000, interval = design$interval, t0 = design$t0,
    bottleneck = 0, shift = 0, role = "diagnosis", stringsAsFactors = FALSE
  )
  # Cross-validation as in as_alfakR(): frequent karyotypes have at least
  # `minobs` cells over both samples.
  config$cv <- list(minobs = 10)
  # Arms. `build` names the alfak2 library (si_library()); `history = "age"`
  # sets the history to the true age of the population at the first sample,
  # which only a simulation knows; `history = "grid"` gives the candidates of
  # si_history_grid(). Build "si1" is 8327282 with the opt-in trend and
  # history options of this study; its defaults reproduce 8327282 exactly.
# Build "si2" adds `history = "estimate"`; with history fixed it agrees with
# 8327282 to about 1e-7 (the history became a mapped TMB parameter).
  sparse <- function(build = "8327282", ...) {
    list(type = "sparse_joint", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L,
         arguments = list(...))
  }
  config$methods <- list(
    zero = list(type = "zero"),
    # Current default: initial_state = "balance", history = dt.
    balance = sparse(),
    # Probe of change A: history = true age of the population.
    balance_history_age = sparse(history = "age"),
    # Reference arm without the balance prior.
    free = sparse(initial_state = "free"),
    # Change A, usable version: history chosen by marginal likelihood.
    history_grid = sparse("si1", history = "grid"),
    # Change B probes. B1a: a lower log-normal centre for the trend scales
    # (probe only), and an exponential penalty pulling them to zero. B1b:
    # nested trend terms chosen by marginal likelihood.
    trend_centre_low = sparse("si1", trend_scale_centre = 0.001),
    trend_exponential = sparse("si1", trend_scale_prior = "exponential"),
    trend_select = sparse("si1", trend_terms = "select"),
    # Changes A and B together, for each version of B.
    grid_exponential = sparse("si1", history = "grid", trend_scale_prior = "exponential"),
    grid_select = sparse("si1", history = "grid", trend_terms = "select"),
    # Change A, final version: history chosen among the interval and longer.
    history_profile = sparse("si1", history = "profile"),
    # The same choice made by estimating the history as an outer parameter,
    # bounded below by the interval (build si2), in one fit.
    history_estimate = sparse("si2", history = "estimate"),
    # The package with the new default (build si4: history chosen as in
    # history_profile, which it reproduces exactly), at its defaults.
    default_new = sparse("si4")
  )
  config
}

# Candidate histories, in days (generations at the assumed birth rate of one
# per day): the interval itself and a range of population ages that does not
# depend on it.
si_history_grid <- function(interval) unique(c(interval, 30, 100, 300, 1000, 3000))
# The same, keeping only histories at least as long as the interval: the
# marginal likelihood prefers 30-100 days even where the interval is
# 180-1080 days and the interval itself predicts better (history_sweep.R).
si_history_profile <- function(interval) unique(c(interval, c(30, 100, 300, 1000, 3000)[c(30, 100, 300, 1000, 3000) > interval]))

si_dataset_plan <- function(config) {
  landscapes <- config$landscapes
  plan <- list()
  for (i in seq_len(nrow(landscapes))) {
    for (s in seq_len(nrow(config$scenarios))) {
      plan[[length(plan) + 1L]] <- list(
        landscape = landscapes[i, ], patient = 1L, scenario = config$scenarios[s, ],
        seed = as.integer(landscapes$landscape_seed[i] %% 100000L * 1000L + 20L + s)
      )
    }
  }
  plan
}

si_fit_path <- function(config, dataset_id, method_name) {
  file.path(config$output_dir, "fits", method_name, paste0(dataset_id, ".rds"))
}

si_index_path <- function(config) file.path(config$output_dir, "dataset_index.csv")

# Count-weighted spread of true fitness times the interval: how much selection
# the two samples can see (handoff, section 4.1).
si_selection_signal <- function(dataset) {
  counts <- dataset$counts
  rate <- stt_true_growth_rate(rownames(counts), dataset$truth)
  weight <- rowSums(counts) / sum(counts)
  centre <- sum(weight * rate)
  sqrt(sum(weight * (rate - centre)^2)) * dataset$interval
}

# Statistics of a held-out prediction `x` against a target, both centred:
# correlation, spread ratio, R2R, and the slope of the target on `x` (below 1:
# `x` is more dispersed than the target supports).
si_centred_stats <- function(target, x) {
  tc <- target - mean(target)
  xc <- x - mean(x)
  c(r = suppressWarnings(stats::cor(target, x)), s = sqrt(sum(xc^2) / sum(tc^2)),
    R2R = 1 - sum((xc - tc)^2) / sum(tc^2), slope = sum(xc * tc) / sum(xc^2))
}

# Cross-validation of as_alfakR() against the truth, as
# alfak2-paper-benchmark/scripts/cv_score_analysis.R computes it.
si_cv_metrics <- function(cv, fq, truth) {
  o <- unname(cv$observations)
  p <- unname(cv$predictions)
  k <- cv$xval_data$k
  ok <- is.finite(o) & is.finite(p)
  out <- c(cv_pairs = sum(ok), R2R = cv$R2R)
  if (sum(ok) < 3L) return(out)
  o <- o[ok]
  p <- p[ok]
  k <- k[ok]
  f <- stt_true_growth_rate(k, truth)
  po <- si_centred_stats(o, p)
  pt <- si_centred_stats(f, p)
  ot <- si_centred_stats(f, o)
  fold <- unlist(lapply(seq_along(fq), function(i) {
    members <- c(fq[i], alfak2:::alfakR_neighbours(fq[i]))
    stats::setNames(rep(i, length(members)), members)
  }))
  fold <- fold[!duplicated(names(fold))]
  fid <- fold[k]
  e <- (p - mean(p)) - (o - mean(o))
  keep <- !is.na(fid)
  between <- if (any(keep)) sum(tapply(e[keep], fid[keep], function(v) length(v) * mean(v)^2)) / sum(e[keep]^2) else NA
  c(out, r_po = po[["r"]], s_po = po[["s"]],
    heldout_r = pt[["r"]], heldout_slope = pt[["slope"]], heldout_R2R = pt[["R2R"]],
    full_r = ot[["r"]], full_slope = ot[["slope"]],
    bias_fq = mean(e[k %in% fq]) / stats::sd(o), n_folds = length(unique(fid[keep])), between = between)
}

# Landscape of as_alfakR() (frequent karyotypes, their neighbours and the
# neighbours' neighbours) against the truth: r and centred R-squared, as in
# the paper's Supplementary Note 3.
si_landscape_metrics <- function(landscape, truth) {
  f <- stt_true_growth_rate(landscape$k, truth)
  m <- landscape$mean
  r2r <- function(obs, pred) {
    obs <- obs - mean(obs)
    pred <- pred - mean(pred)
    1 - sum((pred - obs)^2) / sum(obs^2)
  }
  d2 <- !(landscape$fq | landscape$nn)
  c(landscape_r = stats::cor(m, f), landscape_R2 = r2r(f, m),
    landscape_r_fq = if (sum(landscape$fq) >= 3L) stats::cor(m[landscape$fq], f[landscape$fq]) else NA,
    landscape_r_nn = stats::cor(m[landscape$nn], f[landscape$nn]),
    landscape_r_d2 = stats::cor(m[d2], f[d2]), n_fq = sum(landscape$fq))
}

# One fit with everything the diagnosis reads from it.
si_fit_sparse_joint <- function(config, dataset, nodes, method) {
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  arguments <- method$arguments
  if (identical(arguments$history, "age")) arguments$history <- dataset$t0
  if (identical(arguments$history, "grid")) arguments$history <- si_history_grid(dataset$interval)
  if (identical(arguments$history, "profile")) arguments$history <- si_history_profile(dataset$interval)
  started <- proc.time()[["elapsed"]]
  fit <- do.call(alfak2::fit_sparse_joint, c(
    list(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = method$birth_rate,
         min_cn = method$min_cn, max_cn = method$max_cn),
    arguments
  ))
  fit_seconds <- proc.time()[["elapsed"]] - started
  prediction <- stats::predict(fit, nodes$karyotype)
  started <- proc.time()[["elapsed"]]
  compat <- tryCatch(
    alfak2::as_alfakR(fit, minobs = config$cv$minobs, nboot = 2, n0 = NULL, nb = NULL, correct_efflux = FALSE),
    error = function(e) conditionMessage(e)
  )
  cv_seconds <- proc.time()[["elapsed"]] - started
  fq <- rownames(fit$data$counts)[rowSums(fit$data$counts) >= config$cv$minobs]
  cv <- if (is.character(compat)) c(cv_error = 1) else
    c(si_cv_metrics(compat$xval, fq, dataset$truth), si_landscape_metrics(compat$landscape, dataset$truth))
  d <- fit$diagnostics
  list(
    status = "ok",
    prediction = data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                            fitness_sd = prediction$relative_sd, model_sd = prediction$relative_sd_model,
                            shell = prediction$shell, stringsAsFactors = FALSE),
    majors = fit$summary$karyotype[fit$summary$major],
    domain = fit$summary[c("karyotype", "count_t0", "count_t1", "major", "shell", "relative_fitness",
                           "relative_sd")],
    scales = fit$trend$scales,
    diagnostics = d[c("convergence", "message", "objective", "finite_uncertainty", "n_states", "n_major",
                      "dropped_karyotypes", "fitted_birth_rate", "initial_state", "history")],
    trend_terms = if (is.null(fit$trend$terms)) "full" else fit$trend$terms,
    history_profile = d$history_profile, trend_profile = d$trend_profile,
    cv = cv,
    cv_pairs = if (is.character(compat)) NULL else compat$xval$xval_data,
    cv_error = if (is.character(compat)) compat else NA_character_,
    fit_seconds = fit_seconds, cv_seconds = cv_seconds
  )
}

si_run_method <- function(config, dataset, nodes, method) {
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(
    switch(method$type,
      zero = list(status = "ok", prediction = data.frame(karyotype = nodes$karyotype, fitness = 0,
                                                         fitness_sd = NA_real_, stringsAsFactors = FALSE)),
      sparse_joint = si_fit_sparse_joint(config, dataset, nodes, method),
      stop("Unknown method type: ", method$type, call. = FALSE)
    ),
    error = function(e) list(status = "failed", prediction = NULL, detail = conditionMessage(e))
  )
  result$seconds <- proc.time()[["elapsed"]] - started
  result
}

# Scores of one fit: the development benchmark's per-distance metrics, the
# bias of observed karyotypes split into majors and the rest (hypothesis H1,
# prediction i), and the cross-validation and landscape metrics.
si_score_fit <- function(config, dataset, nodes, fit) {
  truth_rate <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  prediction <- if (is.null(fit$prediction)) {
    data.frame(karyotype = character(), fitness = numeric(), stringsAsFactors = FALSE)
  } else fit$prediction
  score <- stt_score_predictions(nodes, truth_rate, prediction, config$scoring$reference_min_count)
  table <- score$nodes
  observed <- table$distance == 0L
  is_major <- table$karyotype %in% fit$majors
  error <- table$pred_relative - table$truth_relative
  extra <- c(
    bias_major = if (any(observed & is_major)) mean(error[observed & is_major]) else NA,
    bias_observed_nonmajor = if (any(observed & !is_major)) mean(error[observed & !is_major]) else NA,
    n_observed_nonmajor = sum(observed & !is_major)
  )
  list(metrics = score$metrics, extra = extra)
}
