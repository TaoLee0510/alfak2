#!/usr/bin/env Rscript
# Does taking the posterior mode make the predicted drop of unobserved
# karyotypes too shallow? (README.md in this directory, step 2.)
#
# fit_sparse_joint() fixes the scales at their marginal-likelihood estimate
# and reports the joint posterior mode of everything else: initial abundances,
# log birth rates, trend coefficients and the depth ratio. Zero counts make the
# posterior of an unobserved karyotype's fitness one-sided, so its mean may lie
# below its mode. Here the same posterior, with the same fixed scales, is
# sampled by Hamiltonian Monte Carlo and its mean compared with the mode, for
# the domain's karyotypes and for the trend's predictions outside the domain,
# and both are scored against the truth.
#
# The sampler works in coordinates whitened by the model's own expected
# information at the mode (the matrix whose inverse gives the reported
# uncertainty): theta = mode + A z with A A' = information^-1, and identity
# mass. Step size is adapted during warm-up towards an acceptance of 0.75.
#
#   Rscript validation/sparse_joint/unobserved_bias/laplace_mean.R <dataset_id> [--source=rugged|diagnosis|dev] \
#     [--iterations=1500] [--warmup=500] [--steps=16]
#
# Output: results/<source>/laplace_mean/<dataset_id>.rds

if (!exists("ub_script_dir")) ub_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(ub_script_dir, "setup.R"))

# The joint negative log posterior of all latent variables at the fit's
# scales, its mode, and the whitening transform.
ub_joint_posterior <- function(fit) {
  ns <- asNamespace("alfak2")
  rebuilt <- ns$sparse_joint_rebuild(fit)
  problem <- rebuilt$problem
  parameters <- fit$posterior[c("zeta", "g", "coef", "theta1", "log_scale")]
  map <- list(log_scale = factor(rep(NA, length(parameters$log_scale))))
  joint <- TMB::MakeADFun(rebuilt$tmb_data, parameters, map = map, DLL = "alfak2", silent = TRUE)
  solver <- rebuilt$solver
  names_latent <- solver$names
  stopifnot(identical(unname(names(joint$par)), names_latent))
  mode <- solver$mode
  factor <- environment(solver$draw)$factor
  # theta = mode + A z, A = P' L^-T; the gradient in z is A' grad = L^-1 P grad.
  to_theta <- function(z) mode + as.numeric(Matrix::solve(factor, Matrix::solve(factor, z, system = "Lt"), system = "Pt"))
  to_z_gradient <- function(g) as.numeric(Matrix::solve(factor, Matrix::solve(factor, g, system = "P"), system = "L"))
  list(fn = joint$fn, gr = function(theta) as.numeric(joint$gr(theta)), mode = mode, names = names_latent,
       to_theta = to_theta, to_z_gradient = to_z_gradient, n = length(mode), problem = problem)
}

ub_hmc <- function(post, iterations, warmup, steps, seed) {
  set.seed(seed)
  U <- function(z) post$fn(post$to_theta(z))
  grad_U <- function(z) post$to_z_gradient(post$gr(post$to_theta(z)))
  z <- numeric(post$n)
  u <- U(z)
  g <- grad_U(z)
  epsilon <- 0.5 / post$n^0.25
  log_eps_bar <- 0
  h_bar <- 0
  mu <- log(10 * epsilon)
  draws <- matrix(NA_real_, iterations, post$n)
  accepted <- numeric(iterations + warmup)
  for (it in seq_len(iterations + warmup)) {
    p <- stats::rnorm(post$n)
    z_new <- z
    g_new <- g
    p_new <- p - 0.5 * epsilon * g_new
    ok <- TRUE
    for (s in seq_len(steps)) {
      z_new <- z_new + epsilon * p_new
      g_new <- grad_U(z_new)
      if (any(!is.finite(g_new))) { ok <- FALSE; break }
      if (s < steps) p_new <- p_new - epsilon * g_new
    }
    if (ok) {
      p_new <- p_new - 0.5 * epsilon * g_new
      u_new <- U(z_new)
      log_ratio <- u - u_new + 0.5 * (sum(p^2) - sum(p_new^2))
    } else log_ratio <- -Inf
    alpha <- if (is.finite(log_ratio)) min(1, exp(log_ratio)) else 0
    if (stats::runif(1) < alpha) {
      z <- z_new
      u <- u_new
      g <- g_new
    }
    accepted[it] <- alpha
    if (it <= warmup) {
      # Dual averaging of the step size (Hoffman and Gelman 2014).
      eta <- 1 / (it + 10)
      h_bar <- (1 - eta) * h_bar + eta * (0.75 - alpha)
      log_eps <- mu - sqrt(it) / 0.05 * h_bar
      epsilon <- exp(log_eps)
      kappa <- it^-0.75
      log_eps_bar <- kappa * log_eps + (1 - kappa) * log_eps_bar
      if (it == warmup) epsilon <- exp(log_eps_bar)
    } else {
      draws[it - warmup, ] <- post$to_theta(z)
    }
  }
  list(draws = draws, acceptance = mean(accepted[-seq_len(warmup)]), epsilon = epsilon)
}

# Relative fitness (net growth rate minus the count-weighted mean of the
# majors) of the domain's karyotypes and of karyotypes outside it predicted by
# the trend, for a vector of latent variables.
ub_relative <- function(theta, post, fit, outside_design, outside_self) {
  g <- theta[post$names == "g"]
  coef <- theta[post$names == "coef"]
  self <- 2 * exp(post$problem$kernel$log_self) - 1
  net <- self * exp(g)
  s <- fit$summary
  weight <- ifelse(s$major, s$count_t0 + s$count_t1, 0)
  reference <- sum(weight * net) / sum(weight)
  list(domain = net - reference,
       outside = if (is.null(outside_design)) numeric(0) else outside_self * exp(as.numeric(outside_design %*% coef)) - reference)
}

# `base_path`: the path of the source's default fit of a dataset, whose
# chosen history is reused.
ub_laplace_mean <- function(config, dataset_id, iterations, warmup, steps, base_path) {
  dataset <- readRDS(stt_dataset_path(config, dataset_id))
  nodes <- stt_dataset_nodes_si(config, dataset)
  base <- readRDS(base_path(dataset_id))
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  fit <- alfak2::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = 1, min_cn = 1L,
                                  max_cn = 8L, history = base$diagnostics$history)
  ns <- asNamespace("alfak2")
  post <- ub_joint_posterior(fit)
  started <- proc.time()[["elapsed"]]
  chain <- ub_hmc(post, iterations, warmup, steps, seed = dataset$seed)
  seconds <- proc.time()[["elapsed"]] - started
  outside <- setdiff(nodes$karyotype, fit$summary$karyotype)
  cn <- ns$parse_karyotypes(outside)
  design <- ns$sparse_joint_design(cn, fit$trend$centre, fit$trend$major, fit$trend$shell_power)$design
  self_out <- 2 * ns$sparse_joint_self_probability(cn, dataset$beta) - 1
  at_mode <- ub_relative(post$mode, post, fit, design, self_out)
  per_draw <- lapply(seq_len(nrow(chain$draws)), function(i) ub_relative(chain$draws[i, ], post, fit, design, self_out))
  domain_draws <- do.call(cbind, lapply(per_draw, `[[`, "domain"))
  outside_draws <- do.call(cbind, lapply(per_draw, `[[`, "outside"))
  quantiles <- function(m) {
    if (!length(m)) return(matrix(numeric(0), 0, 3))
    t(apply(m, 1L, stats::quantile, probs = c(0.025, 0.5, 0.975), names = FALSE))
  }
  q_domain <- quantiles(domain_draws)
  q_outside <- quantiles(outside_draws)
  mean_domain <- rowMeans(domain_draws)
  mean_outside <- if (length(outside_draws)) rowMeans(outside_draws) else numeric(0)
  sd_domain <- apply(domain_draws, 1L, stats::sd)
  # Effective sample size of the mean relative fitness of the unobserved
  # domain karyotypes, from the autocorrelation of that summary.
  s <- fit$summary
  unobserved <- s$count_t0 + s$count_t1 == 0
  trace <- vapply(per_draw, function(d) mean(d$domain[unobserved]), numeric(1))
  ess <- tryCatch({
    rho <- stats::acf(trace, lag.max = 200, plot = FALSE)$acf[-1]
    cut <- which(rho < 0.05)[1]
    length(trace) / (1 + 2 * sum(rho[seq_len(if (is.na(cut)) length(rho) else cut - 1L)]))
  }, error = function(e) NA_real_)
  truth_domain <- stt_true_growth_rate(s$karyotype, dataset$truth)
  truth_outside <- stt_true_growth_rate(outside, dataset$truth)
  weight <- ifelse(s$major, s$count_t0 + s$count_t1, 0)
  reference <- sum(weight * truth_domain) / sum(weight)
  domain <- data.frame(karyotype = s$karyotype, major = s$major, observed = !unobserved, shell = s$shell,
                       mode = at_mode$domain, mean = mean_domain, median = q_domain[, 2L], low = q_domain[, 1L],
                       high = q_domain[, 3L], sd = sd_domain, reported = s$relative_fitness,
                       reported_sd = s$relative_sd, truth = truth_domain - reference, stringsAsFactors = FALSE)
  outside_table <- data.frame(karyotype = outside, distance = nodes$distance[match(outside, nodes$karyotype)],
                              shell = ns$sparse_joint_major_distance(cn, fit$trend$major),
                              mode = at_mode$outside, mean = mean_outside, median = q_outside[, 2L],
                              low = q_outside[, 1L], high = q_outside[, 3L], truth = truth_outside - reference,
                              stringsAsFactors = FALSE)
  list(dataset_id = dataset_id, n_latent = post$n, acceptance = chain$acceptance, epsilon = chain$epsilon,
       ess_unobserved_mean = ess, seconds = seconds, iterations = iterations, warmup = warmup, steps = steps,
       mode_matches_fit = isTRUE(all.equal(at_mode$domain, s$relative_fitness, tolerance = 1e-6)),
       domain = domain, outside = outside_table, trace = trace)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name, default) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else default
  }
  ids <- args[!grepl("^--", args)]
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- ub_config(value("source", "rugged"))
  .libPaths(c(si_library("af362af", config$library_root), .libPaths()))
  loadNamespace("alfak2")
  out_dir <- file.path(config$output_dir, "laplace_mean")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # The source's datasets and default fits, in the source's own results.
  config_source <- switch(config$source, rugged = lr_config(), diagnosis = si_config(), dev = {
    dev <- si_config()
    dev$output_dir <- config$source_dir
    dev
  })
  base_path <- switch(config$source,
    rugged = function(id) si_fit_path(config_source, id, "base"),
    diagnosis = function(id) si_fit_path(config_source, id, "history_profile"),
    dev = function(id) file.path(config$source_dir, "fits", sprintf("%s__history_profile.rds", id)))
  invisible(si_parallel(ids, function(id) {
    path <- file.path(out_dir, paste0(id, ".rds"))
    if (file.exists(path)) return(TRUE)
    result <- ub_laplace_mean(config_source, id, as.integer(value("iterations", "1500")),
                              as.integer(value("warmup", "500")), as.integer(value("steps", "16")), base_path)
    saveRDS(result, path)
    cat(sprintf("%s: %d latent, acceptance %.2f, step %.3f, ESS %.0f, %.0f s, mode matches fit %s\n", id,
                result$n_latent, result$acceptance, result$epsilon, result$ess_unobserved_mean, result$seconds,
                result$mode_matches_fit))
    TRUE
  }, as.integer(value("workers", "6"))))
}
