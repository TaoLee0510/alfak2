#!/usr/bin/env Rscript
# Step 3 check of the plan: on a small example, the posterior that a fit
# reports (mode, and standard deviations from the expected information) is
# compared with Hamiltonian Monte Carlo draws of the same joint posterior at
# the fitted scales, for the anchored, pairwise and state-space methods.
#
# The example is simulate_series_counts() of the package tests (4 time points,
# 21 states). The sampler is that of ../unobserved_bias/laplace_mean.R: whitened
# by the fit's expected information, identity mass, step size adapted in
# warm-up towards an acceptance of 0.75.
#
#   Rscript validation/sparse_joint/series_methods/hmc_check.R --build=cv9 [--iterations=3000] [--warmup=1000]
#
# Output: results/<build>/hmc_check.rds and a table on stdout.

sm_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(sm_script_dir, "setup.R"))

hc_joint <- function(fit) {
  ns <- asNamespace("alfak2")
  rebuilt <- ns$sparse_joint_rebuild(fit)
  problem <- rebuilt$problem
  map <- problem$map
  map$log_scale <- factor(rep(NA, length(fit$posterior$log_scale)))
  joint <- TMB::MakeADFun(rebuilt$tmb_data, fit$posterior[names(problem$parameters)], map = map, DLL = "alfak2",
                          silent = TRUE)
  solver <- rebuilt$solver
  stopifnot(identical(unname(names(joint$par)), solver$names))
  factor <- environment(solver$draw)$factor
  list(fn = joint$fn, gr = function(theta) as.numeric(joint$gr(theta)), mode = solver$mode, names = solver$names,
       to_theta = function(z) solver$mode + as.numeric(Matrix::solve(factor, Matrix::solve(factor, z, system = "Lt"),
                                                                    system = "Pt")),
       to_z_gradient = function(g) as.numeric(Matrix::solve(factor, Matrix::solve(factor, g, system = "P"),
                                                            system = "L")),
       n = length(solver$mode), solver = solver, problem = problem)
}

hc_hmc <- function(post, iterations, warmup, steps, seed) {
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

# Relative fitness of the domain's states for latent vectors (rows of `theta`).
hc_relative <- function(theta, post, fit) {
  g <- theta[, post$names == "g", drop = FALSE]
  self <- 2 * exp(post$problem$kernel$log_self) - 1
  net <- sweep(exp(g), 2L, self, "*")
  s <- fit$summary
  weight <- ifelse(s$major, rowSums(as.matrix(s[grep("^count_t[0-9]+$", names(s))])), 0)
  net - as.numeric(net %*% weight) / sum(weight)
}

hc_main <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    args[[gsub("-", "_", sub("=.*$", "", body))]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  build <- if (is.null(args$build)) stop("--build is required.", call. = FALSE) else args$build
  iterations <- if (is.null(args$iterations)) 3000L else as.integer(args$iterations)
  warmup <- if (is.null(args$warmup)) 1000L else as.integer(args$warmup)
  .libPaths(c(si_library(build, mt_default_output), .libPaths()))
  loadNamespace("alfak2")
  cat("alfak2 from", dirname(find.package("alfak2")), "\n")
  source(file.path(si_repo, "tests", "testthat", "helper-sparse-joint.R"), local = environment())
  environment(simulate_series_counts) <- asNamespace("alfak2")
  sim <- simulate_series_counts()
  rows <- list()
  out <- list()
  for (method in c("anchored", "pairwise", "state_space")) {
    fit <- alfak2::fit_sparse_joint(sim$counts, times = sim$times, beta = sim$beta, max_cn = 4, history = 60,
                                    series = method)
    post <- hc_joint(fit)
    chain <- hc_hmc(post, iterations, warmup, steps = 16L, seed = 11L)
    g <- post$names == "g"
    hmc_rel <- hc_relative(chain$draws, post, fit)
    mode_rel <- fit$summary$relative_fitness
    model_sd <- fit$summary$relative_sd_model
    hmc_mean <- colMeans(hmc_rel)
    hmc_sd <- apply(hmc_rel, 2L, stats::sd)
    major <- fit$summary$major
    observed <- rowSums(as.matrix(fit$summary[grep("^count_t[0-9]+$", names(fit$summary))])) > 0
    summarise <- function(use, group) {
      data.frame(method = method, group = group, n = sum(use), acceptance = chain$acceptance,
                 # Shift of the posterior mean from the reported mode, in
                 # reported (model) standard deviations.
                 median_abs_shift = stats::median(abs(hmc_mean[use] - mode_rel[use]) / model_sd[use]),
                 max_abs_shift = max(abs(hmc_mean[use] - mode_rel[use]) / model_sd[use]),
                 # HMC standard deviation over the reported one.
                 median_sd_ratio = stats::median(hmc_sd[use] / model_sd[use]),
                 min_sd_ratio = min(hmc_sd[use] / model_sd[use]), max_sd_ratio = max(hmc_sd[use] / model_sd[use]),
                 stringsAsFactors = FALSE)
    }
    rows[[length(rows) + 1L]] <- rbind(summarise(major & model_sd > 0, "major"),
                                       summarise(!major & observed, "observed"),
                                       summarise(!major & !observed, "unobserved"))
    out[[method]] <- list(fit = fit, mode = mode_rel, model_sd = model_sd, hmc_mean = hmc_mean, hmc_sd = hmc_sd,
                          acceptance = chain$acceptance)
  }
  table <- do.call(rbind, rows)
  print(table, digits = 3)
  dir.create(file.path(sm_default_output, build), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(table = table, methods = out), file.path(sm_default_output, build, "hmc_check.rds"))
}

if (sys.nframe() == 0L) hc_main()
