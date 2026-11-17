#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "simulation_backends.R"))
args <- commandArgs(TRUE)
output <- if (length(args)) args[[1L]] else file.path(script_dir, "results", "simulator_validation.csv")
replicates <- 400L; n0 <- 100L; birth <- 1; dt <- 0.5
expected <- n0 * exp(birth * dt)
theoretical_sd <- sqrt(n0 * exp(birth * dt) * (exp(birth * dt) - 1))
exact <- vapply(seq_len(replicates), function(i) {
  sim <- simulate_exact_division_abm(matrix(2L, 1, 1), birth, 0, dt, 1, n0, 900000L + i)
  sum(sim$pop1)
}, numeric(1))
rows <- list(data.frame(backend = "exact_gillespie", time_step = NA_real_, mean_population = mean(exact),
                        mean_se = stats::sd(exact) / sqrt(replicates), theory_exact = expected,
                        theory_tau = expected, standardized_error = (mean(exact) - expected) / (theoretical_sd / sqrt(replicates))))
for (step in c(0.25, 0.05, 0.01)) {
  n_steps <- ceiling(dt / step); actual_step <- dt / n_steps
  population <- vapply(seq_len(replicates), function(i) {
    set.seed(910000L + i)
    n <- n0
    for (j in seq_len(n_steps)) n <- n + min(n, stats::rpois(1, n * birth * actual_step))
    n
  }, numeric(1))
  tau_expected <- n0 * (1 + birth * actual_step)^n_steps
  rows[[length(rows) + 1L]] <- data.frame(backend = "tau_leap", time_step = actual_step,
    mean_population = mean(population), mean_se = stats::sd(population) / sqrt(replicates),
    theory_exact = expected, theory_tau = tau_expected,
    standardized_error = (mean(population) - tau_expected) / (stats::sd(population) / sqrt(replicates)))
}
result <- do.call(rbind, rows)
stopifnot(abs(result$standardized_error[1]) < 5, all(abs(result$standardized_error[-1]) < 5),
          all(diff(abs(result$theory_tau[-1] - expected)) < 0))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE)
cat("PASS: exact beta-zero Yule mean within Monte Carlo error; small-step tau expectations converge.\n")
