#!/usr/bin/env Rscript
# Independent first-moment checks for the fixed-initial-population simulator.
script <- grep('^--file=', commandArgs(FALSE), value = TRUE)
directory <- if (length(script)) dirname(normalizePath(sub('^--file=', '', script[[1L]]))) else 'validation/accuracy_round2'
source(file.path(directory, 'process_prepare_data.R'))
output <- commandArgs(TRUE)
output <- if (length(output)) output[[1L]] else '/private/tmp/alfak2_iteration2_20261104_0518/process-data/generator_validation.csv'
rows <- list()
repeats <- 400L; initial <- 50L; time <- .5; rate <- 1
draws <- vapply(seq_len(repeats), function(i)
  process_exact_gillespie(matrix(2L, 1, 1), rate, 0, time, initial, 771000L + i)$population,
  numeric(1))
p <- exp(-rate * time); q <- 1 - p
mean_expected <- initial / p
variance_expected <- initial * q / p^2
kappa4 <- initial * q * (1 + 4 * q + q^2) / p^4
mean_se <- sqrt(variance_expected / repeats)
variance_se <- sqrt((kappa4 + 3 * variance_expected^2 -
  (repeats - 3) / (repeats - 1) * variance_expected^2) / repeats)
rows[[1L]] <- data.frame(case = 'zero_mutation_Yule', metric = c('mean', 'variance'),
  observed = c(mean(draws), var(draws)), expected = c(mean_expected, variance_expected),
  mcse = c(mean_se, variance_se))
k <- matrix(1:3, ncol = 1L); initial_counts <- c(5, 40, 5)
fitness <- c(.8, 1.2, 1); beta <- .03; time <- .3
operator <- process_independent_operator(k, fitness, beta)
expected <- process_deterministic_mean(operator$generator, initial_counts, time)
draws <- vapply(seq_len(repeats), function(i)
  process_exact_gillespie(k, fitness, beta, time, initial_counts, 772000L + i)$population,
  numeric(3))
rows[[2L]] <- data.frame(case = 'bounded_nonzero_mutation', metric = paste0('endpoint_mean_cn', 1:3),
  observed = rowMeans(draws), expected = expected, mcse = apply(draws, 1L, sd) / sqrt(repeats))
result <- do.call(rbind, rows)
result$standardized_difference <- (result$observed - result$expected) / result$mcse
result$passed <- abs(result$standardized_difference) < 5
write.csv(result, output, row.names = FALSE)
print(result, row.names = FALSE)
stopifnot(all(result$passed))

# Check that the fitted v1 numerical mean dynamics agree with the independently
# generated deterministic arm, rather than attributing integration bias to noise.
data_directory <- dirname(output)
frozen <- file.path(dirname(data_directory), 'v1-source')
if (file.exists(file.path(data_directory, 'dataset_manifest.csv')) &&
    file.exists(file.path(frozen, 'R/division-process.R'))) {
  inference <- new.env(parent = globalenv())
  for (path in c('data-validate.R', 'simulate-l1-gp-landscape.R', 'division-process.R'))
    sys.source(file.path(frozen, 'R', path), inference)
  manifest <- read.csv(file.path(data_directory, 'dataset_manifest.csv'), stringsAsFactors = FALSE)
  manifest <- manifest[manifest$mode == 'deterministic_mean' & manifest$sampling_replicate == 1L, ]
  solver_checks <- lapply(manifest$dataset_file, function(path) {
    dataset <- readRDS(path)
    P <- process_independent_operator(dataset$target_karyotypes, dataset$truth, dataset$config$beta)$probability
    nonzero <- which(P > 0, arr.ind = TRUE)
    operator <- list(from0 = nonzero[, 1L] - 1L, to0 = nonzero[, 2L] - 1L,
                      weight = as.numeric(P[nonzero]))
    predicted <- inference$division_forward_r(dataset$truth_pi0, dataset$truth, operator,
      dataset$config$dt, steps = ceiling(dataset$config$dt / .1))
    growth_truth <- log(sum(dataset$latent_population_t1) / sum(dataset$latent_population_t0))
    data.frame(landscape_id = dataset$landscape_id,
      maximum_frequency_error = max(abs(predicted$pi1 - dataset$truth_pi1)),
      log_growth_error = abs(predicted$log_total_growth - growth_truth))
  })
  solver_checks <- do.call(rbind, solver_checks)
  solver_checks$passed <- pmax(solver_checks$maximum_frequency_error, solver_checks$log_growth_error) < 1e-5
  write.csv(solver_checks, file.path(data_directory, 'deterministic_solver_validation.csv'), row.names = FALSE)
  print(solver_checks, row.names = FALSE)
  stopifnot(all(solver_checks$passed))
}
