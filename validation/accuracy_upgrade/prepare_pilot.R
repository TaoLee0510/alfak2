#!/usr/bin/env Rscript
# Independent finite-domain division ABM; uses no inference code or true-fitness calibration.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "config.R"))
source(file.path(script_dir, "simulation_backends.R"))
args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}
profile <- get_arg("profile", "pilot")
if (!profile %in% c("pilot", "farfield")) stop("--profile must be pilot or farfield")
out <- get_arg("output", if (length(args) && !startsWith(args[[1L]], "--")) args[[1L]] else file.path(script_dir, "data"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
config <- accuracy_upgrade_config()
backend <- get_arg("backend", "tau_leap")
if (!backend %in% c("tau_leap", "gillespie")) stop("--backend must be tau_leap or gillespie")
config$simulator_backend <- backend
config$protocol_version <- if (backend == "gillespie") 2L else 1L
config$generating_protocol <- if (backend == "gillespie") "v2_exact_continuous_division" else "v1_tau_leap_division"
config$run_split <- get_arg("split", "train")
if (!config$run_split %in% names(config$split_replicates)) stop("Unknown split")
if (profile == "farfield") {
  config$n_chr <- 5L; config$max_nodes <- 243L; config$beta <- 0
  config$families <- "interaction"
  config$initial_neighbour_mass <- 0
  config$abm_replicates <- 1L; config$sampling_replicates <- 1L
  config$split_replicates <- list(train = 1L, calibration = 3L, test = 5L)
  config$probe_role <- "support-limit stress: zero mutation, one anchor, guaranteed d0-d5; no identification claim"
}
config$profile <- profile
manifest <- accuracy_seed_manifest(config)
utils::write.csv(manifest, file.path(out, "landscape_seed_manifest.csv"), row.names = FALSE)
saveRDS(config, file.path(out, "config.rds"))
k <- as.matrix(expand.grid(rep(list(config$min_cn:config$max_cn), config$n_chr)))
storage.mode(k) <- "integer"
labels <- apply(k, 1L, paste, collapse = ".")
rownames(k) <- labels
founder <- which(rowSums(abs(k - 2L)) == 0L)
distance_matrix <- as.matrix(stats::dist(k, method = "manhattan"))
edges <- which(distance_matrix == 1 & row(distance_matrix) < col(distance_matrix), arr.ind = TRUE)

make_truth <- function(family, seed) {
  set.seed(seed)
  x <- k - 2L
  if (family == "smooth") {
    f <- 1 + as.numeric(x %*% stats::runif(config$n_chr, -0.16, 0.16)) -
      0.06 * rowSums(x^2)
  } else if (family == "rough") {
    centroids <- matrix(stats::runif(16L * config$n_chr, -5, 8), ncol = config$n_chr)
    raw <- numeric(nrow(k))
    for (j in seq_len(nrow(centroids))) {
      d <- sqrt(rowSums(sweep(k, 2L, centroids[j, ], "-")^2))
      raw <- raw + sin(d / 0.22)
    }
    f <- 1 + 0.18 * (raw - raw[founder]) / stats::sd(raw)
  } else {
    co <- stats::runif(3, 0.12, 0.25) * sample(c(-1, 1), 3, replace = TRUE)
    f <- 1 + 0.06 * x[, 1] + co[1] * x[, 1] * x[, 2] +
      co[2] * x[, 3] * x[, 4] + co[3] * x[, 1] * x[, 2] * x[, 3]
  }
  stopifnot(all(f > 0), all(is.finite(f)))
  as.numeric(f)
}

simulate_population <- function(fitness, seed) {
  set.seed(seed)
  p0 <- rep(0, nrow(k))
  neighbours <- which(distance_matrix[founder, ] == 1L)
  p0[founder] <- 1 - config$initial_neighbour_mass
  p0[neighbours] <- config$initial_neighbour_mass / length(neighbours)
  if (config$simulator_backend == "gillespie") {
    return(simulate_exact_division_abm(k, fitness, config$beta, config$dt, p0,
                                      config$initial_population, seed))
  }
  pop <- as.numeric(stats::rmultinom(1, config$initial_population, p0))
  pop0 <- pop
  n_steps <- ceiling(config$dt / config$time_step)
  step <- config$dt / n_steps
  for (s in seq_len(n_steps)) {
    next_pop <- pop
    for (i in which(pop > 0)) {
      divisions <- min(pop[i], stats::rpois(1, pop[i] * fitness[i] * step))
      if (!divisions) next
      next_pop[i] <- next_pop[i] - divisions
      x <- k[i, ]
      faithful <- stats::rbinom(1, divisions, (1 - config$beta)^sum(x))
      next_pop[i] <- next_pop[i] + 2 * faithful
      for (j in seq_len(divisions - faithful)) {
        n_error <- 0L
        while (n_error == 0L) n_error <- stats::rbinom(1, sum(x), config$beta)
        copies <- rep(seq_along(x), x)
        picked <- copies[sample.int(length(copies), n_error)]
        signs <- sample(c(-1L, 1L), n_error, replace = TRUE)
        delta <- tabulate(picked[signs > 0], length(x)) - tabulate(picked[signs < 0], length(x))
        for (daughter in list(x + delta, x - delta)) {
          label <- paste(daughter, collapse = ".")
          ii <- match(label, labels)
          if (!is.na(ii)) next_pop[ii] <- next_pop[ii] + 1
        }
      }
    }
    pop <- next_pop
    if (!sum(pop)) stop("ABM extinction")
  }
  list(pop0 = pop0, pop1 = pop, n_steps = n_steps)
}

run_manifest <- list()
landscapes <- manifest[manifest$split == config$run_split, , drop = FALSE]
for (li in seq_len(nrow(landscapes))) {
  landscape <- landscapes[li, , drop = FALSE]
  truth <- make_truth(landscape$family, landscape$landscape_seed)
  for (abm_rep in config$abm_replicates) {
    abm_seed <- landscape$landscape_seed * 10L + abm_rep
    latent <- simulate_population(truth, abm_seed)
    for (sample_rep in config$sampling_replicates) {
      sampling_seed <- abm_seed * 10L + sample_rep
      set.seed(sampling_seed)
      y0 <- as.numeric(stats::rmultinom(1, config$n0, latent$pop0))
      y1 <- as.numeric(stats::rmultinom(1, config$n1, latent$pop1))
      counts <- cbind(t0 = y0, t1 = y1)
      rownames(counts) <- labels
      observed <- which(rowSums(counts) > 0)
      counts <- counts[observed, , drop = FALSE]
      distance <- apply(distance_matrix[, observed, drop = FALSE], 1L, min)
      if (profile == "farfield") stopifnot(max(distance) == 5L)
      gauge_weight <- numeric(nrow(k))
      gauge_weight[observed] <- counts[, "t0"] / sum(counts[, "t0"])
      dataset_id <- sprintf("%s%s_A%d_S%d", if (profile == "farfield") "farfield_" else "", landscape$landscape_id, abm_rep, sample_rep)
      dataset <- list(dataset_id = dataset_id, family = landscape$family,
                      landscape_id = landscape$landscape_id, split = landscape$split,
                      abm_seed = abm_seed, sampling_seed = sampling_seed, counts = counts,
                      target_karyotypes = k, target_labels = labels,
                      target_distance = distance, target_edges = edges,
                      gauge_weight = gauge_weight, truth = truth,
                      truth_relative = truth - sum(gauge_weight * truth),
                      truth_pi0 = latent$pop0 / sum(latent$pop0),
                      truth_pi1 = latent$pop1 / sum(latent$pop1),
                      latent_population_t0 = latent$pop0, latent_population_t1 = latent$pop1,
                      config = config)
      dataset_file <- file.path(out, paste0(dataset_id, ".rds"))
      saveRDS(dataset, dataset_file)
      run_manifest[[length(run_manifest) + 1L]] <- data.frame(
        dataset_id = dataset_id, family = landscape$family, landscape_id = landscape$landscape_id,
        split = landscape$split, abm_seed = abm_seed, sampling_seed = sampling_seed,
        generating_protocol = config$generating_protocol,
        observed_nodes = length(observed), target_nodes = length(labels), max_distance = max(distance),
        dataset_file = normalizePath(dataset_file), stringsAsFactors = FALSE)
    }
  }
}
utils::write.csv(do.call(rbind, run_manifest), file.path(out, "dataset_manifest.csv"), row.names = FALSE)
cat(sprintf("Prepared %d immutable paired datasets in %s\n", length(run_manifest), normalizePath(out)))
