round2_lattice <- function(domain) {
  k <- as.matrix(expand.grid(rep(list(domain$min_cn:domain$max_cn), domain$n_chr)))
  storage.mode(k) <- "integer"
  labels <- apply(k, 1L, paste, collapse = ".")
  rownames(k) <- labels
  list(karyotypes = k, labels = labels,
       distance_matrix = as.matrix(stats::dist(k, method = "manhattan")))
}

round2_truth <- function(domain_name, domain, family, seed, lattice) {
  set.seed(seed)
  k <- lattice$karyotypes; x <- k - 2L; p <- ncol(k)
  founder <- which(rowSums(abs(x)) == 0)
  parameters <- list()
  if (family == "smooth") {
    a <- stats::runif(p, -0.16, 0.16)
    g <- as.numeric(x %*% a) - 0.06 * rowSums(x^2)
    parameters <- list(dosage = a, curvature = -0.06)
  } else if (family == "rough") {
    centroids <- matrix(stats::runif(16L * p, -5, 8), ncol = p)
    raw <- numeric(nrow(k))
    for (j in seq_len(nrow(centroids))) {
      d <- sqrt(rowSums(sweep(k, 2L, centroids[j, ], "-")^2))
      raw <- raw + sin(d / 0.22)
    }
    g <- 0.18 * (raw - raw[founder]) / stats::sd(raw)
    parameters <- list(centroids = centroids, wavelength = 0.22, scale = 0.18)
  } else if (family == "interaction") {
    co <- stats::runif(3, 0.12, 0.25) * sample(c(-1, 1), 3, replace = TRUE)
    second_pair <- if (p == 4L) x[, 3] * x[, 4] else x[, 2] * x[, 3]
    g <- 0.06 * x[, 1] + co[1] * x[, 1] * x[, 2] + co[2] * second_pair + co[3] * x[, 1] * x[, 2] * x[, 3]
    parameters <- list(first_order = 0.06, interaction = co)
  } else stop("Unknown generating family")
  if (domain_name == "cn4_1to3") {
    birth <- 1 + g
    version <- "round1_linear_truth_unchanged"
  } else {
    saturation <- log(2)
    birth <- exp(saturation * tanh(g / saturation))
    version <- "bounded_log_link_v1"
    parameters$saturation <- saturation
  }
  if (any(!is.finite(birth)) || any(birth <= 0)) stop("INVALID_PREDECLARED_TRUTH: no seed replacement or clipping allowed")
  stopifnot(birth[founder] == 1)
  list(birth = as.numeric(birth), family = family, parameters = parameters, version = version)
}

round2_prepare_dataset <- function(config, domain_name, landscape, abm_replicate, sampling_replicate,
                                   profile = "standard", latent = NULL, truth = NULL) {
  domain <- config$domains[[domain_name]]
  lattice <- round2_lattice(domain)
  k <- lattice$karyotypes; labels <- lattice$labels; distance_matrix <- lattice$distance_matrix
  founder <- which(rowSums(abs(k - 2L)) == 0)
  fixed_stencil <- which(distance_matrix[founder, ] <= 1L)
  fixed_distance <- apply(distance_matrix[, fixed_stencil, drop = FALSE], 1L, min)
  local_config <- config
  local_config$n_chr <- domain$n_chr; local_config$min_cn <- domain$min_cn; local_config$max_cn <- domain$max_cn
  local_config$max_nodes <- domain$n_target; local_config$domain <- domain_name
  if (profile == "support_limit") { local_config$beta <- 0; local_config$initial_neighbour_mass <- 0 }
  if (is.null(truth)) truth <- round2_truth(domain_name, domain, landscape$family, landscape$landscape_seed, lattice)
  abm_seed <- as.integer(landscape$landscape_seed * 10L + abm_replicate)
  sampling_seed <- as.integer(abm_seed * 10L + sampling_replicate)
  if (is.null(latent)) {
    p0 <- numeric(nrow(k)); neighbors <- which(distance_matrix[founder, ] == 1L)
    p0[founder] <- 1 - local_config$initial_neighbour_mass
    p0[neighbors] <- local_config$initial_neighbour_mass / length(neighbors)
    latent <- simulate_exact_division_abm(k, truth$birth, local_config$beta, local_config$dt, p0,
                                          local_config$initial_population, abm_seed)
  }
  set.seed(sampling_seed)
  y0 <- as.numeric(stats::rmultinom(1, local_config$n0, latent$pop0))
  y1 <- as.numeric(stats::rmultinom(1, local_config$n1, latent$pop1))
  counts <- cbind(t0 = y0, t1 = y1); rownames(counts) <- labels
  observed <- which(rowSums(counts) > 0); counts <- counts[observed, , drop = FALSE]
  actual_distance <- apply(distance_matrix[, observed, drop = FALSE], 1L, min)
  weights <- numeric(length(labels)); weights[observed] <- counts[, "t0"] / sum(counts[, "t0"])
  edges <- which(distance_matrix == 1 & upper.tri(distance_matrix), arr.ind = TRUE)
  id <- sprintf("%s_%s_A%d_S%d", profile, landscape$landscape_id, abm_replicate, sampling_replicate)
  list(dataset_id = id, domain = domain_name, profile = profile, family = landscape$family,
    landscape_id = landscape$landscape_id, split = landscape$split,
    landscape_seed = landscape$landscape_seed, abm_seed = abm_seed, sampling_seed = sampling_seed,
    counts = counts, target_labels = labels, target_karyotypes = k, target_edges = edges,
    target_distance = actual_distance, actual_target_distance = actual_distance,
    auxiliary_stencil_distance = fixed_distance, auxiliary_stencil_labels = labels[fixed_stencil],
    gauge_weight = weights, truth = truth$birth, truth_relative = truth$birth - sum(weights * truth$birth),
    truth_specification = truth, truth_pi0 = latent$pop0 / sum(latent$pop0), truth_pi1 = latent$pop1 / sum(latent$pop1),
    latent_population_t0 = latent$pop0, latent_population_t1 = latent$pop1,
    config = local_config)
}
