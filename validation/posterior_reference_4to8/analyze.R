#!/usr/bin/env Rscript
repo <- normalizePath(".")
directory <- file.path(repo, "validation/posterior_reference_4to8")
output <- Sys.getenv("ALFAK2_REFERENCE_OUTPUT", file.path(directory, "results"))
source(file.path(directory, "reference.R"))
source(file.path(repo, "validation/accuracy_round2/common.R"))
runtime <- reference_load(repo)
cases <- readRDS(file.path(output, "frozen_cases.rds"))
status <- read.csv(file.path(output, "case_status.csv"))
if (!setequal(status$case_id, names(cases))) stop("Incomplete primary case set; do not silently drop cases")
moments <- covariance <- quantiles <- forward <- scale_shape <- list()
for (ci in seq_along(cases)) {
  case <- cases[[ci]]; n <- case$n
  stored <- readRDS(file.path(output, paste0(case$id, "_summary.rds")))
  fit <- stored$laplace; d <- stored$reference
  actual <- readRDS(stored$raw_chain_file)
  mu_g <- fit$latent_fitness_mean; Cg <- fit$latent_fitness_covariance
  mu_b <- exp(mu_g + .5 * diag(Cg)); Cb <- tcrossprod(mu_b) * expm1(Cg)
  stopifnot(max(abs(mu_b - fit$summary$fitness_mean)) < 1e-9,
    max(abs(Cb - fit$fitness_covariance)) < 1e-9)
  # This decomposition holds the reference's first two log-rate moments fixed.
  # The remaining discrepancy measures the non-Gaussian exponential transform;
  # it is not attributed to a new prior or to information added by an oracle.
  rg <- d$summary[match(paste0("log_birth_", 1:n), d$summary$variable), ]
  reference_gaussian_b <- exp(rg$mean + .5 * diag(d$covariance$log_birth))
  reference_gaussian_Cb <- tcrossprod(reference_gaussian_b) * expm1(d$covariance$log_birth)
  E <- matrix(0, nrow(case$edges), n)
  E[cbind(1:nrow(E), case$edges[, 1])] <- -1
  E[cbind(1:nrow(E), case$edges[, 2])] <- 1
  transforms <- list(log_birth = diag(n), birth = diag(n),
    relative = diag(n) - matrix(case$weight, n, n, byrow = TRUE),
    centered = diag(n) - matrix(1 / n, n, n), edge = E)
  set.seed(730000L + ci)
  gd <- sweep(matrix(rnorm(200000L * n), ncol = n) %*% chol(Cg), 2L, mu_g, "+")
  bd <- exp(gd)
  obs <- lapply(actual$chains, function(ch) reference_observables(ch$draws, case))
  ref_log <- do.call(rbind, lapply(obs, function(x) x[, paste0("log_birth_", 1:n), drop = FALSE]))
  ref_birth <- exp(ref_log)
  alpha <- rowMeans(ref_log)
  cuts <- c(-Inf, as.numeric(quantile(alpha, c(.25, .5, .75))), Inf)
  # Post-hoc mechanism diagnostic on all cases, not a new validation gate.
  for (source in c("reference_MCMC", "Laplace_Gaussian")) {
    logs <- if (source == "reference_MCMC") ref_log else gd
    rates <- if (source == "reference_MCMC") ref_birth else bd
    a <- rowMeans(logs); band <- cut(a, cuts, labels = FALSE)
    for (bin in 1:4) {
      rows <- which(band == bin)
      log_edge <- logs[rows, case$edges[, 2], drop = FALSE] - logs[rows, case$edges[, 1], drop = FALSE]
      birth_edge <- rates[rows, case$edges[, 2], drop = FALSE] - rates[rows, case$edges[, 1], drop = FALSE]
      scale_shape[[length(scale_shape) + 1L]] <- data.frame(case_id = case$id,
        source = source, band = bin, n_draws = length(rows),
        common_log_scale_mean = mean(a[rows]),
        mean_log_edge_sd = mean(apply(log_edge, 2L, stats::sd)),
        mean_birth_edge_sd = mean(apply(birth_edge, 2L, stats::sd)),
        interpretation = "post-hoc binned variance; not exact conditional variance or a promotion gate")
    }
  }
  for (type in names(transforms)) {
    A <- transforms[[type]]
    lm <- if (type == "log_birth") as.numeric(mu_g) else as.numeric(A %*% mu_b)
    LC <- if (type == "log_birth") Cg else A %*% Cb %*% t(A)
    ref <- d$summary[match(paste0(type, "_", 1:nrow(A)), d$summary$variable), ]
    gm <- if (type == "log_birth") rg$mean else as.numeric(A %*% reference_gaussian_b)
    GC <- if (type == "log_birth") d$covariance$log_birth else A %*% reference_gaussian_Cb %*% t(A)
    difference <- lm - ref$mean
    sd_difference <- sqrt(pmax(0, diag(LC))) - ref$sd
    moments[[length(moments) + 1L]] <- data.frame(case_id = case$id, n = n, profile = case$profile,
      beta = case$beta, growth = case$growth, reference_certified = d$convergence$passed,
      type = type, index = 1:nrow(A), reference_mean = ref$mean, laplace_mean = lm,
      reference_sd = ref$sd, laplace_sd = sqrt(pmax(0, diag(LC))),
      mcse_mean = ref$mcse_mean, mcse_sd = ref$mcse_sd,
      standardized_mean_difference = difference / ref$sd,
      posterior_squared_loss_ratio = 1 + (difference / ref$sd)^2,
      sd_ratio = sqrt(pmax(0, diag(LC))) / ref$sd,
      mean_difference_in_mcse = difference / ref$mcse_mean,
      sd_difference_in_mcse = sd_difference / ref$mcse_sd,
      material_mean = abs(difference) > .05 * ref$sd & abs(difference) > 4 * ref$mcse_mean,
      material_sd = abs(sd_difference) > .1 * ref$sd & abs(sd_difference) > 4 * ref$mcse_sd,
      gaussian_from_reference_log_moments_mean = gm,
      gaussian_from_reference_log_moments_sd = sqrt(pmax(0, diag(GC))),
      latent_moment_part_mean_error = lm - gm,
      non_gaussian_transform_part_mean_error = gm - ref$mean)
    stopifnot(max(abs((lm - gm) + (gm - ref$mean) - difference)) < 1e-12,
      max(abs(sqrt(diag(d$covariance[[type]])) - ref$sd)) < 1e-10)
    RC <- d$covariance[[type]]
    covariance[[length(covariance) + 1L]] <- data.frame(case_id = case$id, type = type,
      reference_certified = d$convergence$passed,
      laplace_relative_frobenius_error = norm(LC - RC, "F") / norm(RC, "F"),
      reference_covariance_relative_mcse = d$covariance_relative_mcse[[type]],
      gaussian_from_reference_log_moments_relative_error = norm(GC - RC, "F") / norm(RC, "F"))
    lp_draw <- if (type == "log_birth") gd else bd %*% t(A)
    lpq <- t(apply(lp_draw, 2L, quantile, probs = c(.025, .975)))
    lpq_mcse <- t(apply(lp_draw, 2L, posterior::mcse_quantile, probs = c(.025, .975)))
    rq <- d$quantiles[match(ref$variable, rownames(d$quantiles)), c(1, 3), drop = FALSE]
    for (j in seq_len(nrow(A))) {
      values <- do.call(cbind, lapply(obs, function(x) x[, ref$variable[j]]))
      event <- values >= lpq[j, 1] & values <= lpq[j, 2]
      qse <- posterior::mcse_quantile(values, probs = c(.025, .975))
      quantiles[[length(quantiles) + 1L]] <- data.frame(case_id = case$id, type = type, index = j,
        reference_certified = d$convergence$passed,
        reference_q025 = rq[j, 1], reference_q975 = rq[j, 2],
        reference_q025_mcse = qse[1], reference_q975_mcse = qse[2],
        laplace_q025 = lpq[j, 1], laplace_q975 = lpq[j, 2],
        laplace_q025_iid_mcse = lpq_mcse[j, 1], laplace_q975_iid_mcse = lpq_mcse[j, 2],
        reference_posterior_mass_in_laplace_interval = mean(event),
        posterior_mass_mcse = posterior::mcse_mean(event * 1),
        interpretation = "fixed-data posterior mass, not repeated-sampling coverage")
    }
  }
  # Representative posterior states check numerical propagation independently.
  for (ch in seq_along(actual$chains)) {
    z <- actual$chains[[ch]]$draws
    chosen <- unique(round(seq(1, nrow(z), length.out = 32L)))
    for (j in chosen) {
      state <- z[j, ]; sp <- case$specification
      g <- as.numeric(sp$offset + sp$lower %*% state[1:n])
      eta <- as.numeric(sp$eta_basis %*% state[n + 1:(n - 1)])
      pi0 <- exp(eta - max(eta)); pi0 <- pi0 / sum(pi0)
      generator <- sweep(2 * t(sp$probability) - diag(n), 2L, exp(g), "*")
      mass <- runtime$independent$process_deterministic_mean(generator, pi0, 1.5)
      p1 <- mass / sum(mass); lg <- log(sum(mass))
      independent_ll <- sum(sp$y0 * log(pi0 + 1e-14)) + sum(sp$y1 * log(p1 + 1e-14))
      if (sp$growth_sd > 0) independent_ll <- independent_ll + dnorm(lg, sp$growth, sp$growth_sd, log = TRUE)
      numerical <- runtime$cpp$ref_likelihood_cpp(state, sp)
      forward[[length(forward) + 1L]] <- data.frame(case_id = case$id, chain = ch, draw = j,
        maximum_pi1_difference = max(abs(p1 - numerical$pi1)),
        log_growth_difference = lg - numerical$log_growth,
        exact_minus_split_log_likelihood = independent_ll - numerical$log_likelihood)
    }
  }
  cat("Compared", case$id, "\n")
}
write.csv(do.call(rbind, moments), file.path(output, "moment_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, covariance), file.path(output, "covariance_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, quantiles), file.path(output, "quantile_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, forward), file.path(output, "independent_forward_checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, scale_shape), file.path(output, "common_scale_shape_bins.csv"), row.names = FALSE)
extra <- c(file.path(directory, "analyze.R"),
  file.path(repo, "validation/accuracy_round2/process_prepare_data.R"),
  file.path(repo, "validation/accuracy_round2/oracle_posterior_reference.R"),
  file.path(output, "validation/validation_tests.csv"),
  file.path(directory, "inputs/two_node_reference.rds"))
write.csv(data.frame(file = extra, sha256 = unname(round2_sha256(extra))),
  file.path(output, "independent_analysis_sha256.csv"), row.names = FALSE)
cat("All comparison tables written.\n")
