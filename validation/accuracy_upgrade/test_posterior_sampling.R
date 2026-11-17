#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "posterior_sampling.R"))
labels <- c("target", "reference")
covariance <- diag(c(1, 0.01)); dimnames(covariance) <- list(labels, labels)
components <- lapply(c(-2, 2), function(mu) list(weight = 0.5,
  mean = stats::setNames(c(mu, 0), labels), covariance = covariance,
  forward_model = "selection_transport"))
draws <- draw_alfak2_posterior_components(components, labels, 80000L, 1234L)
contrast <- draws[1, ] - draws[2, ]
quantile <- stats::quantile(contrast, c(0.025, 0.975))
# The mixture has narrower outer quantiles than its moment-matched Gaussian.
stopifnot(abs(mean(contrast)) < 0.04, quantile[1] > -3.8, quantile[1] < -3.5,
          quantile[2] > 3.5, quantile[2] < 3.8,
          stats::qnorm(0.975) * sqrt(5.01) > quantile[2] + 0.5)
identical_draws <- draw_alfak2_posterior_components(components, labels, 80000L, 1234L)
stopifnot(identical(draws, identical_draws))
cat("PASS: actual Gaussian-mixture contrast quantiles, full reference propagation, deterministic draws\n")
