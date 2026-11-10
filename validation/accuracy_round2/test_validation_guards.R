#!/usr/bin/env Rscript
# Regression tests for evaluator completeness, independent of native model fits.
validation_dir <- if (file.exists("validation/accuracy_round2/pooled_calibration.R"))
  "validation/accuracy_round2" else "."
source(file.path(validation_dir, "pooled_calibration.R"))
source(file.path(validation_dir, "joint_power.R"))
library(testthat)

test_that("power planning checks every predeclared gate and rejects duplicate rows", {
  x <- expand.grid(family = c("smooth", "rough", "interaction"), landscape = 1:4,
    distance = 0:1, metric = c("rmse", "interval_score95", "spatially_centered_rmse", "edge_gradient_rmse"),
    stringsAsFactors = FALSE)
  x$landscape_id <- paste(x$family, x$landscape)
  x$candidate_id <- "B_frozen"
  x$reference_loss <- 1 + x$landscape / 10
  x$candidate_loss <- .8 * x$reference_loss
  x$oof_verified <- TRUE
  expect_error(round2_joint_power(x, repetitions = 20L), "required_distances")
  expect_error(round2_joint_power(x[x$metric == "rmse", ], required_distances = 0:1,
    repetitions = 20L), "metric-by-distance")
  expect_error(round2_joint_power(x[x$distance == 0, ], required_distances = 0:1,
    repetitions = 20L), "metric-by-distance")
  expect_error(round2_joint_power(rbind(x, x[1, ]), required_distances = 0:1,
    repetitions = 20L), "unique nonmissing")
  y <- x; y$candidate_loss[1] <- NA_real_
  expect_error(round2_joint_power(y, required_distances = 0:1, repetitions = 20L), "Incomplete OOF")
  y <- x[-1, ]
  expect_error(round2_joint_power(y, required_distances = 0:1, repetitions = 20L), "Missing OOF gate vectors")
  a <- round2_joint_power(x, n_total = 36L, required_distances = 0:1, repetitions = 20L)
  b <- round2_joint_power(x[nrow(x):1, ], n_total = 36L, required_distances = 0:1, repetitions = 20L)
  expect_equal(a$n_gates, 32L)
  expect_equal(a, b)
  expect_false(a$full_promotion_validated)
  expect_true("approximate_core_gate_pass_probability" %in% names(a))
  expect_false("approximate_joint_pass_probability" %in% names(a))
})

test_that("calibration selection cannot silently discard unsupported distances", {
  x <- expand.grid(landscape = 1:5, distance = c(0, 2), node = 1:2)
  x$landscape_id <- paste0("L", x$landscape)
  x$dataset_id <- paste0("D", x$landscape)
  x$split <- "development"
  x$pred_relative <- 0
  x$sd_relative <- .1
  x$truth_relative <- .02 * x$landscape * ifelse(x$node == 1, -1, 1)
  rare <- x[x$distance == 0 | x$landscape <= 2, ]
  expect_error(select_round2_calibration(rare, "v1", "test", gammas = 1, kappas = 0,
    required_distances = c(0, 2)), "at least three")
  expect_error(select_round2_calibration(x[x$distance == 0, ], "v1", "test",
    gammas = 1, kappas = 0, required_distances = c(0, 2)), "do not match")
  invalid <- x; invalid$sd_relative[1] <- NA_real_
  expect_error(select_round2_calibration(invalid, "v1", "test", gammas = 1, kappas = 0), "Incomplete or invalid")
  a <- select_round2_calibration(x, "v1", "test", gammas = 1, kappas = 0,
    required_distances = c(0, 2))
  expect_setequal(unique(a$heldout_scores$distance), c(0, 2))
  expect_equal(nrow(a$distance_risk), 6L)
  expect_equal(unname(a$landscapes_per_distance), c(5L, 5L))
  expect_equal(a$distance_scope, "explicit_predeclared_distances")
})
