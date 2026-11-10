source("validation/accuracy_round2/pooled_calibration.R")
source("validation/accuracy_upgrade/interval_calibration.R")
set.seed(61107)
data <- expand.grid(landscape = 1:6, replicate = 1:2, distance = 0:3, node = 1:5)
data$landscape_id <- paste0("L", data$landscape)
data$dataset_id <- paste0(data$landscape_id, "R", data$replicate)
data$split <- "calibration"
data$pred_relative <- rnorm(nrow(data), 0, .1)
data$sd_relative <- runif(nrow(data), .05, .2)
data$truth_relative <- data$pred_relative + rnorm(nrow(data), -.03, .15)
key <- list(model = "test", sha256 = "fixed")
scenario <- list(domain = "synthetic")
old <- fit_alfak2_interval_calibration(data, key, scenario)
same <- fit_round2_pooled_calibration(data, key, scenario, gamma_far = 1, kappa = 0)
stopifnot(max(abs(old$table$q_lower - same$table$q_lower)) < 1e-12,
          max(abs(old$table$q_upper - same$table$q_upper)) < 1e-12)
pooled <- fit_round2_pooled_calibration(data, key, scenario, gamma_far = .5, kappa = 6)
prediction <- apply_round2_pooled_calibration(data, pooled, key, scenario)
stopifnot(identical(prediction$pred_relative, data$pred_relative),
          identical(prediction$sd_relative, data$sd_relative),
          all(prediction$calibrated_low_95 <= prediction$calibrated_low_80),
          all(prediction$calibrated_low_80 <= prediction$calibrated_low_50),
          all(prediction$calibrated_high_95 >= prediction$calibrated_high_80),
          all(prediction$calibrated_high_80 >= prediction$calibrated_high_50))
unknown <- data[1L, ]; unknown$distance <- 99
unknown <- apply_round2_pooled_calibration(unknown, pooled, key, scenario)
stopifnot(is.na(unknown$calibrated_low_95))
expect_error <- function(expression) stopifnot(inherits(try(force(expression), silent = TRUE), "try-error"))
expect_error(apply_round2_pooled_calibration(data, pooled, list(model = "wrong"), scenario))
expect_error(apply_round2_pooled_calibration(data, pooled, key, list(domain = "wrong")))
data$split <- "test"
expect_error(fit_round2_pooled_calibration(data, key, scenario))
data$split <- "development"
selection <- select_round2_calibration(data, key, scenario, gammas = c(0, 1), kappas = c(0, 6))
stopifnot(selection$selected$all_distance_nonworse,
          identical(selection$selected$grid_id,
            selection$grid_summary$grid_id[which.min(ifelse(selection$grid_summary$all_distance_nonworse,
              selection$grid_summary$macro_interval_score, Inf))]))
cat("Pooled calibration: reference equivalence, isolation, nesting, fingerprints and point/SD preservation passed.\n")
