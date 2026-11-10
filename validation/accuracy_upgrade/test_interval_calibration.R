#!/usr/bin/env Rscript
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
source(file.path(script_dir, "interval_calibration.R"))
data <- data.frame(landscape_id = rep(c("L1", "L2", "L3"), each = 4),
                   dataset_id = rep(c("A", "B", "C"), each = 4), split = "calibration",
                   distance = 0L, pred_relative = 0, sd_relative = 1,
                   truth_relative = rep(c(-2, -1, 1, 3), 3))
model_key <- list(model = "frozen", implementation_hash = "source+dll", arguments_hash = "prior")
scenario <- list(beta = 0.005, depth_t0 = 150L)
fit <- fit_alfak2_interval_calibration(data, model_key, scenario)
duplicated_nodes <- rbind(data, data[data$landscape_id == "L1", ])
fit2 <- fit_alfak2_interval_calibration(duplicated_nodes, model_key, scenario)
stopifnot(identical(fit$table$q_lower, fit2$table$q_lower), identical(fit$table$q_upper, fit2$table$q_upper))
prediction <- data.frame(distance = c(0L, 5L), pred_relative = c(0.5, -0.1), sd_relative = c(2, 1))
applied <- apply_alfak2_interval_calibration(prediction, fit, model_key, scenario)
stopifnot(identical(applied$pred_relative, prediction$pred_relative), identical(applied$sd_relative, prediction$sd_relative),
          all(is.na(applied[2, grep("^calibrated_", names(applied))])), applied$calibration_status[2] == "unvalidated_distance")
must_error <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")
leaked <- data; leaked$split[1] <- "test"
stopifnot(must_error(fit_alfak2_interval_calibration(leaked, model_key, scenario)),
          must_error(apply_alfak2_interval_calibration(prediction, fit, model_key, list(beta = 0.005, depth_t0 = 1000L))),
          must_error(apply_alfak2_interval_calibration(prediction, fit, list(model = "changed"), scenario)))
cat("PASS: calibration split isolation, landscape balancing, unchanged point/SD, unknown distance, strict model/scenario matching\n")
