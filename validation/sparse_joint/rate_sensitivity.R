# Sensitivity of the sparse joint model to the two rates it is given: the
# baseline birth rate and the missegregation rate. Neither is known for a
# patient sample.
#
#   Rscript validation/sparse_joint/rate_sensitivity.R [--scenarios=REGEX] [--grid=typical|ends]
#                                                      [--workers=2] [--arguments='list(...)']
#
# `--grid=typical` runs the full grid; `--grid=ends` only the extremes, for
# use on every scenario. Each fit is cached under results/sweeps/rate.
# Set STT_R_LIB to the library holding the alfak2 build under test.

library_dir <- Sys.getenv("STT_R_LIB", unset = "")
if (nzchar(library_dir)) .libPaths(c(library_dir, .libPaths()))
for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R", "methods.R", "report.R")) {
  source(file.path("benchmark/sparse_two_timepoint", part))
}
args <- commandArgs(TRUE)
option <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else default
}
config <- stt_config()
index <- utils::read.csv(file.path(config$output_dir, "dataset_index_dev.csv"), stringsAsFactors = FALSE)
index <- index[index$scoreboard & grepl(option("scenarios", config$cohort_scenario), index$scenario), ]
extra <- eval(parse(text = option("arguments", "list()")))
tag <- option("tag", "default")

# Factors multiply the benchmark's settings: one division per day assumed
# (the simulated populations divide about 0.9 times), and the true
# missegregation rate.
grid <- if (identical(option("grid", "typical"), "ends")) {
  data.frame(birth = c(0.25, 4, 1, 1), beta = c(1, 1, 0.2, 5))
} else {
  rbind(data.frame(birth = c(0.1, 0.25, 0.5, 1, 2, 4), beta = 1),
        data.frame(birth = 1, beta = c(0.2, 0.5, 2, 5)),
        # Same flux as the reference: only the product should matter.
        data.frame(birth = c(2, 0.5), beta = c(0.5, 2)))
}
out_dir <- file.path(config$output_dir, "sweeps", "rate")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tasks <- merge(index["dataset_id"], grid)
tasks$path <- file.path(out_dir, sprintf("%s__%s_b%g_m%g.rds", tasks$dataset_id, tag, tasks$birth, tasks$beta))

run_task <- function(i) {
  task <- tasks[i, ]
  if (file.exists(task$path)) return(TRUE)
  dataset <- readRDS(stt_dataset_path(config, task$dataset_id))
  nodes <- stt_dataset_nodes_cached(task$dataset_id)
  truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(suppressWarnings(do.call(alfak2::fit_sparse_joint, c(
    list(counts, dt = dataset$interval, beta = dataset$beta * task$beta, birth_rate = task$birth,
         min_cn = 1L, max_cn = 8L), extra))), error = function(e) conditionMessage(e))
  result <- list(dataset_id = task$dataset_id, birth = task$birth, beta = task$beta,
                 seconds = proc.time()[["elapsed"]] - started)
  if (is.character(fit)) {
    result$error <- fit
  } else {
    prediction <- stats::predict(fit, nodes$karyotype)
    scored <- stt_score_predictions(
      nodes, truth,
      data.frame(karyotype = prediction$karyotype, fitness = prediction$fitness_mean,
                 fitness_sd = prediction$fitness_sd, stringsAsFactors = FALSE),
      config$scoring$reference_min_count)$metrics
    result$metrics <- scored[scored$stratum %in% c("d0", "d1", "d2"), ]
    result$objective <- fit$diagnostics$objective
    result$convergence <- fit$diagnostics$convergence
    result$scales <- fit$trend$scales
    result$shell <- unname(fit$trend$coefficients[length(fit$trend$coefficients)])
  }
  saveRDS(result, task$path)
  TRUE
}
stt_dataset_nodes_cached <- function(dataset_id) readRDS(file.path(config$output_dir, "nodes", paste0(dataset_id, ".rds")))

workers <- as.integer(option("workers", "2"))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
todo <- which(!file.exists(tasks$path))
cat("Fits to run:", length(todo), "of", nrow(tasks), "\n")
invisible(parallel::mclapply(todo, run_task, mc.cores = workers, mc.preschedule = FALSE))

results <- lapply(tasks$path[file.exists(tasks$path)], readRDS)
failed <- vapply(results, function(r) !is.null(r$error), logical(1))
cat("Fits:", length(results), " failed:", sum(failed), "\n")
rows <- do.call(rbind, lapply(results[!failed], function(r) {
  cbind(dataset_id = r$dataset_id, birth = r$birth, beta = r$beta, objective = as.numeric(r$objective),
        convergence = as.integer(r$convergence), r$metrics[c("stratum", "rmse", "rmse_zero", "bias", "truth_mean",
                                                 "slope", "spearman", "coverage95")])
}))
rows <- merge(rows, index[c("dataset_id", "scenario", "landscape_id")])
# Mean predicted fitness relative to the reference, per day and per division.
rows$predicted_mean <- rows$truth_mean + rows$bias
rows$predicted_per_division <- rows$predicted_mean / rows$birth
columns <- c("rmse", "rmse_zero", "bias", "predicted_mean", "predicted_per_division", "slope", "spearman", "coverage95")
table <- stats::aggregate(rows[columns], rows[c("birth", "beta", "stratum")], function(v) mean(v, na.rm = TRUE))
table[columns] <- round(table[columns], 3)
options(width = 200)
cat("\n== mean over datasets (", length(unique(rows$dataset_id)), ") ==\n")
print(table[order(table$stratum, table$beta, table$birth), ], row.names = FALSE)

# Does the marginal likelihood locate the birth rate? Lower objective is better.
profile <- rows[rows$stratum == "d0" & rows$beta == 1, c("dataset_id", "landscape_id", "birth", "objective")]
if (length(unique(profile$birth)) > 2L) {
  wide <- stats::reshape(profile[c("dataset_id", "birth", "objective")], idvar = "dataset_id", timevar = "birth",
                         direction = "wide")
  objective <- as.matrix(wide[-1L])
  colnames(objective) <- sub("objective.", "", colnames(objective), fixed = TRUE)
  objective <- objective[, order(as.numeric(colnames(objective))), drop = FALSE]
  cat("\n== marginal negative log-likelihood minus its minimum, by assumed birth rate ==\n")
  print(round(cbind(objective - apply(objective, 1L, min, na.rm = TRUE)), 1))
  cat("best rate per dataset:", colnames(objective)[apply(objective, 1L, which.min)], "\n")
  cat("summed over datasets: ", paste(colnames(objective), round(colSums(objective) - min(colSums(objective)), 1),
                                      sep = "=", collapse = "  "), "\n")
}
saveRDS(rows, file.path(config$output_dir, "sweeps", sprintf("rate_%s.rds", tag)))
