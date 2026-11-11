# How much of the error is sampling, and how much is the model?
#
#   Rscript validation/sparse_joint/depth_limit.R [--depth=20000] [--split=dev] [--workers=1]
#
# For the typical scenario, the two samples are redrawn from the same true
# populations at a much larger depth and the sparse joint model is fitted to
# them. Two references are scored beside it:
#   shell mean   every karyotype predicted at the true mean fitness of its
#                distance stratum (0, 1, 2): what a perfect trend without
#                karyotype detail would give;
#   deep fit     the model given (nearly) the whole population.
# The fit at the benchmark depth is read from the benchmark's fits.
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
depth <- as.numeric(option("depth", "20000"))
split <- option("split", "dev")
workers <- as.integer(option("workers", "1"))
config <- stt_config()
index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)), stringsAsFactors = FALSE)
ids <- index$dataset_id[index$scoreboard & index$scenario == config$cohort_scenario]
out_dir <- file.path(config$output_dir, "sweeps", "depth")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

score_one <- function(id) {
  path <- file.path(out_dir, sprintf("%s__depth%g.rds", id, depth))
  if (file.exists(path)) return(readRDS(path))
  dataset <- readRDS(stt_dataset_path(config, id))
  nodes <- readRDS(file.path(config$output_dir, "nodes", paste0(id, ".rds")))
  truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  score <- function(prediction, name) {
    metrics <- stt_score_predictions(nodes, truth, prediction, config$scoring$reference_min_count)$metrics
    cbind(dataset_id = id, method = name, metrics[metrics$stratum %in% c("d0", "d1", "d2"),
                                                  c("stratum", "rmse", "rmse_zero", "bias", "spearman", "coverage95")])
  }
  # Shell mean: the truth relative to the scoring reference, averaged within
  # stratum. The reference karyotypes form a group of their own, predicted at
  # their weighted mean, zero, so that the scorer's recentring leaves the
  # prediction unchanged.
  reference <- stt_reference_weights(nodes, config$scoring$reference_min_count)
  relative <- truth - sum(reference * truth)
  group <- ifelse(reference > 0, "reference", nodes$distance)
  shell_mean <- stats::ave(relative, group)
  shell_mean[reference > 0] <- 0
  # Deep samples from the same populations; the scored nodes stay those of the
  # benchmark dataset.
  deep <- stt_sample_counts(dataset$truth$frequency_t0, dataset$truth$frequency_t1, depth, depth,
                            dataset$interval, dataset$seed + 7L)
  storage.mode(deep) <- "integer"
  method <- config$methods$alfak2_sparse_joint
  started <- proc.time()[["elapsed"]]
  fit <- alfak2::fit_sparse_joint(deep, dt = dataset$interval, beta = dataset$beta, birth_rate = method$birth_rate,
                                  min_cn = method$min_cn, max_cn = method$max_cn)
  seconds <- proc.time()[["elapsed"]] - started
  prediction <- stats::predict(fit, nodes$karyotype)
  benchmark_fit <- readRDS(file.path(config$output_dir, "fits", sprintf("%s__alfak2_sparse_joint.rds", id)))
  rows <- rbind(
    score(data.frame(karyotype = nodes$karyotype, fitness = shell_mean), "shell mean"),
    score(benchmark_fit$prediction, "benchmark depth"),
    score(data.frame(karyotype = prediction$karyotype, fitness = prediction$relative_fitness,
                     fitness_sd = prediction$relative_sd), sprintf("depth %g", depth))
  )
  rows$seconds <- seconds
  rows$states <- fit$diagnostics$n_states
  saveRDS(rows, path)
  rows
}
rows <- parallel::mclapply(ids, function(id) tryCatch(score_one(id), error = function(e) conditionMessage(e)),
                           mc.cores = workers, mc.preschedule = FALSE)
failed <- vapply(rows, is.character, logical(1))
if (any(failed)) cat("Failed:", paste(ids[failed], unlist(rows[failed]), collapse = "; "), "\n")
rows <- do.call(rbind, rows[!failed])
columns <- c("rmse", "rmse_zero", "bias", "spearman", "coverage95")
table <- stats::aggregate(rows[columns], rows[c("stratum", "method")], function(v) mean(v, na.rm = TRUE))
table[columns] <- round(table[columns], 3)
options(width = 160)
cat(sprintf("%s, %d datasets; deep samples of %g cells per timepoint (median fit %.0f s, %d states)\n",
            config$cohort_scenario, length(unique(rows$dataset_id)), depth, stats::median(rows$seconds),
            as.integer(stats::median(rows$states))))
print(table[order(table$stratum, table$method), ], row.names = FALSE)
