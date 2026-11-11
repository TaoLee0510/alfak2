# Calibration of a method's uncertainty on the benchmark: how large its errors
# are next to the standard deviations it reports, by kind of karyotype and by
# size of the predicted deficit.
#
#   Rscript validation/sparse_joint/calibration.R <method> [split] [scenario regex]

for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R", "methods.R", "report.R")) {
  source(file.path("benchmark/sparse_two_timepoint", part))
}
args <- commandArgs(TRUE)
method <- args[[1L]]
split <- if (length(args) > 1L) args[[2L]] else "dev"
config <- stt_config()
index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)), stringsAsFactors = FALSE)
index <- index[index$scoreboard & index$role == "core", ]
if (length(args) > 2L) index <- index[grepl(args[[3L]], index$scenario), ]

nodes <- do.call(rbind, lapply(index$dataset_id, function(id) {
  path <- file.path(config$output_dir, "fits", sprintf("%s__%s.rds", id, method))
  if (!file.exists(path)) return(NULL)
  fit <- readRDS(path)
  if (is.null(fit$prediction)) return(NULL)
  dataset <- readRDS(stt_dataset_path(config, id))
  scored <- readRDS(file.path(config$output_dir, "nodes", paste0(id, ".rds")))
  truth <- stt_true_growth_rate(scored$karyotype, dataset$truth)
  table <- stt_score_predictions(scored, truth, fit$prediction, config$scoring$reference_min_count)$nodes
  cbind(dataset_id = id, scenario = dataset$scenario, table)
}))
nodes <- nodes[is.finite(nodes$pred_sd), ]
nodes$error <- nodes$pred_relative - nodes$truth_relative
nodes$z <- nodes$error / nodes$pred_sd
total <- nodes$count_t0 + nodes$count_t1
nodes$kind <- ifelse(nodes$distance == 0L, ifelse(total >= 5, "observed, 5+ cells", "observed, <5 cells"),
                     paste("unobserved, distance", nodes$distance))

# Datasets weigh equally: summarise within a dataset, then average.
summarise <- function(part, by) {
  per_dataset <- do.call(rbind, lapply(split(part, list(part$dataset_id, part[[by]]), drop = TRUE), function(x) {
    data.frame(group = x[[by]][1L], n = nrow(x), rmse = sqrt(mean(x$error^2)), reported_sd = sqrt(mean(x$pred_sd^2)),
               z_sd = sqrt(mean(x$z^2)), coverage95 = mean(abs(x$z) <= stats::qnorm(0.975)),
               mean_prediction = mean(x$pred_relative), mean_truth = mean(x$truth_relative))
  }))
  out <- stats::aggregate(per_dataset[-1L], per_dataset["group"], mean)
  out$datasets <- as.integer(table(per_dataset$group)[as.character(out$group)])
  out
}
options(width = 200)
cat("Method:", method, " datasets:", length(unique(nodes$dataset_id)), "\n\n== by kind of karyotype ==\n")
print(format(summarise(nodes, "kind"), digits = 3), row.names = FALSE)

# Does the error grow with the predicted deficit? Unobserved karyotypes only.
unobserved <- nodes[nodes$distance > 0L, ]
unobserved$deficit <- cut(-unobserved$pred_relative, c(-Inf, 0.01, 0.02, 0.04, 0.08, 0.16, Inf))
cat("\n== unobserved karyotypes, by predicted deficit ==\n")
print(format(summarise(unobserved, "deficit"), digits = 3), row.names = FALSE)
cat("\n== by scenario ==\n")
print(format(summarise(nodes, "scenario"), digits = 3), row.names = FALSE)
saveRDS(nodes, file.path(config$output_dir, "report", sprintf("calibration_nodes_%s_%s.rds", method, split)))
