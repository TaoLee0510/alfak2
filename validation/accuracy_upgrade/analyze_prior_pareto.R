#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Usage: analyze_prior_pareto.R RESULTS OUTPUT_CSV")
files <- list.files(args[[1L]], "_metrics\\.csv$", full.names = TRUE)
files <- files[grepl("^(baseline|joint_endpoint_b[0-9])", basename(files))]
x <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
x <- x[x$distance_group %in% paste0("d", 0:3) & x$metric %in% c("relative_rmse", "interval_score_95", "coverage_95"), ]
a <- aggregate(value ~ model + distance_group + metric, x, mean)
a$key <- paste(a$distance_group, a$metric, sep = "_")
wide <- reshape(a[, c("model", "key", "value")], idvar = "model", timevar = "key", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
point <- paste0("d", 0:3, "_relative_rmse")
interval <- paste0("d", 0:3, "_interval_score_95")
base <- wide[wide$model == "baseline", ]
wide$all_distance_rmse_better_than_baseline <- apply(wide[, point], 1, function(z) all(z < unlist(base[, point], use.names = FALSE)))
wide$all_distance_interval_score_better_than_baseline <- apply(wide[, interval], 1, function(z) all(z < unlist(base[, interval], use.names = FALSE)))
v <- as.matrix(wide[, c(point, interval)])
wide$pareto_rmse_and_interval_score <- vapply(seq_len(nrow(v)), function(i) {
  !any(vapply(setdiff(seq_len(nrow(v)), i), function(j) all(v[j, ] <= v[i, ]) && any(v[j, ] < v[i, ]), logical(1)))
}, logical(1))
wide$max_relative_rmse_ratio <- apply(wide[, point], 1, function(z) max(z / unlist(base[, point], use.names = FALSE)))
wide$max_relative_interval_score_ratio <- apply(wide[, interval], 1, function(z) max(z / unlist(base[, interval], use.names = FALSE)))
coverage <- as.matrix(wide[, paste0("d", 0:3, "_coverage_95")])
wide$max_coverage_error_from_95 <- apply(abs(coverage - 0.95), 1, max)
dir.create(dirname(args[[2L]]), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(wide, args[[2L]], row.names = FALSE)
print(wide[, c("model", "all_distance_rmse_better_than_baseline", "all_distance_interval_score_better_than_baseline",
               "pareto_rmse_and_interval_score", "max_relative_rmse_ratio", "max_relative_interval_score_ratio", "max_coverage_error_from_95")], row.names = FALSE)
