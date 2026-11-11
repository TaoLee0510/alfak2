# Compare methods on the datasets where the candidate has a fit.
for (f in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R", "methods.R", "report.R")) source(file.path("benchmark/sparse_two_timepoint", f))
# Usage: Rscript compare.R <candidate> [further methods to compare with ...]
args <- commandArgs(TRUE); candidate <- args[1]; others <- c("zero", "alfak2_two_layer_fixed", "alfakR_minobs20", args[-1])
config <- stt_config(); index <- read.csv(file.path(config$output_dir, "dataset_index_dev.csv"), stringsAsFactors = FALSE)
fit_path <- function(id, m) file.path(config$output_dir, "fits", sprintf("%s__%s.rds", id, m))
ids <- index$dataset_id[index$scoreboard & file.exists(fit_path(index$dataset_id, candidate))]
rows <- do.call(rbind, lapply(ids, function(id) {
  dataset <- readRDS(stt_dataset_path(config, id)); nodes <- readRDS(file.path(config$output_dir, "nodes", paste0(id, ".rds")))
  truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  do.call(rbind, lapply(c(candidate, others), function(m) { fit <- readRDS(fit_path(id, m))
    pred <- if (is.null(fit$prediction)) data.frame(karyotype = character(), fitness = numeric()) else fit$prediction
    cbind(dataset_id = id, method = m, status = fit$status, seconds = fit$seconds, detail = substr(as.character(fit$detail), 1, 90), stt_score_predictions(nodes, truth, pred)$metrics) })) }))
rows <- merge(rows, index[c("dataset_id", "landscape_id", "wavelength_class", "scenario")])
saveRDS(rows, file.path(config$output_dir, "report", paste0("compare_", candidate, ".rds")))
cand <- rows[rows$method == candidate & rows$stratum == "d0", ]
cat("candidate fits:", nrow(cand), " status:", paste(names(table(cand$status)), table(cand$status), collapse = " "), " median seconds:", round(median(cand$seconds)), " max:", round(max(cand$seconds)), "\n")
if (any(cand$status != "ok")) print(cand[cand$status != "ok", c("dataset_id", "detail")], row.names = FALSE)
main <- rows[rows$stratum %in% c("d0", "d1", "d2"), ]
cols <- c("rmse", "rmse_zero", "truth_sd", "bias", "slope", "spearman", "auc_beneficial", "edge_sign_accuracy", "coverage95")
per_land <- aggregate(main[cols], main[c("landscape_id", "stratum", "method")], function(v) mean(v, na.rm = TRUE))
agg <- aggregate(per_land[cols], per_land[c("stratum", "method")], function(v) mean(v, na.rm = TRUE)); agg[cols] <- round(agg[cols], 3)
options(width = 200); cat("\n== mean over landscapes ==\n"); print(agg[order(agg$stratum, agg$method), ], row.names = FALSE)
cat("\n== candidate RMSE / zero RMSE by wavelength class ==\n")
cc <- main[main$method == candidate, ]; cc$ratio <- cc$rmse / cc$rmse_zero
print(round(tapply(cc$ratio, list(cc$wavelength_class, cc$stratum), mean), 3))
cat("\n== candidate Spearman by wavelength class (current two-layer in brackets) ==\n")
tl <- main[main$method == "alfak2_two_layer_fixed", ]
a <- tapply(cc$spearman, list(cc$wavelength_class, cc$stratum), mean, na.rm = TRUE); b <- tapply(tl$spearman, list(tl$wavelength_class, tl$stratum), mean, na.rm = TRUE)
print(noquote(matrix(sprintf("%.2f (%.2f)", a, b), nrow = nrow(a), dimnames = dimnames(a))))
cat("\n== paired RMSE difference, candidate minus other, landscape bootstrap ==\n")
set.seed(1); for (o in others) for (s in c("d0", "d1", "d2")) { x <- per_land[per_land$stratum == s, ]; a <- x$rmse[x$method == candidate][order(x$landscape_id[x$method == candidate])]; b <- x$rmse[x$method == o][order(x$landscape_id[x$method == o])]
  d <- a - b; boot <- replicate(2000, mean(sample(d, replace = TRUE))); cat(sprintf("%-24s %s: %+.4f [%+.4f, %+.4f] (%+.0f%%)\n", o, s, mean(d), quantile(boot, .025), quantile(boot, .975), 100 * mean(d) / mean(b))) }
