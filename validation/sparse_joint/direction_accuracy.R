# Direction accuracy of the main methods, from the saved benchmark fits
# (Section 5.1 of doc/alfak2_method.tex).
#
#   Rscript validation/sparse_joint/direction_accuracy.R [--workers=3]
#
# For each pair of neighbouring karyotypes x and x + e_c (one copy of one
# chromosome apart) among a dataset's scored nodes, the sign of the predicted
# change in relative fitness is compared with the sign of the true change. A
# correct sign scores 1, a wrong one 0 and a tie (equal predictions) 1/2, so
# that the null predictor and random guessing both score 1/2. Pairs with equal
# true fitness are skipped. A pair belongs to distance d when its farther
# member is at distance d from the observed karyotypes. Predictions and
# reference values are those of stt_score_predictions() (missing predictions
# at the reference level); centring cancels in the differences.
#
# Distances 0-2 use the scoreboard's nodes and fits. Distances 3-4 use the
# nodes and fits of far_distance.R, including the intermediate nodes of the
# sampling walks, so that every karyotype at distance 3 or 4 has neighbours one
# step closer. The score stage's edge_sign_accuracy counts ties as wrong and is
# not used here.
#
# Writes benchmark/results/sparse_two_timepoint/sweeps/direction/
# direction_<split>.csv, one row per dataset, method and distance (columns
# dataset_id, method, stratum, n_pairs, direction). A split whose file exists is
# not recomputed. One split takes about ten minutes with three workers.

for (part in c("landscape.R", "config.R", "score.R")) source(file.path("benchmark/sparse_two_timepoint", part))
args <- commandArgs(TRUE)
workers <- as.integer(sub("^--workers=", "", c(grep("^--workers=", args, value = TRUE), "--workers=3")[[1L]]))
config <- stt_config()
res <- config$output_dir
out <- file.path(res, "sweeps", "direction")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
methods <- c("alfak2_sparse_joint_balance", "alfak2_sparse_joint", "alfak2_two_layer_fixed", "alfakR_minobs20")

direction <- function(nodes, truth, prediction, distances) {
  table <- stt_score_predictions(nodes, truth, prediction, config$scoring$reference_min_count)$nodes
  cn <- stt_parse_karyotypes(table$karyotype)
  index_of <- stats::setNames(seq_len(nrow(table)), table$karyotype)
  edges <- do.call(rbind, lapply(seq_len(ncol(cn)), function(chr) {
    moved <- cn
    moved[, chr] <- moved[, chr] + 1
    child <- index_of[stt_format_karyotypes(moved)]
    keep <- which(!is.na(child))
    cbind(keep, unname(child[keep]))
  }))
  d <- pmax(table$distance[edges[, 1]], table$distance[edges[, 2]])
  dt <- table$truth_relative[edges[, 2]] - table$truth_relative[edges[, 1]]
  dp <- table$pred_relative[edges[, 2]] - table$pred_relative[edges[, 1]]
  dp[abs(dp) < 1e-12] <- 0
  score <- ifelse(dp == 0, 0.5, as.numeric(sign(dp) == sign(dt)))
  informative <- dt != 0
  do.call(rbind, lapply(distances, function(k) {
    use <- informative & d == k
    data.frame(stratum = paste0("d", k), n_pairs = sum(use), direction = if (any(use)) mean(score[use]) else NA_real_)
  }))
}

prediction_of <- function(path) {
  empty <- data.frame(karyotype = character(), fitness = numeric(), stringsAsFactors = FALSE)
  if (!file.exists(path)) return(empty)
  p <- readRDS(path)$prediction
  if (is.null(p)) empty else p
}

one_dataset <- function(id) {
  dataset <- readRDS(file.path(res, "datasets", paste0(id, ".rds")))
  near <- readRDS(file.path(res, "nodes", paste0(id, ".rds")))
  far <- readRDS(file.path(res, "sweeps", "far", "nodes", paste0(id, ".rds")))
  truth_near <- stt_true_growth_rate(near$karyotype, dataset$truth)
  truth_far <- stt_true_growth_rate(far$karyotype, dataset$truth)
  do.call(rbind, lapply(methods, function(m) {
    a <- direction(near, truth_near, prediction_of(file.path(res, "fits", sprintf("%s__%s.rds", id, m))), 0:2)
    b <- direction(far, truth_far, prediction_of(file.path(res, "sweeps", "far", "fits", sprintf("%s__%s.rds", id, m))), 3:4)
    cbind(dataset_id = id, method = m, rbind(a, b), stringsAsFactors = FALSE)
  }))
}

for (split in c("heldout", "dev")) {
  file <- file.path(out, sprintf("direction_%s.csv", split))
  if (file.exists(file)) {
    cat("Exists:", file, "\n")
    next
  }
  index <- utils::read.csv(file.path(res, sprintf("dataset_index_%s.csv", split)), stringsAsFactors = FALSE)
  ids <- index$dataset_id[index$scoreboard & index$role == "core"]
  rows <- parallel::mclapply(ids, one_dataset, mc.cores = workers, mc.preschedule = FALSE)
  failed <- !vapply(rows, is.data.frame, logical(1))
  if (any(failed)) stop("Failed datasets: ", paste(ids[failed], collapse = ", "))
  utils::write.csv(do.call(rbind, rows), file, row.names = FALSE)
  cat("Wrote", file, "\n")
}
