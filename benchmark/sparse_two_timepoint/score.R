# Scoring for the sparse two-timepoint benchmark.
#
# Every method is scored on the same nodes, which are defined from the counts
# alone: observed karyotypes (distance 0) and their one- and two-step lattice
# neighbourhoods. Truth and predictions are compared as relative fitness, each
# centred with the same reference weights, so a method's absolute scale and
# zero point do not enter.

stt_one_step_neighbours <- function(karyotypes, min_cn, max_cn) {
  karyotypes <- as.matrix(karyotypes)
  n_chr <- ncol(karyotypes)
  out <- vector("list", 2L * n_chr)
  for (chr in seq_len(n_chr)) {
    for (direction in c(-1, 1)) {
      moved <- karyotypes
      moved[, chr] <- moved[, chr] + direction
      keep <- moved[, chr] >= min_cn & moved[, chr] <= max_cn
      out[[2L * (chr - 1L) + (direction + 3L) / 2L]] <- moved[keep, , drop = FALSE]
    }
  }
  out <- do.call(rbind, out)
  out[!duplicated(stt_format_karyotypes(out)), , drop = FALSE]
}

# Karyotypes three and four steps from the observed ones, for the
# extrapolation check (validation/sparse_joint/far_distance.R). The shells are
# far too large to enumerate at 22 chromosomes, so they are sampled by random
# walks of three or four single-copy steps from observed karyotypes, each
# observed karyotype being an equally likely start, keeping endpoints whose
# distance to the nearest observed karyotype is exactly the number of steps.
# Within a shell this sampling is not uniform. The intermediate nodes of the
# kept walks are returned as unscored rows (`scored = FALSE`): graph-based
# methods need a path from the observed karyotypes to a target. The observed
# karyotypes are included as the reference for centring.
stt_far_nodes <- function(counts, min_cn = 1L, max_cn = 8L, cap = 2000L, seed = 1L, max_batches = 60L) {
  counts <- as.matrix(counts)
  observed <- rownames(counts)[rowSums(counts) > 0]
  obs <- stt_parse_karyotypes(observed)
  n_chr <- ncol(obs)
  distance_to_observed <- function(k) {
    out <- rep(Inf, nrow(k))
    for (i in seq_len(nrow(obs))) out <- pmin(out, rowSums(abs(sweep(k, 2L, obs[i, ], "-"))))
    out
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  walk <- function(steps, n) {
    position <- obs[sample.int(nrow(obs), n, replace = TRUE), , drop = FALSE]
    path <- vector("list", steps)
    for (s in seq_len(steps)) {
      chr <- sample.int(n_chr, n, replace = TRUE)
      move <- sample(c(-1L, 1L), n, replace = TRUE)
      index <- cbind(seq_len(n), chr)
      candidate <- position[index] + move
      blocked <- candidate < min_cn | candidate > max_cn
      move[blocked] <- -move[blocked]
      position[index] <- position[index] + move
      path[[s]] <- position
    }
    list(end = position, path = path)
  }
  keep <- list(`3` = character(), `4` = character())
  paths <- character()
  for (steps in 3:4) {
    batch <- 0L
    while (length(keep[[as.character(steps)]]) < cap && batch < max_batches) {
      batch <- batch + 1L
      w <- walk(steps, 4000L)
      d <- distance_to_observed(w$end)
      hit <- which(d == steps)
      labels <- stt_format_karyotypes(w$end[hit, , drop = FALSE])
      new <- !duplicated(labels) & !labels %in% keep[[as.character(steps)]]
      hit <- hit[new]
      labels <- labels[new]
      room <- cap - length(keep[[as.character(steps)]])
      if (length(labels) > room) {
        hit <- hit[seq_len(room)]
        labels <- labels[seq_len(room)]
      }
      keep[[as.character(steps)]] <- c(keep[[as.character(steps)]], labels)
      for (s in seq_len(steps - 1L)) {
        paths <- c(paths, stt_format_karyotypes(w$path[[s]][hit, , drop = FALSE]))
      }
    }
  }
  if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  scored <- c(keep[["3"]], keep[["4"]])
  paths <- setdiff(unique(paths), c(observed, scored))
  path_distance <- if (length(paths)) distance_to_observed(stt_parse_karyotypes(paths)) else numeric()
  nodes <- data.frame(
    karyotype = c(observed, scored, paths),
    distance = c(rep(0L, length(observed)), rep(3L, length(keep[["3"]])), rep(4L, length(keep[["4"]])),
                 as.integer(path_distance)),
    scored = c(rep(TRUE, length(observed) + length(scored)), rep(FALSE, length(paths))),
    stringsAsFactors = FALSE
  )
  nodes$count_t0 <- 0
  nodes$count_t1 <- 0
  nodes$count_t0[seq_along(observed)] <- counts[observed, 1L]
  nodes$count_t1[seq_along(observed)] <- counts[observed, 2L]
  nodes
}

# Observed karyotypes plus their one- and two-step neighbourhoods. A shell
# larger than its cap is subsampled uniformly; the draw depends only on `seed`.
# The two-step shell is built from the complete one-step shell before either
# is subsampled.
stt_scoring_nodes <- function(counts, min_cn = 1L, max_cn = 8L, d1_cap = 2000L, d2_cap = 2000L, seed = 1L) {
  counts <- as.matrix(counts)
  observed <- rownames(counts)[rowSums(counts) > 0]
  d0 <- stt_parse_karyotypes(observed)
  shell1 <- stt_one_step_neighbours(d0, min_cn, max_cn)
  labels1 <- setdiff(stt_format_karyotypes(shell1), observed)
  d1 <- stt_parse_karyotypes(labels1)
  shell2 <- stt_one_step_neighbours(d1, min_cn, max_cn)
  labels2 <- setdiff(stt_format_karyotypes(shell2), c(observed, labels1))
  shell_sizes <- c(d1 = length(labels1), d2 = length(labels2))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  subsample <- function(labels, cap) {
    if (length(labels) <= cap) return(labels)
    labels <- sort(labels)
    labels[sort(sample.int(length(labels), cap))]
  }
  labels1 <- subsample(labels1, d1_cap)
  labels2 <- subsample(labels2, d2_cap)
  if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  nodes <- data.frame(
    karyotype = c(observed, labels1, labels2),
    distance = rep(0:2, c(length(observed), length(labels1), length(labels2))),
    stringsAsFactors = FALSE
  )
  nodes$count_t0 <- 0
  nodes$count_t1 <- 0
  nodes$count_t0[seq_along(observed)] <- counts[observed, 1L]
  nodes$count_t1[seq_along(observed)] <- counts[observed, 2L]
  attr(nodes, "shell_sizes") <- shell_sizes
  nodes
}

# Reference for relative fitness: count-weighted mean over the well-observed
# karyotypes. Falls back to all observed karyotypes when none reaches
# `min_count`.
stt_reference_weights <- function(nodes, min_count = 5) {
  total <- nodes$count_t0 + nodes$count_t1
  eligible <- nodes$distance == 0L & total >= min_count
  if (!any(eligible)) eligible <- nodes$distance == 0L
  weights <- ifelse(eligible, total, 0)
  weights / sum(weights)
}

stt_auc <- function(score, positive) {
  ok <- is.finite(score) & !is.na(positive)
  score <- score[ok]
  positive <- positive[ok]
  n_pos <- sum(positive)
  n_neg <- sum(!positive)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  (sum(rank(score)[positive]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

stt_interval_score <- function(lower, upper, truth, level = 0.95) {
  alpha <- 1 - level
  (upper - lower) + (2 / alpha) * pmax(lower - truth, 0) + (2 / alpha) * pmax(truth - upper, 0)
}

stt_safe_cor <- function(x, y, method) {
  if (length(x) < 3L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  stats::cor(x, y, method = method)
}

# `prediction` has columns karyotype, fitness and optionally fitness_sd. Nodes
# without a finite prediction are scored as relative fitness zero (the
# no-information answer) and counted in `prediction_fraction`.
stt_score_predictions <- function(nodes, truth_rate, prediction, reference_min_count = 5) {
  stopifnot(length(truth_rate) == nrow(nodes))
  weights <- stt_reference_weights(nodes, reference_min_count)
  match_index <- match(nodes$karyotype, as.character(prediction$karyotype))
  fitness <- as.numeric(prediction$fitness)[match_index]
  fitness_sd <- if ("fitness_sd" %in% names(prediction)) as.numeric(prediction$fitness_sd)[match_index] else rep(NA_real_, nrow(nodes))
  available <- is.finite(fitness)
  reference_ok <- all(available[weights > 0])
  truth_relative <- truth_rate - sum(weights * truth_rate)
  if (reference_ok) {
    pred_relative <- fitness - sum((weights * fitness)[weights > 0])
    pred_relative[!available] <- 0
  } else {
    pred_relative <- rep(0, nrow(nodes))
    available[] <- FALSE
  }
  has_sd <- available & is.finite(fitness_sd) & fitness_sd > 0
  lower <- pred_relative - stats::qnorm(0.975) * fitness_sd
  upper <- pred_relative + stats::qnorm(0.975) * fitness_sd

  node_table <- data.frame(
    karyotype = nodes$karyotype, distance = nodes$distance,
    count_t0 = nodes$count_t0, count_t1 = nodes$count_t1,
    truth_relative = truth_relative, pred_relative = pred_relative,
    pred_sd = ifelse(has_sd, fitness_sd, NA_real_), available = available,
    stringsAsFactors = FALSE
  )

  # One-step moves into each shell: parent one shell closer, child in the shell.
  karyotypes <- stt_parse_karyotypes(nodes$karyotype)
  index_of <- stats::setNames(seq_len(nrow(nodes)), nodes$karyotype)
  edge_rows <- lapply(seq_len(ncol(karyotypes)), function(chr) {
    moved <- karyotypes
    moved[, chr] <- moved[, chr] + 1
    child <- index_of[stt_format_karyotypes(moved)]
    keep <- which(!is.na(child))
    cbind(from = keep, to = unname(child[keep]))
  })
  edges <- do.call(rbind, edge_rows)
  edge_distance <- pmax(nodes$distance[edges[, "from"]], nodes$distance[edges[, "to"]])
  edge_truth <- truth_relative[edges[, "to"]] - truth_relative[edges[, "from"]]
  edge_pred <- pred_relative[edges[, "to"]] - pred_relative[edges[, "from"]]

  # Primary strata are the three distances. Observed karyotypes are also split
  # by support, because most of them are seen in only a few cells.
  total <- nodes$count_t0 + nodes$count_t1
  no_edges <- rep(FALSE, length(edge_distance))
  # One stratum per distance present (0 to 2 for the scoreboard, also 3 and 4
  # for the extrapolation check), then observed karyotypes by support.
  strata <- c(
    lapply(sort(unique(nodes$distance)), function(d) {
      list(name = paste0("d", d), use = nodes$distance == d, edges = edge_distance == d)
    }),
    list(
      list(name = "d0_ge5", use = nodes$distance == 0L & total >= 5, edges = no_edges),
      list(name = "d0_lt5", use = nodes$distance == 0L & total < 5, edges = no_edges)
    )
  )
  by_stratum <- lapply(strata, function(stratum) {
    use <- stratum$use
    if (!any(use)) return(NULL)
    truth <- truth_relative[use]
    pred <- pred_relative[use]
    sd_ok <- has_sd[use]
    edge_use <- stratum$edges
    data.frame(
      stratum = stratum$name,
      n_nodes = sum(use),
      prediction_fraction = mean(available[use]),
      rmse = sqrt(mean((pred - truth)^2)),
      rmse_zero = sqrt(mean(truth^2)),
      # Median absolute error is insensitive to a few extreme predictions.
      median_abs_error = stats::median(abs(pred - truth)),
      median_abs_error_zero = stats::median(abs(truth)),
      extreme_fraction = mean(abs(pred) > 1),
      # Decomposition: rmse^2 = bias^2 + error_sd^2 (population SD). `truth_mean`
      # is how much less fit the shell is than the reference on average.
      bias = mean(pred - truth),
      error_sd = sqrt(mean((pred - truth - mean(pred - truth))^2)),
      truth_mean = mean(truth),
      truth_sd = sqrt(mean((truth - mean(truth))^2)),
      skill = 1 - mean((pred - truth)^2) / mean(truth^2),
      slope = if (isTRUE(stats::var(truth) > 0)) stats::cov(pred, truth) / stats::var(truth) else NA_real_,
      pearson = stt_safe_cor(pred, truth, "pearson"),
      spearman = stt_safe_cor(pred, truth, "spearman"),
      auc_beneficial = stt_auc(pred, truth > 0),
      n_edges = sum(edge_use),
      edge_rmse = if (any(edge_use)) sqrt(mean((edge_pred[edge_use] - edge_truth[edge_use])^2)) else NA_real_,
      edge_rmse_zero = if (any(edge_use)) sqrt(mean(edge_truth[edge_use]^2)) else NA_real_,
      edge_sign_accuracy = if (any(edge_use)) {
        informative <- edge_use & edge_truth != 0
        mean(sign(edge_pred[informative]) == sign(edge_truth[informative]))
      } else NA_real_,
      interval_fraction = mean(sd_ok),
      coverage95 = if (any(sd_ok)) mean(truth[sd_ok] >= lower[use][sd_ok] & truth[sd_ok] <= upper[use][sd_ok]) else NA_real_,
      interval_score95 = if (any(sd_ok)) mean(stt_interval_score(lower[use][sd_ok], upper[use][sd_ok], truth[sd_ok])) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  list(metrics = do.call(rbind, by_stratum), nodes = node_table, reference_available = reference_ok)
}
