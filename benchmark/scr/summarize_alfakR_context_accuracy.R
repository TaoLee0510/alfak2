#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: summarize_alfakR_context_accuracy.R <benchmark_result_dir>", call. = FALSE)
}

base_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- file.path(base_dir, "tables")
repeat_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
repeat_dirs <- repeat_dirs[grepl("repeat_[0-9]+$", basename(repeat_dirs))]
if (!length(repeat_dirs)) {
  stop("No repeat_XX directories found under: ", base_dir, call. = FALSE)
}

repeat_dirs <- repeat_dirs[order(basename(repeat_dirs))]
node_files <- file.path(repeat_dirs, "tables", "node_accuracy.tsv")
if (!all(file.exists(node_files))) {
  stop("Missing node_accuracy.tsv in: ",
       paste(repeat_dirs[!file.exists(node_files)], collapse = ", "),
       call. = FALSE)
}

read_one <- function(path, repeat_id) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  x$fit_repeat <- repeat_id
  x
}

nodes <- do.call(rbind, Map(read_one, node_files, seq_along(node_files)))
nodes <- subset(nodes, status == "ok" & engine %in% c("alfak2", "alfakR"))

cond_keys <- c(
  "fit_repeat", "simulation_id", "lambda_label", "time_start",
  "time_gap", "time_delta", "minobs", "sim_pm", "fit_beta_label"
)
node_keys <- c(cond_keys, "k")

recode_scope <- function(x) {
  ifelse(x == "direct", "fq", ifelse(x == "nn", "NN", x))
}

metric_row <- function(df, error_col = "estimation_error",
                       estimate_col = "estimated_fitness",
                       truth_col = "true_fitness") {
  e <- df[[error_col]]
  est <- df[[estimate_col]]
  truth <- df[[truth_col]]
  centered_error <- (est - mean(est, na.rm = TRUE)) -
    (truth - mean(truth, na.rm = TRUE))
  data.frame(
    n_nodes = length(e),
    mae = mean(abs(e), na.rm = TRUE),
    rmse = sqrt(mean(e^2, na.rm = TRUE)),
    centered_mae = mean(abs(centered_error), na.rm = TRUE),
    centered_rmse = sqrt(mean(centered_error^2, na.rm = TRUE)),
    signed_bias = mean(e, na.rm = TRUE),
    pearson = suppressWarnings(
      if (length(e) >= 3 && sd(est, na.rm = TRUE) > 0 &&
          sd(truth, na.rm = TRUE) > 0) {
        cor(est, truth, use = "complete.obs", method = "pearson")
      } else {
        NA_real_
      }
    ),
    spearman = suppressWarnings(
      if (length(e) >= 3 && sd(est, na.rm = TRUE) > 0 &&
          sd(truth, na.rm = TRUE) > 0) {
        cor(est, truth, use = "complete.obs", method = "spearman")
      } else {
        NA_real_
      }
    )
  )
}

summarise_split <- function(df, group_vars, fun = metric_row) {
  groups <- interaction(df[, group_vars], drop = TRUE, lex.order = TRUE)
  out <- do.call(rbind, lapply(split(df, groups), function(z) {
    cbind(z[1, group_vars, drop = FALSE], fun(z))
  }))
  rownames(out) <- NULL
  out
}

qfun <- function(x) {
  c(
    median = median(x, na.rm = TRUE),
    q25 = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
    q75 = as.numeric(quantile(x, 0.75, na.rm = TRUE))
  )
}

flatten_aggregate <- function(x, id_cols) {
  y <- data.frame(x[, seq_len(id_cols)], check.names = FALSE)
  for (nm in names(x)[-seq_len(id_cols)]) {
    mat <- if (is.matrix(x[[nm]])) x[[nm]] else do.call(rbind, x[[nm]])
    colnames(mat) <- paste0(nm, "_", colnames(mat))
    y <- cbind(y, mat)
  }
  y
}

scope_map <- unique(subset(
  nodes,
  engine == "alfakR",
  select = c(node_keys, "support_scope")
))
names(scope_map)[names(scope_map) == "support_scope"] <- "alfakR_scope"
scope_map$alfakR_scope <- recode_scope(scope_map$alfakR_scope)

scoped_nodes <- merge(nodes, scope_map, by = node_keys)
scoped_nodes$method_label <- ifelse(
  scoped_nodes$engine == "alfakR",
  "alfakR_empirical",
  paste0("alfak2_", scoped_nodes$input_policy)
)
scoped_nodes <- subset(
  scoped_nodes,
  method_label %in% c("alfak2_full", "alfak2_minobs_matched",
                      "alfakR_empirical") &
    alfakR_scope %in% c("fq", "NN", "other")
)

condition_own <- summarise_split(
  scoped_nodes,
  c(cond_keys, "alfakR_scope", "method_label")
)
aggregate_own <- aggregate(
  cbind(n_nodes, mae, rmse, centered_mae, centered_rmse, signed_bias,
        pearson, spearman) ~ minobs + alfakR_scope + method_label,
  condition_own,
  qfun
)
aggregate_own <- flatten_aggregate(aggregate_own, 3)
aggregate_own <- aggregate_own[order(
  match(aggregate_own$alfakR_scope, c("fq", "NN", "other")),
  aggregate_own$minobs,
  match(aggregate_own$method_label,
        c("alfak2_full", "alfak2_minobs_matched", "alfakR_empirical"))
), ]

write.table(
  aggregate_own,
  file.path(out_dir, "alfakR_context_accuracy_each_method.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

own_wide <- reshape(
  condition_own[, c(cond_keys, "alfakR_scope", "method_label", "n_nodes",
                    "mae", "rmse", "centered_mae", "centered_rmse")],
  idvar = c(cond_keys, "alfakR_scope"),
  timevar = "method_label",
  direction = "wide"
)

make_own_delta <- function(policy_label) {
  data.frame(
    own_wide[, c(cond_keys, "alfakR_scope")],
    alfak2_policy = sub("alfak2_", "", policy_label),
    n_nodes_alfak2 = own_wide[[paste0("n_nodes.", policy_label)]],
    n_nodes_alfakR = own_wide[["n_nodes.alfakR_empirical"]],
    alfak2_mae = own_wide[[paste0("mae.", policy_label)]],
    alfakR_mae = own_wide[["mae.alfakR_empirical"]],
    delta_mae = own_wide[[paste0("mae.", policy_label)]] -
      own_wide[["mae.alfakR_empirical"]],
    alfak2_rmse = own_wide[[paste0("rmse.", policy_label)]],
    alfakR_rmse = own_wide[["rmse.alfakR_empirical"]],
    delta_rmse = own_wide[[paste0("rmse.", policy_label)]] -
      own_wide[["rmse.alfakR_empirical"]],
    alfak2_centered_rmse = own_wide[[paste0("centered_rmse.", policy_label)]],
    alfakR_centered_rmse = own_wide[["centered_rmse.alfakR_empirical"]],
    delta_centered_rmse = own_wide[[paste0("centered_rmse.", policy_label)]] -
      own_wide[["centered_rmse.alfakR_empirical"]]
  )
}

own_delta <- rbind(
  make_own_delta("alfak2_full"),
  make_own_delta("alfak2_minobs_matched")
)
aggregate_own_delta <- aggregate(
  cbind(n_nodes_alfak2, n_nodes_alfakR, alfak2_mae, alfakR_mae,
        delta_mae, alfak2_rmse, alfakR_rmse, delta_rmse,
        alfak2_centered_rmse, alfakR_centered_rmse,
        delta_centered_rmse) ~ minobs + alfakR_scope + alfak2_policy,
  own_delta,
  qfun
)
aggregate_own_delta <- flatten_aggregate(aggregate_own_delta, 3)
aggregate_own_delta <- aggregate_own_delta[order(
  match(aggregate_own_delta$alfakR_scope, c("fq", "NN", "other")),
  aggregate_own_delta$minobs,
  match(aggregate_own_delta$alfak2_policy, c("full", "minobs_matched"))
), ]

write.table(
  aggregate_own_delta,
  file.path(out_dir, "alfakR_context_accuracy_each_method_delta_vs_alfakR.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

strict_nodes <- scoped_nodes[, c(
  node_keys, "alfakR_scope", "method_label", "estimated_fitness",
  "true_fitness", "estimation_error"
)]
alfakR_nodes <- subset(strict_nodes, method_label == "alfakR_empirical")
names(alfakR_nodes)[names(alfakR_nodes) == "estimated_fitness"] <-
  "alfakR_estimated"
names(alfakR_nodes)[names(alfakR_nodes) == "estimation_error"] <-
  "alfakR_error"

make_strict <- function(policy_label) {
  alfak2_nodes <- subset(strict_nodes, method_label == policy_label)
  names(alfak2_nodes)[names(alfak2_nodes) == "estimated_fitness"] <-
    "alfak2_estimated"
  names(alfak2_nodes)[names(alfak2_nodes) == "estimation_error"] <-
    "alfak2_error"
  matched <- merge(
    alfak2_nodes[, c(node_keys, "alfakR_scope", "alfak2_estimated",
                     "true_fitness", "alfak2_error")],
    alfakR_nodes[, c(node_keys, "alfakR_scope", "alfakR_estimated",
                     "alfakR_error")],
    by = c(node_keys, "alfakR_scope")
  )
  matched$alfak2_policy <- sub("alfak2_", "", policy_label)
  matched
}

strict <- rbind(
  make_strict("alfak2_full"),
  make_strict("alfak2_minobs_matched")
)

strict_metric <- function(df) {
  e2 <- df$alfak2_error
  er <- df$alfakR_error
  truth <- df$true_fitness
  ce2 <- (df$alfak2_estimated - mean(df$alfak2_estimated, na.rm = TRUE)) -
    (truth - mean(truth, na.rm = TRUE))
  cer <- (df$alfakR_estimated - mean(df$alfakR_estimated, na.rm = TRUE)) -
    (truth - mean(truth, na.rm = TRUE))
  data.frame(
    n_common = length(e2),
    alfak2_mae = mean(abs(e2), na.rm = TRUE),
    alfakR_mae = mean(abs(er), na.rm = TRUE),
    delta_mae = mean(abs(e2), na.rm = TRUE) - mean(abs(er), na.rm = TRUE),
    alfak2_rmse = sqrt(mean(e2^2, na.rm = TRUE)),
    alfakR_rmse = sqrt(mean(er^2, na.rm = TRUE)),
    delta_rmse = sqrt(mean(e2^2, na.rm = TRUE)) -
      sqrt(mean(er^2, na.rm = TRUE)),
    alfak2_centered_rmse = sqrt(mean(ce2^2, na.rm = TRUE)),
    alfakR_centered_rmse = sqrt(mean(cer^2, na.rm = TRUE)),
    delta_centered_rmse = sqrt(mean(ce2^2, na.rm = TRUE)) -
      sqrt(mean(cer^2, na.rm = TRUE)),
    alfak2_abs_better_rate = mean(abs(e2) < abs(er), na.rm = TRUE),
    alfakR_abs_better_rate = mean(abs(er) < abs(e2), na.rm = TRUE)
  )
}

strict_condition <- do.call(rbind, lapply(
  split(
    strict,
    interaction(strict[, c(cond_keys, "alfakR_scope", "alfak2_policy")],
                drop = TRUE, lex.order = TRUE)
  ),
  function(z) {
    cbind(z[1, c(cond_keys, "alfakR_scope", "alfak2_policy"),
            drop = FALSE], strict_metric(z))
  }
))
rownames(strict_condition) <- NULL

aggregate_strict <- aggregate(
  cbind(n_common, alfak2_mae, alfakR_mae, delta_mae, alfak2_rmse,
        alfakR_rmse, delta_rmse, alfak2_centered_rmse,
        alfakR_centered_rmse, delta_centered_rmse,
        alfak2_abs_better_rate, alfakR_abs_better_rate) ~
    minobs + alfakR_scope + alfak2_policy,
  strict_condition,
  qfun
)
aggregate_strict <- flatten_aggregate(aggregate_strict, 3)
aggregate_strict <- aggregate_strict[order(
  match(aggregate_strict$alfakR_scope, c("fq", "NN", "other")),
  aggregate_strict$minobs,
  match(aggregate_strict$alfak2_policy, c("full", "minobs_matched"))
), ]

write.table(
  aggregate_strict,
  file.path(out_dir, "alfakR_context_accuracy_strict_common_karyotypes_delta.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message("Wrote alfakR-context accuracy summaries to: ", out_dir)
