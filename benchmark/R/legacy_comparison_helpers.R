read_landscape_csv_or_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("rds", "rda")) {
    readRDS(path)
  } else {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

normalize_landscape_for_comparison <- function(x, patient_id, method) {
  if (is.null(x) || !nrow(x)) return(data.frame())
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!"k" %in% names(x) && "karyotype" %in% names(x)) x$k <- x$karyotype
  if (!"mean" %in% names(x) && "fitness_mean" %in% names(x)) x$mean <- x$fitness_mean
  if (!"median" %in% names(x)) x$median <- x$mean
  if (!"sd" %in% names(x) && "fitness_sd" %in% names(x)) x$sd <- x$fitness_sd
  if (!"sd" %in% names(x)) x$sd <- NA_real_
  if (!"fq" %in% names(x)) x$fq <- identical(method, "alfak2") && x$support_tier %in% "directly_informed"
  if (!"nn" %in% names(x)) {
    x$nn <- identical(method, "alfak2") &&
      x$support_tier %in% c("local_borrowed", "weakly_supported", "graph_borrowed")
  }
  if (!"support_tier" %in% names(x)) {
    x$support_tier <- ifelse(x$fq %in% TRUE, "directly_informed",
      ifelse(x$nn %in% TRUE, "legacy_nn", "legacy_other")
    )
  }
  data.frame(
    patient_id = patient_id,
    method = method,
    k = as.character(x$k),
    mean = as.numeric(x$mean),
    median = as.numeric(x$median),
    sd = as.numeric(x$sd),
    fq = x$fq %in% TRUE,
    nn = x$nn %in% TRUE,
    support_tier = as.character(x$support_tier),
    stringsAsFactors = FALSE
  )
}

read_old_alfakR_patient <- function(old_fit_dir, patient_id) {
  landscape <- read_landscape_csv_or_rds(file.path(old_fit_dir, patient_id, "landscape.Rds"))
  xval <- read_landscape_csv_or_rds(file.path(old_fit_dir, patient_id, "xval.Rds"))
  diag <- read_landscape_csv_or_rds(file.path(old_fit_dir, patient_id, "nn_prior_diagnostics.Rds"))
  list(
    patient_id = patient_id,
    landscape = normalize_landscape_for_comparison(landscape, patient_id, "alfakR_nn_prior_empirical_two_shell"),
    xval = if (length(xval) == 1L) as.numeric(xval) else NA_real_,
    diagnostics = diag
  )
}

read_new_alfak2_patient <- function(new_fit_dir, patient_id) {
  landscape <- read_landscape_csv_or_rds(file.path(new_fit_dir, patient_id, "landscape.rds"))
  bundle <- read_landscape_csv_or_rds(file.path(new_fit_dir, patient_id, "bundle.rds"))
  xval <- read_landscape_csv_or_rds(file.path(new_fit_dir, patient_id, "xval.Rds"))
  diag <- if (!is.null(bundle$fit$local$diagnostics)) bundle$fit$local$diagnostics else list()
  list(
    patient_id = patient_id,
    landscape = normalize_landscape_for_comparison(landscape, patient_id, "alfak2"),
    xval = if (length(xval) == 1L) as.numeric(xval) else NA_real_,
    diagnostics = diag,
    input = bundle$input
  )
}

common_landscape_metrics <- function(old, new, top_k = c(10L, 25L, 50L)) {
  old_l <- old$landscape
  new_l <- new$landscape
  common <- merge(
    old_l,
    new_l,
    by = c("patient_id", "k"),
    suffixes = c("_old", "_new")
  )
  ok <- is.finite(common$mean_old) & is.finite(common$mean_new)
  top_values <- lapply(top_k, function(k) {
    requested_k <- as.integer(k)
    k <- min(as.integer(k), sum(ok))
    if (k <= 0L) {
      return(NA_real_)
    }
    idx <- which(ok)
    top_old <- idx[order(common$mean_old[idx], decreasing = TRUE)[seq_len(k)]]
    top_new <- idx[order(common$mean_new[idx], decreasing = TRUE)[seq_len(k)]]
    out <- length(intersect(top_old, top_new)) / k
    names(out) <- paste0("top", requested_k, "_overlap")
    out
  })
  top_values <- as.data.frame(as.list(unlist(top_values)), stringsAsFactors = FALSE)
  old_direct <- old_l$fq %in% TRUE
  new_direct <- new_l$fq %in% TRUE
  data.frame(
    patient_id = old$patient_id,
    n_old = nrow(old_l),
    n_new = nrow(new_l),
    n_common = nrow(common),
    old_direct = sum(old_direct),
    new_direct = sum(new_direct),
    common_fraction_of_old = nrow(common) / max(1L, nrow(old_l)),
    common_fraction_of_new = nrow(common) / max(1L, nrow(new_l)),
    old_xval_r2 = old$xval,
    new_xval_r2 = new$xval,
    new_convergence = if (!is.null(new$diagnostics$convergence)) new$diagnostics$convergence else NA_integer_,
    new_gradient_norm = if (!is.null(new$diagnostics$gradient_norm)) new$diagnostics$gradient_norm else NA_real_,
    spearman_common = suppressWarnings(stats::cor(common$mean_old[ok], common$mean_new[ok], method = "spearman")),
    pearson_common = suppressWarnings(stats::cor(scale(common$mean_old[ok])[, 1], scale(common$mean_new[ok])[, 1], method = "pearson")),
    old_common_mean_range = diff(range(common$mean_old[ok], na.rm = TRUE)),
    new_common_mean_range = diff(range(common$mean_new[ok], na.rm = TRUE)),
    mean_old_sd = mean(old_l$sd, na.rm = TRUE),
    mean_new_sd = mean(new_l$sd, na.rm = TRUE),
    top_values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

patient_ids_from_fit_dirs <- function(new_fit_dir, old_fit_dir, patient_subset = NULL) {
  new_ids <- basename(list.dirs(new_fit_dir, recursive = FALSE, full.names = TRUE))
  old_ids <- basename(list.dirs(old_fit_dir, recursive = FALSE, full.names = TRUE))
  ids <- intersect(new_ids, old_ids)
  ids <- ids[grepl("^P[0-9]+$", ids)]
  if (!is.null(patient_subset) && length(patient_subset)) {
    ids <- intersect(as.character(patient_subset), ids)
  }
  sort_pid_levels_bench(ids)
}

profile_from_landscape <- function(landscape) {
  beneficial_cnv_profile_alfak2(landscape)
}

compare_beneficial_profiles <- function(old, new) {
  old_prof <- profile_from_landscape(old$landscape)
  new_prof <- profile_from_landscape(new$landscape)
  moves <- union(names(old_prof$proportion), names(new_prof$proportion))
  old_prop <- old_prof$proportion[moves]
  new_prop <- new_prof$proportion[moves]
  old_valid <- old_prof$valid_n[moves]
  new_valid <- new_prof$valid_n[moves]
  ok <- is.finite(old_prop) & is.finite(new_prop)
  edge_old <- old_prof$edges
  edge_new <- new_prof$edges
  common_edge <- merge(edge_old, edge_new, by = c("parent_k", "child_k", "chr", "direction", "move"), suffixes = c("_old", "_new"))
  edge_ok <- is.finite(common_edge$delta_old) & is.finite(common_edge$delta_new)
  list(
    profile = data.frame(
      patient_id = old$patient_id,
      move = moves,
      old_proportion = as.numeric(old_prop),
      new_proportion = as.numeric(new_prop),
      delta_proportion = as.numeric(new_prop - old_prop),
      old_valid_n = as.integer(old_valid),
      new_valid_n = as.integer(new_valid),
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      patient_id = old$patient_id,
      n_common_moves = sum(ok),
      beneficial_prop_spearman = suppressWarnings(stats::cor(old_prop[ok], new_prop[ok], method = "spearman")),
      beneficial_prop_pearson = suppressWarnings(stats::cor(old_prop[ok], new_prop[ok], method = "pearson")),
      beneficial_prop_mae = mean(abs(new_prop[ok] - old_prop[ok]), na.rm = TRUE),
      old_valid_edges = sum(old_prof$valid_n, na.rm = TRUE),
      new_valid_edges = sum(new_prof$valid_n, na.rm = TRUE),
      common_edges = sum(edge_ok),
      common_edge_delta_spearman = suppressWarnings(stats::cor(common_edge$delta_old[edge_ok], common_edge$delta_new[edge_ok], method = "spearman")),
      common_edge_beneficial_agreement = mean(common_edge$beneficial_old[edge_ok] == common_edge$beneficial_new[edge_ok], na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    old_profile = old_prof,
    new_profile = new_prof
  )
}

run_alfakR_alfak2_real_comparison <- function(repo_dir,
                                              new_fit_dir = file.path(repo_dir, "benchmark", "results", "real_samples", "fits"),
                                              old_fit_dir = "/Users/4482173/Documents/GitHub/alfakR/benchmark/results/fits/nn_prior_empirical_two_shell/pm_0.00005/MINOBS_20",
                                              patient_subset = NULL,
                                              output_subdir = "legacy_comparison") {
  dirs <- ensure_benchmark_dirs(repo_dir, output_subdir)
  ids <- patient_ids_from_fit_dirs(new_fit_dir, old_fit_dir, patient_subset)
  if (!length(ids)) stop("No shared patient ids found between old and new fit directories.", call. = FALSE)
  old <- lapply(ids, read_old_alfakR_patient, old_fit_dir = old_fit_dir)
  new <- lapply(ids, read_new_alfak2_patient, new_fit_dir = new_fit_dir)
  names(old) <- names(new) <- ids

  metric_tbl <- do.call(rbind, Map(common_landscape_metrics, old, new))
  profile_list <- Map(compare_beneficial_profiles, old, new)
  beneficial_profile_tbl <- do.call(rbind, lapply(profile_list, `[[`, "profile"))
  beneficial_summary_tbl <- do.call(rbind, lapply(profile_list, `[[`, "summary"))
  landscape_long <- do.call(rbind, c(lapply(old, `[[`, "landscape"), lapply(new, `[[`, "landscape")))

  utils::write.csv(metric_tbl, file.path(dirs$tables, "alfakR_vs_alfak2_landscape_metrics.csv"), row.names = FALSE)
  utils::write.csv(beneficial_profile_tbl, file.path(dirs$tables, "alfakR_vs_alfak2_beneficial_profile_by_move.csv"), row.names = FALSE)
  utils::write.csv(beneficial_summary_tbl, file.path(dirs$tables, "alfakR_vs_alfak2_beneficial_summary.csv"), row.names = FALSE)
  utils::write.csv(landscape_long, file.path(dirs$tables, "alfakR_vs_alfak2_landscape_long.csv"), row.names = FALSE)

  saveRDS(
    list(
      patients = ids,
      landscape_metrics = metric_tbl,
      beneficial_profile = beneficial_profile_tbl,
      beneficial_summary = beneficial_summary_tbl,
      landscape_long = landscape_long,
      old = old,
      new = new,
      dirs = dirs
    ),
    file.path(dirs$results, "alfakR_vs_alfak2_real_comparison.rds")
  )

  list(
    patients = ids,
    landscape_metrics = metric_tbl,
    beneficial_profile = beneficial_profile_tbl,
    beneficial_summary = beneficial_summary_tbl,
    landscape_long = landscape_long,
    dirs = dirs
  )
}

plot_node_count_comparison <- function(metric_tbl) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  df <- rbind(
    data.frame(patient_id = metric_tbl$patient_id, method = "alfakR", n_nodes = metric_tbl$n_old),
    data.frame(patient_id = metric_tbl$patient_id, method = "alfak2", n_nodes = metric_tbl$n_new)
  )
  df$patient_id <- factor(df$patient_id, levels = sort_pid_levels_bench(df$patient_id))
  ggplot2::ggplot(df, ggplot2::aes(.data$patient_id, .data$n_nodes, fill = .data$method)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Patient", y = "Landscape nodes, log10 scale", fill = "Method") +
    ggplot2::theme_bw()
}

plot_common_rank_comparison <- function(metric_tbl) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  metric_tbl$patient_id <- factor(metric_tbl$patient_id, levels = sort_pid_levels_bench(metric_tbl$patient_id))
  ggplot2::ggplot(metric_tbl, ggplot2::aes(.data$patient_id, .data$spearman_common)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70") +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::coord_cartesian(ylim = c(-1, 1)) +
    ggplot2::labs(x = "Patient", y = "Spearman correlation on common karyotypes") +
    ggplot2::theme_bw()
}

plot_beneficial_delta_heatmap <- function(profile_tbl) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  chr <- suppressWarnings(as.integer(sub("[+-]$", "", profile_tbl$move)))
  profile_tbl$move <- factor(profile_tbl$move, levels = one_step_move_levels(max(chr, na.rm = TRUE)))
  profile_tbl$patient_id <- factor(profile_tbl$patient_id, levels = sort_pid_levels_bench(profile_tbl$patient_id))
  ggplot2::ggplot(profile_tbl, ggplot2::aes(.data$move, .data$patient_id, fill = .data$delta_proportion)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#2B6CB0", mid = "#F7F7F7", high = "#C53030", midpoint = 0, na.value = "#BDBDBD") +
    ggplot2::labs(x = "One-step CN event", y = "Patient", fill = "alfak2 - alfakR\nbeneficial proportion") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}
