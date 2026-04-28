resolve_alfak2_repo <- function() {
  candidates <- unique(c(
    tryCatch(dirname(knitr::current_input(dir = TRUE)), error = function(e) NA_character_),
    getwd()
  ))
  for (cand in candidates) {
    if (is.na(cand) || !nzchar(cand)) next
    cand <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    probes <- c(cand, file.path(cand, ".."), file.path(cand, "..", ".."))
    for (probe in probes) {
      probe <- normalizePath(probe, winslash = "/", mustWork = FALSE)
      if (file.exists(file.path(probe, "DESCRIPTION")) &&
          dir.exists(file.path(probe, "benchmark"))) {
        return(probe)
      }
    }
  }
  stop("Could not locate the alfak2 repository root.", call. = FALSE)
}

load_alfak2_for_benchmark <- function(repo_dir) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_dir, quiet = TRUE)
  } else {
    library(alfak2)
  }
  invisible(TRUE)
}

ensure_benchmark_dirs <- function(repo_dir, subdir = "real_samples") {
  out <- list(
    results = file.path(repo_dir, "benchmark", "results", subdir),
    tables = file.path(repo_dir, "benchmark", "results", subdir, "tables"),
    figures = file.path(repo_dir, "benchmark", "results", subdir, "figures"),
    fits = file.path(repo_dir, "benchmark", "results", subdir, "fits")
  )
  for (path in out) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  out
}

sort_pid_levels_bench <- function(x) {
  x <- unique(as.character(x))
  num <- suppressWarnings(as.integer(sub("^P", "", x)))
  x[order(ifelse(is.na(num), Inf, num), x)]
}

available_real_patient_ids <- function(data_dir) {
  dd <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
  sort_pid_levels_bench(dd[grepl("^P[0-9]+$", dd)])
}

read_real_sample_metadata <- function(data_dir) {
  meta_path <- file.path(data_dir, "meta_data.xlsx")
  if (!file.exists(meta_path)) return(NULL)
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package `readxl` is required to read benchmark/data/meta_data.xlsx.", call. = FALSE)
  }
  as.data.frame(readxl::read_excel(meta_path), stringsAsFactors = FALSE)
}

patient_manifest <- function(data_dir,
                             patient_id,
                             stage_order = c("Primary", "Recurrent")) {
  meta <- read_real_sample_metadata(data_dir)
  if (!is.null(meta)) {
    needed <- c("ID", "Stage", "pid", "ID_alt", "Age")
    miss <- setdiff(needed, names(meta))
    if (length(miss)) {
      stop("Missing metadata columns: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    tab <- meta[as.character(meta$pid) == patient_id, , drop = FALSE]
    if (!nrow(tab)) stop("No metadata rows for patient ", patient_id, call. = FALSE)
    effective_id <- ifelse(!is.na(tab$ID_alt) & nzchar(tab$ID_alt), tab$ID_alt, tab$ID)
    age <- suppressWarnings(as.numeric(tab$Age))
    out <- data.frame(
      patient_id = patient_id,
      effective_id = as.character(effective_id),
      stage = as.character(tab$Stage),
      age = age,
      stringsAsFactors = FALSE
    )
  } else {
    rds <- list.files(
      file.path(data_dir, patient_id, "perspectives"),
      pattern = "TranscriptomePerspective[.]rds$",
      recursive = TRUE,
      full.names = TRUE
    )
    if (length(rds) < 2L) stop("Need at least two perspective RDS files for ", patient_id, call. = FALSE)
    effective_id <- basename(dirname(rds))
    out <- data.frame(
      patient_id = patient_id,
      effective_id = effective_id,
      stage = paste0("t", seq_along(effective_id) - 1L),
      age = NA_real_,
      stringsAsFactors = FALSE
    )
    stage_order <- out$stage
  }

  out$rds_path <- file.path(data_dir, patient_id, "perspectives", out$effective_id, "TranscriptomePerspective.rds")
  out$exists <- file.exists(out$rds_path)
  out$stage_rank <- match(out$stage, stage_order)
  out$stage_rank[is.na(out$stage_rank)] <- seq_len(sum(is.na(out$stage_rank))) + length(stage_order)
  out <- out[order(out$stage_rank, out$effective_id), , drop = FALSE]
  out <- out[out$exists, , drop = FALSE]
  if (nrow(out) < 2L) stop("Need two available samples for ", patient_id, call. = FALSE)
  out <- out[seq_len(2L), , drop = FALSE]

  delta_years <- suppressWarnings(abs(diff(out$age)))
  if (!is.finite(delta_years) || delta_years <= 0) delta_years <- 1
  out$dt_years <- c(0, delta_years)
  out
}

read_transcriptome_perspective <- function(rds_path) {
  obj <- readRDS(rds_path)
  if (!is.list(obj) || !"profile" %in% names(obj)) {
    stop("Invalid TranscriptomePerspective object: ", rds_path, call. = FALSE)
  }
  profile <- as.matrix(obj$profile)
  storage.mode(profile) <- "numeric"
  list(
    passaging_id = obj$passaging_id,
    perspective_type = obj$perspective_type,
    profile = round(profile)
  )
}

profile_to_karyotype_strings <- function(profile_mat) {
  x <- t(round(as.matrix(profile_mat)))
  apply(x, 1, paste, collapse = ".")
}

diploid_label <- function(n_chr) paste(rep(2L, n_chr), collapse = ".")

build_real_count_matrix <- function(data_dir,
                                    patient_id,
                                    stage_order = c("Primary", "Recurrent"),
                                    min_total_count = 5,
                                    drop_diploid = TRUE) {
  mf <- patient_manifest(data_dir, patient_id, stage_order = stage_order)
  long <- do.call(rbind, lapply(seq_len(nrow(mf)), function(i) {
    obj <- read_transcriptome_perspective(mf$rds_path[i])
    data.frame(
      karyotype = profile_to_karyotype_strings(obj$profile),
      timepoint = paste0("t", i - 1L),
      stage = mf$stage[i],
      effective_id = mf$effective_id[i],
      stringsAsFactors = FALSE
    )
  }))
  counts <- as.matrix(table(long$karyotype, factor(long$timepoint, levels = c("t0", "t1"))))
  storage.mode(counts) <- "integer"
  counts <- counts[rowSums(counts) >= min_total_count, , drop = FALSE]
  if (drop_diploid && nrow(counts)) {
    dlab <- diploid_label(length(strsplit(rownames(counts)[1], ".", fixed = TRUE)[[1]]))
    counts <- counts[rownames(counts) != dlab, , drop = FALSE]
  }
  if (!nrow(counts)) stop("No karyotypes remain after filtering for ", patient_id, call. = FALSE)
  list(
    patient_id = patient_id,
    manifest = mf,
    counts = counts,
    dt = max(mf$dt_years, na.rm = TRUE),
    long = long
  )
}

landscape_from_alfak2_fit <- function(fit, patient_id) {
  s <- alfak2::summarize_alfak2(fit, layer = "global")
  data.frame(
    patient_id = patient_id,
    k = as.character(s$karyotype),
    mean = as.numeric(s$fitness_mean),
    median = as.numeric(s$fitness_mean),
    sd = as.numeric(s$fitness_sd),
    conf_low = as.numeric(s$conf_low),
    conf_high = as.numeric(s$conf_high),
    support_tier = as.character(s$support_tier),
    fq = as.character(s$support_tier) == "directly_informed",
    nn = as.character(s$support_tier) %in% c("local_borrowed", "weakly_supported"),
    stringsAsFactors = FALSE
  )
}

xval_r2r <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- as.numeric(obs[ok])
  pred <- as.numeric(pred[ok])
  if (length(obs) < 2L) return(NA_real_)
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  denom <- sum((obs - mean(obs))^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  1 - sum((pred - obs)^2) / denom
}

one_step_neighbor_labels_bench <- function(label, min_cn = 0) {
  k <- as.integer(strsplit(as.character(label), ".", fixed = TRUE)[[1]])
  if (anyNA(k)) return(character(0))
  out <- character(0)
  for (chr in seq_along(k)) {
    for (direction in c(1L, -1L)) {
      child <- k
      child[chr] <- child[chr] + direction
      if (any(child < min_cn)) next
      out <- c(out, paste(child, collapse = "."))
    }
  }
  unique(out)
}

alfak2_heldout_xval <- function(fit,
                                 min_cn = 0,
                                 max_folds = Inf,
                                 seed = 1) {
  if (!inherits(fit, "alfak2_fit")) {
    stop("`fit` must be an alfak2_fit object.", call. = FALSE)
  }
  local_summary <- fit$local$summary
  global_graph <- fit$global$graph
  if (is.null(local_summary) || !nrow(local_summary)) {
    return(list(R2R = NA_real_, detail = data.frame(), status = "no_local_summary"))
  }
  direct_labels <- as.character(local_summary$karyotype[local_summary$support_tier == "directly_informed"])
  anchor_labels <- as.character(local_summary$karyotype[is.finite(local_summary$fitness_mean)])
  if (!length(direct_labels) || length(anchor_labels) < 3L) {
    return(list(R2R = NA_real_, detail = data.frame(), status = "too_few_anchors"))
  }

  ids <- unlist(lapply(seq_along(direct_labels), function(i) {
    ki <- unique(c(direct_labels[i], one_step_neighbor_labels_bench(direct_labels[i], min_cn = min_cn)))
    ki <- ki[ki %in% anchor_labels]
    out <- rep(i, length(ki))
    names(out) <- ki
    out
  }))
  ids <- ids[!duplicated(names(ids))]
  fold_ids <- unique(ids)
  if (is.finite(max_folds) && max_folds > 0 && length(fold_ids) > max_folds) {
    set.seed(as.integer(seed))
    fold_ids <- sort(sample(fold_ids, size = as.integer(max_folds)))
  }
  hp <- fit$global$hyperparameters
  detail <- vector("list", length(fold_ids))

  for (i in seq_along(fold_ids)) {
    fold <- fold_ids[i]
    test_labels <- names(ids)[ids == fold]
    heldout_local <- fit$local
    heldout_local$summary <- local_summary[!(as.character(local_summary$karyotype) %in% test_labels), , drop = FALSE]
    if (nrow(heldout_local$summary) < 2L) {
      detail[[i]] <- data.frame(
        fold = fold,
        k = test_labels,
        state_class = ifelse(test_labels %in% direct_labels, "fq", "nn"),
        validation = local_summary$fitness_mean[match(test_labels, local_summary$karyotype)],
        estimate = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    pred_fit <- try(
      alfak2::fit_graph_posterior(
        heldout_local,
        global_graph,
        lambda_l_grid = hp$lambda_l,
        lambda_e_grid = hp$lambda_e,
        sigma_obs_grid = hp$sigma_obs
      ),
      silent = TRUE
    )
    validation <- local_summary$fitness_mean[match(test_labels, local_summary$karyotype)]
    estimate <- rep(NA_real_, length(test_labels))
    if (!inherits(pred_fit, "try-error")) {
      idx <- match(test_labels, as.character(pred_fit$summary$karyotype))
      ok <- !is.na(idx)
      estimate[ok] <- pred_fit$summary$fitness_mean[idx[ok]]
    }
    detail[[i]] <- data.frame(
      fold = fold,
      k = test_labels,
      state_class = ifelse(test_labels %in% direct_labels, "fq", "nn"),
      validation = as.numeric(validation),
      estimate = as.numeric(estimate),
      stringsAsFactors = FALSE
    )
  }

  detail <- do.call(rbind, detail)
  detail <- detail[is.finite(detail$validation) & is.finite(detail$estimate), , drop = FALSE]
  list(
    R2R = xval_r2r(detail$validation, detail$estimate),
    detail = detail,
    status = if (nrow(detail) >= 2L) "ok" else "too_few_complete_predictions",
    max_folds = max_folds,
    evaluated_folds = length(fold_ids)
  )
}

attach_alfak2_xval <- function(fit, xval) {
  if (is.null(xval)) return(fit)
  fit$xval <- xval
  if (is.null(fit$diagnostics)) fit$diagnostics <- list()
  fit$diagnostics$xval <- list(
    R2R = xval$R2R,
    status = xval$status,
    n = if (!is.null(xval$detail)) nrow(xval$detail) else NA_integer_,
    evaluated_folds = xval$evaluated_folds,
    max_folds = xval$max_folds
  )
  fit
}

fit_real_patient_alfak2 <- function(repo_dir,
                                    patient_id,
                                    min_total_count = 20,
                                    drop_diploid = TRUE,
                                    beta = 0.00005,
                                    min_cn = 0,
                                    max_cn = 5,
                                    local_shell_depth = 0,
                                    global_extra_shell = 1,
                                    max_nodes = 4000,
                                    compute_xval = TRUE,
                                    xval_max_folds = Inf,
                                    xval_seed = 1,
                                    force = FALSE,
                                    output_subdir = "real_samples",
                                    control = list(eval.max = 250, iter.max = 250)) {
  dirs <- ensure_benchmark_dirs(repo_dir, output_subdir)
  fit_dir <- file.path(dirs$fits, patient_id)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  bundle_path <- file.path(fit_dir, "bundle.rds")
  fit_path <- file.path(fit_dir, "alfak2_fit.rds")
  landscape_path <- file.path(fit_dir, "landscape.rds")
  counts_path <- file.path(fit_dir, "counts.csv")
  xval_path <- file.path(fit_dir, "xval.Rds")
  xval_detail_path <- file.path(fit_dir, "xval_detail.csv")

  if (!force && file.exists(fit_path) && file.exists(landscape_path) && file.exists(bundle_path)) {
    bundle <- readRDS(bundle_path)
    bundle$min_total_count <- min_total_count
    bundle$output_subdir <- output_subdir
    if (isTRUE(compute_xval) && !file.exists(xval_path)) {
      xval <- alfak2_heldout_xval(bundle$fit, min_cn = min_cn, max_folds = xval_max_folds, seed = xval_seed)
      saveRDS(xval$R2R, xval_path)
      utils::write.csv(xval$detail, xval_detail_path, row.names = FALSE)
      bundle$fit <- attach_alfak2_xval(bundle$fit, xval)
      bundle$xval <- xval
      bundle$paths$xval <- xval_path
      bundle$paths$xval_detail <- xval_detail_path
      saveRDS(bundle$fit, fit_path)
      saveRDS(bundle, bundle_path)
    } else if (file.exists(xval_path)) {
      bundle$xval <- list(
        R2R = readRDS(xval_path),
        detail = if (file.exists(xval_detail_path)) utils::read.csv(xval_detail_path, stringsAsFactors = FALSE) else data.frame(),
        status = "loaded"
      )
      bundle$fit <- attach_alfak2_xval(bundle$fit, bundle$xval)
      bundle$paths$xval <- xval_path
      bundle$paths$xval_detail <- xval_detail_path
      saveRDS(bundle$fit, fit_path)
      saveRDS(bundle, bundle_path)
    }
    return(bundle)
  }

  built <- build_real_count_matrix(
    file.path(repo_dir, "benchmark", "data"),
    patient_id = patient_id,
    min_total_count = min_total_count,
    drop_diploid = drop_diploid
  )
  utils::write.csv(
    data.frame(karyotype = rownames(built$counts), built$counts, check.names = FALSE),
    counts_path,
    row.names = FALSE
  )
  fit <- alfak2::fit_alfak2(
    built$counts,
    dt = built$dt,
    beta = beta,
    min_cn = min_cn,
    max_cn = max_cn,
    local_shell_depth = local_shell_depth,
    global_extra_shell = global_extra_shell,
    max_nodes = max_nodes,
    control = control
  )
  landscape <- landscape_from_alfak2_fit(fit, patient_id)
  xval <- if (isTRUE(compute_xval)) {
    alfak2_heldout_xval(fit, min_cn = min_cn, max_folds = xval_max_folds, seed = xval_seed)
  } else {
    NULL
  }
  utils::write.csv(landscape, file.path(fit_dir, "landscape.csv"), row.names = FALSE)
  if (!is.null(xval)) {
    saveRDS(xval$R2R, xval_path)
    utils::write.csv(xval$detail, xval_detail_path, row.names = FALSE)
    fit <- attach_alfak2_xval(fit, xval)
  }
  saveRDS(fit, fit_path)
  saveRDS(landscape, landscape_path)
  bundle <- list(
    patient_id = patient_id,
    min_total_count = min_total_count,
    output_subdir = output_subdir,
    input = built,
    fit = fit,
    landscape = landscape,
    xval = xval,
    paths = list(
      fit = fit_path,
      landscape = landscape_path,
      counts = counts_path,
      xval = xval_path,
      xval_detail = xval_detail_path
    )
  )
  saveRDS(bundle, bundle_path)
  bundle
}

selected_real_patient_ids <- function(repo_dir, patient_subset = NULL) {
  data_dir <- file.path(repo_dir, "benchmark", "data")
  pids <- available_real_patient_ids(data_dir)
  if (!is.null(patient_subset) && length(patient_subset)) {
    pids <- intersect(sort_pid_levels_bench(patient_subset), pids)
  }
  if (!length(pids)) stop("No patient ids selected.", call. = FALSE)
  pids
}

run_real_sample_landscapes <- function(repo_dir,
                                       patient_subset = NULL,
                                       continue_on_error = TRUE,
                                       ...) {
  pids <- selected_real_patient_ids(repo_dir, patient_subset = patient_subset)
  names(pids) <- pids
  out <- lapply(pids, function(pid) {
    tryCatch(
      fit_real_patient_alfak2(repo_dir, pid, ...),
      error = function(e) {
        if (!continue_on_error) stop(e)
        warning("Skipping ", pid, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  })
  out[!vapply(out, is.null, logical(1))]
}

run_real_sample_landscape_grid <- function(repo_dir,
                                           patient_subset = NULL,
                                           min_total_count_grid = c(5L, 10L, 20L),
                                           parallel_workers = 1L,
                                           output_root = "real_samples",
                                           continue_on_error = TRUE,
                                           ...) {
  pids <- selected_real_patient_ids(repo_dir, patient_subset = patient_subset)
  min_total_count_grid <- as.integer(min_total_count_grid)
  if (!length(min_total_count_grid) || anyNA(min_total_count_grid) || any(min_total_count_grid < 1L)) {
    stop("`min_total_count_grid` must contain positive integers.", call. = FALSE)
  }
  min_total_count_grid <- unique(min_total_count_grid)
  tasks <- expand.grid(
    patient_id = pids,
    min_total_count = min_total_count_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tasks$task_id <- paste0(tasks$patient_id, "_MINOBS_", tasks$min_total_count)

  run_one <- function(i) {
    task <- tasks[i, , drop = FALSE]
    out_subdir <- file.path(output_root, paste0("MINOBS_", task$min_total_count))
    tryCatch(
      {
        bundle <- fit_real_patient_alfak2(
          repo_dir = repo_dir,
          patient_id = task$patient_id,
          min_total_count = task$min_total_count,
          output_subdir = out_subdir,
          ...
        )
        bundle$task_id <- task$task_id
        bundle$min_total_count <- task$min_total_count
        bundle$output_subdir <- out_subdir
        bundle
      },
      error = function(e) {
        if (!continue_on_error) stop(e)
        warning("Skipping ", task$task_id, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  }

  idx <- seq_len(nrow(tasks))
  parallel_workers <- suppressWarnings(as.numeric(parallel_workers)[1])
  if (is.na(parallel_workers) || parallel_workers < 1) parallel_workers <- 1
  if (is.infinite(parallel_workers)) parallel_workers <- length(idx)
  parallel_workers <- as.integer(min(parallel_workers, length(idx)))
  if (parallel_workers > 1L && .Platform$OS.type != "windows") {
    out <- parallel::mclapply(idx, run_one, mc.cores = parallel_workers)
  } else {
    if (parallel_workers > 1L && .Platform$OS.type == "windows") {
      warning("Forked parallelism is unavailable on Windows; running tasks serially.", call. = FALSE)
    }
    out <- lapply(idx, run_one)
  }
  names(out) <- tasks$task_id
  out <- out[!vapply(out, is.null, logical(1))]
  out
}

parse_karyotype_matrix_bench <- function(labels) {
  parts <- strsplit(as.character(labels), ".", fixed = TRUE)
  n_chr <- unique(lengths(parts))
  if (length(n_chr) != 1L) stop("Karyotypes have inconsistent chromosome counts.", call. = FALSE)
  matrix(as.integer(unlist(parts, use.names = FALSE)), ncol = n_chr, byrow = TRUE)
}

one_step_move_levels <- function(n_chr) {
  as.vector(rbind(paste0(seq_len(n_chr), "+"), paste0(seq_len(n_chr), "-")))
}

one_step_edge_table <- function(landscape_df,
                                start_tiers = "directly_informed",
                                min_cn = 0) {
  labels <- as.character(landscape_df$k)
  k_mat <- parse_karyotype_matrix_bench(labels)
  n_chr <- ncol(k_mat)
  fitness <- as.numeric(landscape_df$mean)
  names(fitness) <- labels
  start <- which(as.character(landscape_df$support_tier) %in% start_tiers)
  if (!length(start)) start <- which(landscape_df$fq %in% TRUE)
  if (!length(start)) return(data.frame())

  rows <- vector("list", length(start) * n_chr * 2L)
  idx <- 0L
  for (from in start) {
    parent <- k_mat[from, ]
    for (chr in seq_len(n_chr)) {
      for (direction in c(1L, -1L)) {
        child <- parent
        child[chr] <- child[chr] + direction
        if (any(child < min_cn)) next
        child_label <- paste(child, collapse = ".")
        to <- match(child_label, labels)
        if (is.na(to)) next
        delta <- fitness[to] - fitness[from]
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          parent_k = labels[from],
          child_k = labels[to],
          chr = chr,
          direction = if (direction > 0L) "gain" else "loss",
          move = paste0(chr, if (direction > 0L) "+" else "-"),
          parent_fitness = fitness[from],
          child_fitness = fitness[to],
          delta = delta,
          beneficial = is.finite(delta) && delta > 0,
          child_support_tier = landscape_df$support_tier[to],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!idx) return(data.frame())
  do.call(rbind, rows[seq_len(idx)])
}

beneficial_cnv_profile_alfak2 <- function(landscape_df) {
  n_chr <- ncol(parse_karyotype_matrix_bench(landscape_df$k))
  levels <- one_step_move_levels(n_chr)
  edge <- one_step_edge_table(landscape_df)
  prop <- rep(NA_real_, length(levels))
  valid_n <- beneficial_n <- integer(length(levels))
  names(prop) <- names(valid_n) <- names(beneficial_n) <- levels
  if (nrow(edge)) {
    for (mv in levels) {
      ii <- edge$move == mv & is.finite(edge$delta)
      valid_n[mv] <- sum(ii)
      beneficial_n[mv] <- sum(edge$beneficial[ii] %in% TRUE)
      if (valid_n[mv] > 0L) prop[mv] <- beneficial_n[mv] / valid_n[mv]
    }
  }
  list(proportion = prop, valid_n = valid_n, beneficial_n = beneficial_n, edges = edge)
}

beneficial_matrices_from_bundles <- function(bundle_list) {
  profiles <- lapply(bundle_list, function(x) beneficial_cnv_profile_alfak2(x$landscape))
  prop <- do.call(rbind, lapply(profiles, `[[`, "proportion"))
  valid_n <- do.call(rbind, lapply(profiles, `[[`, "valid_n"))
  beneficial_n <- do.call(rbind, lapply(profiles, `[[`, "beneficial_n"))
  rownames(prop) <- rownames(valid_n) <- rownames(beneficial_n) <- names(bundle_list)
  list(
    proportion = prop,
    valid_n = valid_n,
    beneficial_n = beneficial_n,
    edges = do.call(rbind, Map(function(pid, prof) {
      if (!nrow(prof$edges)) return(data.frame())
      cbind(patient_id = pid, prof$edges, stringsAsFactors = FALSE)
    }, names(profiles), profiles))
  )
}

matrix_to_plot_df <- function(mat, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(df) <- c("patient_id", "move", value_name)
  move <- as.character(df$move)
  chr <- suppressWarnings(as.integer(sub("[+-]$", "", move)))
  sign <- sub("^[0-9]+", "", move)
  df$move <- factor(move, levels = one_step_move_levels(max(chr, na.rm = TRUE)))
  df$patient_id <- factor(as.character(df$patient_id), levels = sort_pid_levels_bench(df$patient_id))
  df
}

plot_beneficial_matrix <- function(beneficial_mat) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  df <- matrix_to_plot_df(beneficial_mat, "beneficial_proportion")
  ggplot2::ggplot(df, ggplot2::aes(.data$move, .data$patient_id, fill = .data$beneficial_proportion)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient2(
      low = "#2B6CB0", mid = "#F7F7F7", high = "#C53030",
      midpoint = 0.5, limits = c(0, 1), na.value = "#BDBDBD"
    ) +
    ggplot2::labs(x = "One-step CN event", y = "Patient", fill = "Beneficial\nproportion") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

plot_valid_matrix <- function(valid_mat) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  df <- matrix_to_plot_df(valid_mat, "valid_n")
  ggplot2::ggplot(df, ggplot2::aes(.data$move, .data$patient_id, fill = .data$valid_n)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "#BDBDBD") +
    ggplot2::labs(x = "One-step CN event", y = "Patient", fill = "Valid edges") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

plot_landscape_distribution <- function(bundle_list) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package `ggplot2` is required.", call. = FALSE)
  df <- do.call(rbind, lapply(bundle_list, function(x) x$landscape))
  ggplot2::ggplot(df, ggplot2::aes(.data$mean, fill = .data$support_tier)) +
    ggplot2::geom_histogram(bins = 35, alpha = 0.75, position = "identity") +
    ggplot2::facet_wrap(~patient_id, scales = "free_y") +
    ggplot2::labs(x = "Posterior mean fitness", y = "Node count", fill = "Support tier") +
    ggplot2::theme_bw()
}

write_matrix_csv <- function(mat, path) {
  utils::write.csv(data.frame(patient_id = rownames(mat), mat, check.names = FALSE), path, row.names = FALSE)
  invisible(path)
}
