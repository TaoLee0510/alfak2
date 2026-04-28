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

fit_real_patient_alfak2 <- function(repo_dir,
                                    patient_id,
                                    min_total_count = 20,
                                    drop_diploid = TRUE,
                                    beta = 0.01,
                                    min_cn = 0,
                                    max_cn = 5,
                                    local_shell_depth = 0,
                                    global_extra_shell = 1,
                                    max_nodes = 4000,
                                    force = FALSE,
                                    control = list(eval.max = 250, iter.max = 250)) {
  dirs <- ensure_benchmark_dirs(repo_dir, "real_samples")
  fit_dir <- file.path(dirs$fits, patient_id)
  dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)
  fit_path <- file.path(fit_dir, "alfak2_fit.rds")
  landscape_path <- file.path(fit_dir, "landscape.rds")
  counts_path <- file.path(fit_dir, "counts.csv")

  if (!force && file.exists(fit_path) && file.exists(landscape_path)) {
    return(readRDS(file.path(fit_dir, "bundle.rds")))
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
  utils::write.csv(landscape, file.path(fit_dir, "landscape.csv"), row.names = FALSE)
  saveRDS(fit, fit_path)
  saveRDS(landscape, landscape_path)
  bundle <- list(
    patient_id = patient_id,
    input = built,
    fit = fit,
    landscape = landscape,
    paths = list(fit = fit_path, landscape = landscape_path, counts = counts_path)
  )
  saveRDS(bundle, file.path(fit_dir, "bundle.rds"))
  bundle
}

run_real_sample_landscapes <- function(repo_dir,
                                       patient_subset = NULL,
                                       continue_on_error = TRUE,
                                       ...) {
  data_dir <- file.path(repo_dir, "benchmark", "data")
  pids <- available_real_patient_ids(data_dir)
  if (!is.null(patient_subset) && length(patient_subset)) {
    pids <- intersect(sort_pid_levels_bench(patient_subset), pids)
  }
  if (!length(pids)) stop("No patient ids selected.", call. = FALSE)
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
