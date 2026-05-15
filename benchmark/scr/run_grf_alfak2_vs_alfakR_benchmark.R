#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Run a cross-repo GRF benchmark comparing alfak2 against alfakR nn_prior modes.\n\n",
    "This script is intentionally standalone and is not part of the alfak2 package API.\n",
    "It loads both repositories with pkgload::load_all() so fits use the current local source.\n\n",
    "Usage:\n",
    "  Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R [options]\n\n",
    "Core options:\n",
    "  --mode=all|prepare|fit-task|summarize\n",
    "  --task-index=1                 # 1-based row from task_table.tsv for --mode=fit-task\n",
    "  --alfak2-repo=/share/lab_crd/lab_crd/taoli/Project/alfak2\n",
    "  --alfakR-repo=/share/lab_crd/lab_crd/taoli/Project/alfakR\n",
    "  --output-dir=benchmark/results/grf_alfak2_vs_alfakR\n",
    "  --source-input-dir=benchmark/results/grf_alfak2_vs_alfakR_shared_inputs  # optional prepared input reuse\n",
    "  --path-map=/old/root=/new/root  # optional comma-separated runtime path remapping\n",
    "  --write-node-accuracy=true      # set false for memory-light final summary tables\n",
    "  --methods=none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step\n",
    "  --minobs=5,10,20\n",
    "  --n-sim=1\n",
    "  --lambdas=0.8\n",
    "  --time-gaps=2,4\n",
    "  --time-starts=0\n",
    "  --nboot=5\n",
    "  --sample-depth=2000\n",
    "  --n-cores=1\n\n",
    "GRF/ABM options matching the alfakR GRF benchmark:\n",
    "  --seed=424242\n",
    "  --pm=5e-05                    # data-generating beta used in the ABM simulation\n",
    "  --beta-levels=1e-05,5e-05,1e-04,1e-03,1e-02  # fitted beta grid\n",
    "  --k-dim=22\n",
    "  --n-centroids=64\n",
    "  --time-max=360\n",
    "  --passage-interval=45\n",
    "  --abm-pop-size=50000\n",
    "  --abm-delta-t=1\n",
    "  --abm-max-pop=2000000\n",
    "  --abm-culling-survival=0.01\n\n",
    "alfak2 options:\n",
    "  --alfak2-input-policies=full,minobs_matched,soft_minobs\n",
    "  --alfak2-input-depth=effective\n",
    "  --alfak2-effective-depth-mode=min\n",
    "  --alfak2-graph-edge-weight=mutation|unit|normalized\n",
    "  --alfak2-anchor-count-reference=minobs|none|<number>\n",
    "  --alfak2-local-shell-depth=0\n",
    "  --alfak2-global-extra-shell=1\n",
    "  --alfak2-max-nodes=150000\n\n",
    "Cache options:\n",
    "  --force-refit=false\n",
    "  --force-sim=false\n",
    "  --reuse-dirty-cache=false\n",
    "  --recompile-dll=false          # force native-code rebuild before loading source trees\n",
    sep = ""
  )
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (identical(arg, "--help") || identical(arg, "-h")) {
      out$help <- TRUE
      next
    }
    if (!grepl("^--", arg)) stop("Unexpected positional argument: ", arg, call. = FALSE)
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      val <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      val <- "true"
    }
    key <- gsub("-", "_", key, fixed = TRUE)
    out[[key]] <- val
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || !length(value) || !nzchar(as.character(value[[1L]]))) default else value[[1L]]
}

arg_logical <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  value <- tolower(trimws(as.character(value)))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Expected a boolean for --", gsub("_", "-", name), call. = FALSE)
}

arg_numeric <- function(args, name, default) {
  value <- suppressWarnings(as.numeric(arg_value(args, name, default)))
  if (!is.finite(value)) default else value
}

arg_integer <- function(args, name, default) {
  value <- suppressWarnings(as.integer(arg_value(args, name, default)))
  if (!is.finite(value)) default else value
}

arg_numeric_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  out <- suppressWarnings(as.numeric(strsplit(as.character(value), ",", fixed = TRUE)[[1L]]))
  out <- out[is.finite(out)]
  if (!length(out)) default else out
}

arg_integer_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  out <- suppressWarnings(as.integer(strsplit(as.character(value), ",", fixed = TRUE)[[1L]]))
  out <- out[is.finite(out)]
  if (!length(out)) default else out
}

arg_character_vec <- function(args, name, default) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  out <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  out <- out[nzchar(out)]
  if (!length(out)) default else out
}

resolve_repo_dir <- function(start = getwd()) {
  candidates <- unique(normalizePath(
    file.path(start, c(".", "..", "../..", "../../..")),
    winslash = "/",
    mustWork = FALSE
  ))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) &&
        dir.exists(file.path(cand, "benchmark"))) {
      return(cand)
    }
  }
  stop("Could not locate repository root.", call. = FALSE)
}

normalize_output_dir <- function(repo_dir, output_dir) {
  if (grepl("^/", output_dir)) normalizePath(output_dir, winslash = "/", mustWork = FALSE) else
    normalizePath(file.path(repo_dir, output_dir), winslash = "/", mustWork = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  invisible(path)
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing TSV file: ", path, call. = FALSE)
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )
}

parse_path_map <- function(args) {
  value <- arg_value(args, "path_map", Sys.getenv("ALFAK_BENCH_PATH_MAP", ""))
  if (is.null(value) || !nzchar(as.character(value))) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  entries <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  entries <- entries[nzchar(entries)]
  if (!length(entries)) {
    return(data.frame(from = character(), to = character(), stringsAsFactors = FALSE))
  }
  rows <- lapply(entries, function(entry) {
    eq <- regexpr("=", entry, fixed = TRUE)[[1L]]
    if (eq < 2L) stop("Invalid --path-map entry: ", entry, call. = FALSE)
    data.frame(
      from = substr(entry, 1L, eq - 1L),
      to = substr(entry, eq + 1L, nchar(entry)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(nchar(out$from), decreasing = TRUE), , drop = FALSE]
}

apply_path_map_chr <- function(x, path_map) {
  if (is.null(path_map) || !nrow(path_map) || !length(x)) return(x)
  out <- as.character(x)
  na <- is.na(out)
  out[na] <- ""
  for (i in seq_len(nrow(path_map))) {
    from <- path_map$from[[i]]
    to <- path_map$to[[i]]
    hit <- startsWith(out, from)
    out[hit] <- paste0(to, substring(out[hit], nchar(from) + 1L))
  }
  out[na] <- NA_character_
  out
}

path_columns <- function(x) {
  intersect(
    c(
      "repo_dir", "alfak2_repo", "alfakR_repo", "output_dir", "source_input_dir",
      "outdir", "grf_rds", "input_rds", "fit_path", "landscape_path",
      "bootstrap_path", "posterior_path", "xval_path"
    ),
    names(x)
  )
}

apply_path_map_df <- function(x, path_map) {
  if (is.null(x) || !length(x) || is.null(path_map) || !nrow(path_map)) return(x)
  for (nm in path_columns(x)) {
    if (is.character(x[[nm]]) || is.factor(x[[nm]])) {
      x[[nm]] <- apply_path_map_chr(x[[nm]], path_map)
    }
  }
  x
}

apply_path_map_list <- function(x, path_map) {
  if (is.null(x) || !length(x) || is.null(path_map) || !nrow(path_map)) return(x)
  for (nm in path_columns(x)) {
    if (is.character(x[[nm]]) && length(x[[nm]]) == 1L) {
      x[[nm]] <- apply_path_map_chr(x[[nm]], path_map)
    }
  }
  x
}

read_rds_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

scalar_or_na <- function(x, default = NA_real_) {
  if (is.null(x) || !length(x)) return(default)
  x[[1L]]
}

arg_anchor_count_reference <- function(args, name = "alfak2_anchor_count_reference", default = "minobs") {
  value <- arg_value(args, name, default)
  if (is.null(value)) return(NULL)
  value <- trimws(as.character(value))
  if (!nzchar(value) || tolower(value) %in% c("none", "null", "na", "false")) return(NULL)
  if (tolower(value) %in% c("minobs", "auto")) return("minobs")
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!is.finite(numeric_value) || numeric_value <= 0) {
    stop("`--", gsub("_", "-", name), "` must be minobs, none, or a positive number.", call. = FALSE)
  }
  numeric_value
}

resolve_anchor_count_reference <- function(cfg, task) {
  value <- cfg$alfak2_anchor_count_reference
  if (is.null(value)) return(NULL)
  if (identical(as.character(task$input_policy[[1L]]), "soft_minobs")) {
    if (identical(value, "minobs")) return(as.numeric(task$minobs[[1L]]))
    return(as.numeric(value))
  }
  NULL
}

row_field <- function(row, name, default = NA) {
  if (!name %in% names(row)) return(default)
  value <- row[[name]]
  if (is.null(value) || !length(value)) return(default)
  value[[1L]]
}

validate_beta_grid <- function(x, name = "beta_levels") {
  x <- sort(unique(as.numeric(x)))
  x <- x[is.finite(x)]
  if (!length(x) || any(x < 0 | x > 1)) {
    stop("`", name, "` must contain values in [0, 1].", call. = FALSE)
  }
  x
}

system_text <- function(cmd, args, cwd = NULL) {
  old <- NULL
  if (!is.null(cwd)) {
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    setwd(cwd)
  }
  out <- tryCatch(
    suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = TRUE)),
    error = function(e) character(0)
  )
  paste(out, collapse = "\n")
}

repo_state <- function(repo_dir, package) {
  repo_dir <- normalizePath(repo_dir, winslash = "/", mustWork = FALSE)
  branch <- system_text("git", c("branch", "--show-current"), cwd = repo_dir)
  head <- system_text("git", c("rev-parse", "HEAD"), cwd = repo_dir)
  status <- system_text("git", c("status", "--short"), cwd = repo_dir)
  data.frame(
    package = package,
    repo_dir = repo_dir,
    branch = if (nzchar(branch)) branch else NA_character_,
    head = if (nzchar(head)) head else NA_character_,
    dirty = nzchar(status),
    status_short = if (nzchar(status)) gsub("\n", " | ", status, fixed = TRUE) else "",
    stringsAsFactors = FALSE
  )
}

compile_repo_dll <- function(repo_dir, package) {
  if (!requireNamespace("pkgbuild", quietly = TRUE)) {
    stop("Package `pkgbuild` is required to recompile native code for ", package, ".", call. = FALSE)
  }
  if (!pkgbuild::pkg_has_src(repo_dir)) return(invisible(FALSE))
  message("Recompiling native code for ", package, " on this machine.")
  pkgbuild::clean_dll(repo_dir)
  pkgbuild::compile_dll(repo_dir, force = TRUE, quiet = TRUE)
  invisible(TRUE)
}

load_current_repos <- function(alfakR_repo, alfak2_repo, recompile_dll = FALSE) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package `pkgload` is required so benchmark uses current source trees.", call. = FALSE)
  }
  if (isTRUE(recompile_dll)) {
    compile_repo_dll(alfakR_repo, "alfakR")
    compile_repo_dll(alfak2_repo, "alfak2")
  }
  pkgload::load_all(alfakR_repo, quiet = TRUE, compile = FALSE, recompile = FALSE)
  pkgload::load_all(alfak2_repo, quiet = TRUE, compile = FALSE, recompile = FALSE)
  invisible(TRUE)
}

format_grf_label <- function(x) {
  x <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  gsub("[^A-Za-z0-9]+", "p", x)
}

pm_to_label <- function(x) {
  x <- format(as.numeric(x), scientific = TRUE, digits = 12)
  gsub("[^A-Za-z0-9]+", "_", x)
}

path_token <- function(x) {
  x <- as.character(x)
  x <- gsub("-", "m", x, fixed = TRUE)
  gsub("[^A-Za-z0-9]+", "p", x)
}

parse_karyotype_ids_base <- function(ids) {
  ids <- as.character(ids)
  if (!length(ids) || any(!nzchar(ids))) stop("Karyotype ids must be non-empty.", call. = FALSE)
  pieces <- strsplit(ids, ".", fixed = TRUE)
  lens <- lengths(pieces)
  if (length(unique(lens)) != 1L) stop("Karyotypes have inconsistent dimensions.", call. = FALSE)
  mat <- matrix(as.integer(unlist(pieces, use.names = FALSE)), nrow = length(ids), byrow = TRUE)
  if (anyNA(mat)) stop("Karyotype ids contain non-integer fields.", call. = FALSE)
  rownames(mat) <- ids
  mat
}

format_karyotypes_base <- function(mat) {
  apply(as.matrix(mat), 1L, paste, collapse = ".")
}

compute_grf_fitness_truth <- function(karyotypes, centroids, lambda) {
  karyotypes <- unique(as.character(karyotypes))
  karyotypes <- karyotypes[nzchar(karyotypes)]
  if (!length(karyotypes)) return(stats::setNames(numeric(0), character(0)))
  if (!is.matrix(centroids) || !is.numeric(centroids) || !nrow(centroids)) {
    stop("`centroids` must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("`lambda` must be a positive finite scalar.", call. = FALSE)
  }
  k_mat <- parse_karyotype_ids_base(karyotypes)
  if (ncol(k_mat) != ncol(centroids)) {
    stop("Karyotype dimension does not match GRF centroid dimension.", call. = FALSE)
  }
  out <- vapply(seq_len(nrow(k_mat)), function(i) {
    diffs <- sweep(centroids, 2L, as.numeric(k_mat[i, ]), FUN = "-")
    distances <- sqrt(rowSums(diffs^2))
    sum(sin(distances / lambda)) / (pi * sqrt(nrow(centroids)))
  }, numeric(1))
  stats::setNames(out, rownames(k_mat))
}

make_nn_grf_initial_x0 <- function(k_dim = 22L) {
  k_dim <- as.integer(k_dim)
  if (!is.finite(k_dim) || k_dim < 2L) stop("`k_dim` must be an integer >= 2.", call. = FALSE)
  base <- rep(2L, k_dim)
  states <- list()
  for (chr in seq_len(min(k_dim, 8L))) {
    state <- base
    state[[chr]] <- if (chr %% 3L == 0L) 1L else 3L
    states[[length(states) + 1L]] <- state
  }
  if (k_dim >= 4L) {
    state <- base
    state[1:2] <- c(3L, 3L)
    states[[length(states) + 1L]] <- state
    state <- base
    state[3:4] <- c(1L, 3L)
    states[[length(states) + 1L]] <- state
  }
  ids <- unique(vapply(states, function(v) paste(v, collapse = "."), character(1)))
  weights <- rev(seq_along(ids))
  weights <- weights / sum(weights)
  stats::setNames(as.numeric(weights), ids)
}

make_nn_grf_centroids <- function(lambda, initial_karyotypes, n_centroids, jitter_sd = NULL) {
  initial_mat <- parse_karyotype_ids_base(initial_karyotypes)
  k_dim <- ncol(initial_mat)
  n_centroids <- as.integer(n_centroids)
  if (!is.finite(n_centroids) || n_centroids < 1L) {
    stop("`n_centroids` must be a positive integer.", call. = FALSE)
  }
  if (is.null(jitter_sd)) jitter_sd <- min(0.15, lambda * 0.10)
  centroids <- matrix(NA_real_, nrow = n_centroids, ncol = k_dim)
  for (i in seq_len(n_centroids)) {
    source_idx <- ((i - 1L) %% nrow(initial_mat)) + 1L
    state <- as.numeric(initial_mat[source_idx, ])
    chr <- sample.int(k_dim, 1L)
    direction <- sample(c(-1, 1), 1L)
    shift <- lambda * pi / 2
    if (state[[chr]] + direction * shift < 0.25) direction <- 1
    state[[chr]] <- state[[chr]] + direction * shift
    if (is.finite(jitter_sd) && jitter_sd > 0) {
      state <- state + stats::rnorm(k_dim, mean = 0, sd = jitter_sd)
    }
    centroids[i, ] <- pmax(0.25, state)
  }
  centroids
}

simulate_nn_prior_grf_abm <- function(seed,
                                      lambda,
                                      p,
                                      k_dim,
                                      n_centroids,
                                      time_max,
                                      passage_interval,
                                      abm_pop_size,
                                      abm_delta_t,
                                      abm_max_pop,
                                      abm_culling_survival) {
  set.seed(seed)
  x0 <- make_nn_grf_initial_x0(k_dim = k_dim)
  centroids <- make_nn_grf_centroids(
    lambda = lambda,
    initial_karyotypes = names(x0),
    n_centroids = n_centroids
  )
  sim_times <- seq(0, time_max, by = passage_interval)
  if (length(sim_times) < 2L) stop("GRF simulation needs at least two requested timepoints.", call. = FALSE)
  record_interval <- max(1L, as.integer(round(passage_interval / abm_delta_t)))
  sim_wide <- suppressMessages(
    alfakR::run_abm_simulation_grf(
      centroids = centroids,
      lambda = lambda,
      p = p,
      times = sim_times,
      x0 = x0,
      abm_pop_size = abm_pop_size,
      abm_delta_t = abm_delta_t,
      abm_max_pop = abm_max_pop,
      abm_culling_survival = abm_culling_survival,
      abm_record_interval = record_interval,
      abm_seed = seed,
      normalize_freq = FALSE
    )
  )
  list(
    sim_wide = as.data.frame(sim_wide, check.names = FALSE),
    centroids = centroids,
    lambda = lambda,
    p = p,
    x0 = x0,
    seed = seed
  )
}

select_closest_time_row <- function(sim_times, target) {
  which.min(abs(as.numeric(sim_times) - as.numeric(target)))
}

build_two_timepoint_yi_from_abm <- function(sim_wide,
                                            time_start,
                                            time_gap,
                                            passage_interval,
                                            sample_depth,
                                            seed) {
  if (is.null(sim_wide) || !is.data.frame(sim_wide) || !"time" %in% names(sim_wide)) {
    stop("`sim_wide` must be a data frame with a `time` column.", call. = FALSE)
  }
  sample_depth <- as.integer(sample_depth)
  if (!is.finite(sample_depth) || sample_depth < 1L) {
    stop("`sample_depth` must be a positive integer.", call. = FALSE)
  }
  requested_times <- c(as.numeric(time_start), as.numeric(time_start) + as.numeric(time_gap) * as.numeric(passage_interval))
  row_idx <- vapply(requested_times, function(tt) select_closest_time_row(sim_wide$time, tt), integer(1))
  if (length(unique(row_idx)) != 2L) {
    stop("Selected two-timepoint input collapsed to one ABM row; increase `time_gap`.", call. = FALSE)
  }
  selected_times <- as.numeric(sim_wide$time[row_idx])
  count_cols <- setdiff(names(sim_wide), "time")
  if (!length(count_cols)) stop("ABM output contains no karyotype columns.", call. = FALSE)
  set.seed(seed)
  count_mat <- matrix(
    0L,
    nrow = length(count_cols),
    ncol = 2L,
    dimnames = list(count_cols, format(selected_times, scientific = FALSE, trim = TRUE))
  )
  for (j in seq_along(row_idx)) {
    counts <- suppressWarnings(as.numeric(sim_wide[row_idx[[j]], count_cols, drop = TRUE]))
    counts[!is.finite(counts) | counts < 0] <- 0
    if (sum(counts) <= 0) stop("ABM output has zero mass at a selected time.", call. = FALSE)
    count_mat[, j] <- as.integer(stats::rmultinom(1L, size = sample_depth, prob = counts / sum(counts))[, 1L])
  }
  count_mat <- count_mat[rowSums(count_mat, na.rm = TRUE) > 0, , drop = FALSE]
  if (!nrow(count_mat)) stop("Sampled two-timepoint count matrix is empty.", call. = FALSE)
  list(
    x = count_mat,
    dt = 1,
    metadata = list(
      requested_times = requested_times,
      selected_times = selected_times,
      time_gap = as.numeric(time_gap),
      time_delta = diff(selected_times),
      sample_depth = sample_depth,
      seed = seed
    )
  )
}

drop_diploid_counts <- function(counts) {
  if (!nrow(counts)) return(counts)
  k_dim <- length(strsplit(rownames(counts)[[1]], ".", fixed = TRUE)[[1]])
  diploid <- paste(rep(2L, k_dim), collapse = ".")
  counts[rownames(counts) != diploid, , drop = FALSE]
}

prepare_alfakR_yi <- function(yi, drop_diploid = TRUE) {
  counts <- as.matrix(yi$x)
  storage.mode(counts) <- "integer"
  if (isTRUE(drop_diploid)) counts <- drop_diploid_counts(counts)
  if (!nrow(counts)) stop("No rows remain after diploid filtering.", call. = FALSE)
  out <- yi
  out$x <- counts
  out
}

prepare_alfak2_counts <- function(yi, minobs, input_policy, drop_diploid = TRUE) {
  counts <- as.matrix(yi$x)
  storage.mode(counts) <- "integer"
  if (isTRUE(drop_diploid)) counts <- drop_diploid_counts(counts)
  if (identical(input_policy, "minobs_matched")) {
    counts <- counts[rowSums(counts, na.rm = TRUE) >= as.integer(minobs), , drop = FALSE]
  } else if (!identical(input_policy, "full") && !identical(input_policy, "soft_minobs")) {
    stop("Unsupported alfak2 input policy `", input_policy, "`.", call. = FALSE)
  }
  if (!nrow(counts)) stop("No rows remain for alfak2 input policy `", input_policy, "`.", call. = FALSE)
  if (identical(input_policy, "soft_minobs")) {
    row_totals <- rowSums(counts, na.rm = TRUE)
    weights <- pmin(1, pmax(0, row_totals) / as.numeric(minobs))
    weights[!is.finite(weights)] <- 0
    weights <- cbind(t0 = weights, t1 = weights)
    rownames(weights) <- rownames(counts)
    attr(counts, "observation_weights") <- weights
    attr(counts, "soft_minobs") <- list(minobs = as.integer(minobs), rule = "row_total_over_minobs")
  }
  counts
}

summarize_input_rows <- function(yi, minobs_values, drop_diploid = TRUE) {
  counts <- as.matrix(yi$x)
  raw_rows <- nrow(counts)
  if (isTRUE(drop_diploid)) counts <- drop_diploid_counts(counts)
  row_totals <- rowSums(counts, na.rm = TRUE)
  data.frame(
    minobs = as.integer(minobs_values),
    raw_input_rows = as.integer(raw_rows),
    input_rows_after_drop = as.integer(length(row_totals)),
    input_rows_minobs = as.integer(vapply(minobs_values, function(m) sum(row_totals >= m), integer(1))),
    stringsAsFactors = FALSE
  )
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 2L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

safe_centered_r2 <- function(pred, truth) {
  ok <- is.finite(pred) & is.finite(truth)
  if (sum(ok) < 2L) return(NA_real_)
  truth_c <- truth[ok] - mean(truth[ok])
  pred_c <- pred[ok] - mean(pred[ok])
  den <- sum(truth_c^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  1 - sum((pred_c - truth_c)^2) / den
}

support_scope_for_node <- function(engine, support_tier, fq, nn) {
  if (identical(engine, "alfak2")) {
    ifelse(
      support_tier == "directly_informed",
      "direct",
      ifelse(support_tier %in% c("local_borrowed", "weakly_supported", "graph_borrowed", "prior_dominated"), "nn", "other")
    )
  } else {
    ifelse(fq %in% TRUE, "direct", ifelse(nn %in% TRUE, "nn", "other"))
  }
}

coerce_alfakR_landscape_nodes <- function(landscape,
                                          method,
                                          nn_prior,
                                          fit_row,
                                          centroids,
                                          lambda) {
  if (is.null(landscape) || !nrow(landscape)) return(data.frame())
  x <- as.data.frame(landscape, stringsAsFactors = FALSE)
  if (!"k" %in% names(x) && "karyotype" %in% names(x)) x$k <- x$karyotype
  if (!"mean" %in% names(x) && "fitness_mean" %in% names(x)) x$mean <- x$fitness_mean
  if (!"sd" %in% names(x) && "fitness_sd" %in% names(x)) x$sd <- x$fitness_sd
  if (!"sd" %in% names(x)) x$sd <- NA_real_
  if (!"fq" %in% names(x)) x$fq <- FALSE
  if (!"nn" %in% names(x)) x$nn <- FALSE
  truth <- compute_grf_fitness_truth(x$k, centroids, lambda)
  support_scope <- support_scope_for_node("alfakR", NA_character_, x$fq, x$nn)
  data.frame(
    simulation_id = fit_row$simulation_id,
    lambda = fit_row$lambda,
    lambda_label = fit_row$lambda_label,
    time_start = fit_row$time_start,
    time_gap = fit_row$time_gap,
    time_delta = fit_row$time_delta,
    minobs = fit_row$minobs,
    sim_pm = row_field(fit_row, "sim_pm", NA_real_),
    pm = row_field(fit_row, "pm", NA_real_),
    fit_beta_label = row_field(fit_row, "fit_beta_label", NA_character_),
    graph_edge_weight = row_field(fit_row, "graph_edge_weight", NA_character_),
    anchor_count_reference = row_field(fit_row, "anchor_count_reference", NA_real_),
    anchor_count_power = row_field(fit_row, "anchor_count_power", NA_real_),
    method = method,
    engine = "alfakR",
    input_policy = "alfakR_minobs_internal",
    nn_prior = nn_prior,
    k = as.character(x$k),
    estimated_fitness = as.numeric(x$mean),
    estimated_sd = as.numeric(x$sd),
    true_fitness = as.numeric(truth[match(as.character(x$k), names(truth))]),
    support_tier = support_scope,
    support_scope = support_scope,
    status = fit_row$status,
    stringsAsFactors = FALSE
  )
}

coerce_alfak2_summary_nodes <- function(s,
                                        method,
                                        input_policy,
                                        fit_row,
                                        centroids,
                                        lambda,
                                        use_legacy_scale = TRUE) {
  est_col <- if (isTRUE(use_legacy_scale) && "fitness_mean_alfakR_scale" %in% names(s)) {
    "fitness_mean_alfakR_scale"
  } else {
    "fitness_mean"
  }
  sd_col <- if (isTRUE(use_legacy_scale) && "fitness_sd_alfakR_scale" %in% names(s)) {
    "fitness_sd_alfakR_scale"
  } else {
    "fitness_sd"
  }
  truth <- compute_grf_fitness_truth(s$karyotype, centroids, lambda)
  support_scope <- support_scope_for_node("alfak2", as.character(s$support_tier), FALSE, FALSE)
  data.frame(
    simulation_id = fit_row$simulation_id,
    lambda = fit_row$lambda,
    lambda_label = fit_row$lambda_label,
    time_start = fit_row$time_start,
    time_gap = fit_row$time_gap,
    time_delta = fit_row$time_delta,
    minobs = fit_row$minobs,
    sim_pm = row_field(fit_row, "sim_pm", NA_real_),
    pm = row_field(fit_row, "pm", NA_real_),
    fit_beta_label = row_field(fit_row, "fit_beta_label", NA_character_),
    graph_edge_weight = row_field(fit_row, "graph_edge_weight", NA_character_),
    anchor_count_reference = row_field(fit_row, "anchor_count_reference", NA_real_),
    anchor_count_power = row_field(fit_row, "anchor_count_power", NA_real_),
    method = method,
    engine = "alfak2",
    input_policy = input_policy,
    nn_prior = "alfak2",
    k = as.character(s$karyotype),
    estimated_fitness = as.numeric(s[[est_col]]),
    estimated_sd = as.numeric(s[[sd_col]]),
    true_fitness = as.numeric(truth[match(as.character(s$karyotype), names(truth))]),
    support_tier = as.character(s$support_tier),
    support_scope = support_scope,
    status = fit_row$status,
    stringsAsFactors = FALSE
  )
}

coerce_alfak2_nodes <- function(fit,
                                method,
                                input_policy,
                                fit_row,
                                centroids,
                                lambda,
                                use_legacy_scale = TRUE) {
  s <- alfak2::summarize_alfak2(fit, layer = "global")
  coerce_alfak2_summary_nodes(
    s = s,
    method = method,
    input_policy = input_policy,
    fit_row = fit_row,
    centroids = centroids,
    lambda = lambda,
    use_legacy_scale = use_legacy_scale
  )
}

summarize_accuracy <- function(node_tbl) {
  if (is.null(node_tbl) || !nrow(node_tbl)) return(data.frame())
  group_cols <- c(
    "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "minobs",
    intersect(c("sim_pm", "pm", "fit_beta_label"), names(node_tbl)),
    "method", "engine", "input_policy", "nn_prior"
  )
  scopes <- c("all", "direct", "nn", "other")
  rows <- list()
  idx <- 0L
  split_key <- interaction(node_tbl[group_cols], drop = TRUE, lex.order = TRUE)
  for (key in levels(split_key)) {
    df0 <- node_tbl[split_key == key, , drop = FALSE]
    if (!nrow(df0)) next
    for (scope in scopes) {
      df <- if (identical(scope, "all")) df0 else df0[df0$support_scope == scope, , drop = FALSE]
      ok <- is.finite(df$estimated_fitness) & is.finite(df$true_fitness)
      err <- df$estimated_fitness[ok] - df$true_fitness[ok]
      pred_c <- df$estimated_fitness[ok] - mean(df$estimated_fitness[ok])
      truth_c <- df$true_fitness[ok] - mean(df$true_fitness[ok])
      centered_err <- pred_c - truth_c
      idx <- idx + 1L
      base <- df0[1L, group_cols, drop = FALSE]
      rows[[idx]] <- data.frame(
        base,
        support_scope = scope,
        n_nodes = nrow(df),
        n_scored = sum(ok),
        observed_node_fraction = mean(df$support_scope == "direct", na.rm = TRUE),
        pearson = safe_cor(df$estimated_fitness, df$true_fitness, "pearson"),
        spearman = safe_cor(df$estimated_fitness, df$true_fitness, "spearman"),
        rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
        mae = if (length(err)) mean(abs(err)) else NA_real_,
        signed_bias = if (length(err)) mean(err) else NA_real_,
        centered_rmse = if (length(centered_err)) sqrt(mean(centered_err^2)) else NA_real_,
        centered_mae = if (length(centered_err)) mean(abs(centered_err)) else NA_real_,
        centered_r2 = safe_centered_r2(df$estimated_fitness, df$true_fitness),
        sign_accuracy = if (length(centered_err)) mean(sign(pred_c) == sign(truth_c), na.rm = TRUE) else NA_real_,
        false_high_rate = if (length(centered_err)) mean(pred_c > 0 & truth_c <= 0, na.rm = TRUE) else NA_real_,
        mean_estimated_sd = mean(df$estimated_sd, na.rm = TRUE),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

make_delta_vs_alfakR <- function(summary_tbl) {
  if (is.null(summary_tbl) || !nrow(summary_tbl)) return(data.frame())
  alfak2_rows <- summary_tbl[summary_tbl$engine == "alfak2", , drop = FALSE]
  alfakR_rows <- summary_tbl[summary_tbl$engine == "alfakR", , drop = FALSE]
  if (!nrow(alfak2_rows) || !nrow(alfakR_rows)) return(data.frame())
  keys <- c(
    "lambda", "lambda_label", "time_start", "time_gap", "time_delta", "minobs",
    intersect(c("sim_pm", "pm", "fit_beta_label"), names(summary_tbl)),
    "support_scope"
  )
  rows <- list()
  idx <- 0L
  for (i in seq_len(nrow(alfak2_rows))) {
    a2 <- alfak2_rows[i, , drop = FALSE]
    candidates <- alfakR_rows
    for (k in keys) {
      candidates <- candidates[candidates[[k]] == a2[[k]], , drop = FALSE]
    }
    if (!nrow(candidates)) next
    for (j in seq_len(nrow(candidates))) {
      ar <- candidates[j, , drop = FALSE]
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        a2[, keys, drop = FALSE],
        alfak2_method = a2$method,
        alfak2_input_policy = a2$input_policy,
        alfakR_method = ar$method,
        alfakR_nn_prior = ar$nn_prior,
        delta_centered_rmse = a2$centered_rmse - ar$centered_rmse,
        delta_centered_mae = a2$centered_mae - ar$centered_mae,
        delta_spearman = a2$spearman - ar$spearman,
        delta_pearson = a2$pearson - ar$pearson,
        delta_centered_r2 = a2$centered_r2 - ar$centered_r2,
        delta_false_high_rate = a2$false_high_rate - ar$false_high_rate,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

empty_node_accuracy_table <- function() {
  data.frame(
    simulation_id = integer(),
    lambda = numeric(),
    lambda_label = character(),
    time_start = numeric(),
    time_gap = numeric(),
    time_delta = numeric(),
    minobs = integer(),
    sim_pm = numeric(),
    pm = numeric(),
    fit_beta_label = character(),
    graph_edge_weight = character(),
    anchor_count_reference = numeric(),
    anchor_count_power = numeric(),
    method = character(),
    engine = character(),
    input_policy = character(),
    nn_prior = character(),
    k = character(),
    estimated_fitness = numeric(),
    estimated_sd = numeric(),
    true_fitness = numeric(),
    support_tier = character(),
    support_scope = character(),
    status = character(),
    estimation_error = numeric(),
    stringsAsFactors = FALSE
  )
}

new_accuracy_scope_accumulator <- function() {
  list(
    n_nodes = 0L,
    n_scored = 0L,
    n_direct = 0L,
    estimated = list(),
    truth = list(),
    sd_sum = 0,
    sd_n = 0L
  )
}

new_accuracy_group_accumulator <- function(base, scopes) {
  scope_acc <- stats::setNames(vector("list", length(scopes)), scopes)
  for (scope in scopes) scope_acc[[scope]] <- new_accuracy_scope_accumulator()
  list(base = base, scopes = scope_acc)
}

append_accuracy_scope <- function(scope_acc, df) {
  n <- nrow(df)
  scope_acc$n_nodes <- scope_acc$n_nodes + n
  if (!n) return(scope_acc)
  ok <- is.finite(df$estimated_fitness) & is.finite(df$true_fitness)
  scope_acc$n_scored <- scope_acc$n_scored + sum(ok)
  scope_acc$n_direct <- scope_acc$n_direct + sum(df$support_scope == "direct", na.rm = TRUE)
  if (any(ok)) {
    idx <- length(scope_acc$estimated) + 1L
    scope_acc$estimated[[idx]] <- df$estimated_fitness[ok]
    scope_acc$truth[[idx]] <- df$true_fitness[ok]
  }
  sd_ok <- is.finite(df$estimated_sd)
  if (any(sd_ok)) {
    scope_acc$sd_sum <- scope_acc$sd_sum + sum(df$estimated_sd[sd_ok])
    scope_acc$sd_n <- scope_acc$sd_n + sum(sd_ok)
  }
  scope_acc
}

accuracy_group_key <- function(base) {
  paste(vapply(base, function(x) as.character(x[[1L]]), character(1)), collapse = "\r")
}

accumulate_accuracy_nodes <- function(acc, node_tbl) {
  if (is.null(node_tbl) || !nrow(node_tbl)) return(acc)
  group_cols <- c(
    "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "minobs",
    intersect(c("sim_pm", "pm", "fit_beta_label"), names(node_tbl)),
    "method", "engine", "input_policy", "nn_prior"
  )
  scopes <- c("all", "direct", "nn", "other")
  split_key <- interaction(node_tbl[group_cols], drop = TRUE, lex.order = TRUE)
  for (key in levels(split_key)) {
    df0 <- node_tbl[split_key == key, , drop = FALSE]
    if (!nrow(df0)) next
    base <- df0[1L, group_cols, drop = FALSE]
    group_key <- accuracy_group_key(base)
    if (is.null(acc[[group_key]])) {
      acc[[group_key]] <- new_accuracy_group_accumulator(base, scopes)
    }
    for (scope in scopes) {
      df <- if (identical(scope, "all")) df0 else df0[df0$support_scope == scope, , drop = FALSE]
      acc[[group_key]]$scopes[[scope]] <- append_accuracy_scope(acc[[group_key]]$scopes[[scope]], df)
    }
  }
  acc
}

accuracy_accumulators_to_summary <- function(acc) {
  if (!length(acc)) return(data.frame())
  scopes <- c("all", "direct", "nn", "other")
  rows <- list()
  idx <- 0L
  for (group_key in names(acc)) {
    group <- acc[[group_key]]
    for (scope in scopes) {
      s <- group$scopes[[scope]]
      estimated <- unlist(s$estimated, use.names = FALSE)
      truth <- unlist(s$truth, use.names = FALSE)
      err <- estimated - truth
      pred_c <- estimated - mean(estimated)
      truth_c <- truth - mean(truth)
      centered_err <- pred_c - truth_c
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        group$base,
        support_scope = scope,
        n_nodes = s$n_nodes,
        n_scored = s$n_scored,
        observed_node_fraction = if (s$n_nodes) s$n_direct / s$n_nodes else NaN,
        pearson = safe_cor(estimated, truth, "pearson"),
        spearman = safe_cor(estimated, truth, "spearman"),
        rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
        mae = if (length(err)) mean(abs(err)) else NA_real_,
        signed_bias = if (length(err)) mean(err) else NA_real_,
        centered_rmse = if (length(centered_err)) sqrt(mean(centered_err^2)) else NA_real_,
        centered_mae = if (length(centered_err)) mean(abs(centered_err)) else NA_real_,
        centered_r2 = safe_centered_r2(estimated, truth),
        sign_accuracy = if (length(centered_err)) mean(sign(pred_c) == sign(truth_c), na.rm = TRUE) else NA_real_,
        false_high_rate = if (length(centered_err)) mean(pred_c > 0 & truth_c <= 0, na.rm = TRUE) else NA_real_,
        mean_estimated_sd = if (s$sd_n) s$sd_sum / s$sd_n else NaN,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

read_grf_sim_cached <- function(path, grf_cache = NULL) {
  path <- as.character(path)
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  if (is.environment(grf_cache)) {
    key <- path
    if (exists(key, envir = grf_cache, inherits = FALSE)) {
      return(get(key, envir = grf_cache, inherits = FALSE))
    }
    value <- readRDS(path)
    assign(key, value, envir = grf_cache)
    return(value)
  }
  readRDS(path)
}

read_nodes_for_fit_row <- function(fr, grf_cache = NULL) {
  grf_path <- as.character(fr$grf_rds[[1]])
  grf_sim <- read_grf_sim_cached(grf_path, grf_cache = grf_cache)
  if (is.null(grf_sim)) return(data.frame())
  if (identical(as.character(fr$engine[[1]]), "alfakR")) {
    landscape <- read_rds_if_exists(as.character(fr$landscape_path[[1]]))
    if (is.null(landscape)) return(data.frame())
    coerce_alfakR_landscape_nodes(
      landscape = landscape,
      method = as.character(fr$method[[1]]),
      nn_prior = as.character(fr$nn_prior[[1]]),
      fit_row = fr,
      centroids = grf_sim$centroids,
      lambda = as.numeric(fr$lambda[[1]])
    )
  } else {
    summary <- read_rds_if_exists(as.character(fr$landscape_path[[1]]))
    if (is.data.frame(summary) && "karyotype" %in% names(summary)) {
      return(coerce_alfak2_summary_nodes(
        s = summary,
        method = as.character(fr$method[[1]]),
        input_policy = as.character(fr$input_policy[[1]]),
        fit_row = fr,
        centroids = grf_sim$centroids,
        lambda = as.numeric(fr$lambda[[1]]),
        use_legacy_scale = TRUE
      ))
    }
    fit <- read_rds_if_exists(as.character(fr$fit_path[[1]]))
    if (is.null(fit)) return(data.frame())
    coerce_alfak2_nodes(
      fit = fit,
      method = as.character(fr$method[[1]]),
      input_policy = as.character(fr$input_policy[[1]]),
      fit_row = fr,
      centroids = grf_sim$centroids,
      lambda = as.numeric(fr$lambda[[1]]),
      use_legacy_scale = TRUE
    )
  }
}

summarize_fit_results_streaming <- function(fit_tbl) {
  acc <- list()
  grf_cache <- new.env(parent = emptyenv())
  for (i in seq_len(nrow(fit_tbl))) {
    fr <- fit_tbl[i, , drop = FALSE]
    if (!identical(as.character(fr$status[[1]]), "ok")) next
    if (i %% 25L == 0L) {
      message("  processed fit rows: ", i, " / ", nrow(fit_tbl))
    }
    nodes <- read_nodes_for_fit_row(fr, grf_cache = grf_cache)
    acc <- accumulate_accuracy_nodes(acc, nodes)
  }
  accuracy_accumulators_to_summary(acc)
}

capture_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

run_alfakR_fit <- function(task, cfg, repo_versions) {
  outdir <- task$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  allow_cache <- !isTRUE(cfg$force_refit) &&
    (isTRUE(cfg$reuse_dirty_cache) || !isTRUE(repo_versions$dirty[repo_versions$package == "alfakR"]))
  cached <- if (allow_cache) read_rds_if_exists(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$landscape_path)) {
    cached$cached <- TRUE
    saveRDS(cached, result_path)
    return(cached)
  }
  yi0 <- readRDS(task$input_rds)
  yi <- prepare_alfakR_yi(yi0, drop_diploid = cfg$drop_diploid)
  if (max(rowSums(yi$x), na.rm = TRUE) < task$minobs) {
    res <- as.list(task)
    res$status <- "error"
    res$error_message <- paste0("No frequent karyotypes reach minobs=", task$minobs)
    res$elapsed_sec <- NA_real_
    res$warning_count <- 0L
    res$warning_messages <- NA_character_
    saveRDS(res, result_path)
    return(res)
  }
  warning_log_path <- file.path(outdir, "fit_warnings.log")
  started <- Sys.time()
  res <- tryCatch({
    set.seed(task$benchmark_seed)
    fit_beta <- as.numeric(row_field(task, "pm", cfg$pm))
    cap <- capture_warnings(
      alfakR::alfak(
        yi = yi,
        outdir = outdir,
        passage_times = NULL,
        minobs = task$minobs,
        nboot = cfg$nboot,
        n0 = cfg$n0,
        nb = cfg$nb,
        pm = fit_beta,
        correct_efflux = cfg$correct_efflux,
        nn_prior = task$nn_prior,
        nn_prior_grid_n = cfg$grid_n,
        nn_prior_fit_subset = cfg$nn_prior_fit_subset,
        nn_prior_zero_exposure_quantile = cfg$nn_prior_zero_exposure_quantile,
        nn_prior_zero_weight_scale = cfg$nn_prior_zero_weight_scale,
        nn_prior_zero_weight_cap_ratio = cfg$nn_prior_zero_weight_cap_ratio,
        nn_prior_zero_birth_fallback_weight = cfg$nn_prior_zero_birth_fallback_weight,
        nn_prior_zero_birth_child_floor = cfg$nn_prior_zero_birth_child_floor,
        nn_prior_zero_birth_child_shape = cfg$nn_prior_zero_birth_child_shape,
        nn_prior_zero_birth_replicate_floor = cfg$nn_prior_zero_birth_replicate_floor,
        nn_prior_zero_birth_replicate_shape = cfg$nn_prior_zero_birth_replicate_shape,
        nn_prior_two_step_support = cfg$nn_prior_two_step_support,
        nn_prior_two_step_support_min = cfg$nn_prior_two_step_support_min,
        nn_prior_two_step_cap_floor = cfg$nn_prior_two_step_cap_floor
      )
    )
    if (length(cap$warnings)) writeLines(cap$warnings, warning_log_path)
    xval_path <- file.path(outdir, "xval.Rds")
    xv <- read_rds_if_exists(xval_path)
    c(
      as.list(task),
      list(
        status = "ok",
        cached = FALSE,
        error_message = NA_character_,
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        warning_count = length(cap$warnings),
        warning_messages = if (length(cap$warnings)) paste(cap$warnings, collapse = " || ") else NA_character_,
        landscape_path = file.path(outdir, "landscape.Rds"),
        bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
        posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"),
        xval_path = xval_path,
        xval = if (length(xv) == 1L) as.numeric(xv) else NA_real_
      )
    )
  }, error = function(e) {
    c(
      as.list(task),
      list(
        status = "error",
        cached = FALSE,
        error_message = conditionMessage(e),
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        warning_count = 0L,
        warning_messages = NA_character_,
        landscape_path = file.path(outdir, "landscape.Rds"),
        bootstrap_path = file.path(outdir, "bootstrap_res.Rds"),
        posterior_path = file.path(outdir, "landscape_posterior_samples.Rds"),
        xval_path = file.path(outdir, "xval.Rds"),
        xval = NA_real_
      )
    )
  })
  saveRDS(res, result_path)
  res
}

run_alfak2_fit <- function(task, cfg, repo_versions) {
  outdir <- task$outdir
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(outdir, "fit_result.rds")
  allow_cache <- !isTRUE(cfg$force_refit) &&
    (isTRUE(cfg$reuse_dirty_cache) || !isTRUE(repo_versions$dirty[repo_versions$package == "alfak2"]))
  cached <- if (allow_cache) read_rds_if_exists(result_path) else NULL
  if (is.list(cached) && identical(cached$status, "ok") && file.exists(cached$fit_path)) {
    cached$cached <- TRUE
    saveRDS(cached, result_path)
    return(cached)
  }
  yi <- readRDS(task$input_rds)
  counts <- prepare_alfak2_counts(yi, minobs = task$minobs, input_policy = task$input_policy, drop_diploid = cfg$drop_diploid)
  selected_times <- as.numeric(colnames(counts))
  dt <- diff(selected_times)
  if (!is.finite(dt) || dt <= 0) dt <- as.numeric(task$time_delta)
  max_cn <- cfg$alfak2_max_cn
  if (!is.finite(max_cn)) {
    k_mat <- parse_karyotype_ids_base(rownames(counts))
    max_cn <- max(k_mat, na.rm = TRUE) + cfg$alfak2_local_shell_depth + cfg$alfak2_global_extra_shell
  }
  anchor_count_reference <- resolve_anchor_count_reference(cfg, task)
  started <- Sys.time()
  res <- tryCatch({
    fit_beta <- as.numeric(row_field(task, "pm", cfg$pm))
    fit <- alfak2::fit_alfak2(
      counts,
      dt = dt,
      beta = fit_beta,
      min_cn = cfg$alfak2_min_cn,
      max_cn = as.integer(max_cn),
      local_shell_depth = cfg$alfak2_local_shell_depth,
      global_extra_shell = cfg$alfak2_global_extra_shell,
      max_nodes = cfg$alfak2_max_nodes,
      lambda_l_grid = cfg$alfak2_lambda_l_grid,
      lambda_e_grid = cfg$alfak2_lambda_e_grid,
      sigma_obs_grid = cfg$alfak2_sigma_obs_grid,
      graph_edge_weight = cfg$alfak2_graph_edge_weight,
      anchor_support_tiers = "directly_informed",
      anchor_count_reference = anchor_count_reference,
      anchor_count_power = cfg$alfak2_anchor_count_power,
      input_depth = cfg$alfak2_input_depth,
      effective_depth = cfg$alfak2_effective_depth,
      effective_depth_mode = cfg$alfak2_effective_depth_mode,
      observation_model = cfg$alfak2_observation_model,
      dm_concentration = cfg$alfak2_dm_concentration,
      alfakR_scale = TRUE,
      n0 = cfg$n0,
      nb = cfg$nb,
      correct_efflux = cfg$correct_efflux,
      legacy_weight = cfg$alfak2_legacy_weight,
      control = list(eval.max = cfg$alfak2_eval_max, iter.max = cfg$alfak2_iter_max),
      retry_control = list(eval.max = cfg$alfak2_retry_max, iter.max = cfg$alfak2_retry_max)
    )
    fit_path <- file.path(outdir, "alfak2_fit.rds")
    summary_path <- file.path(outdir, "landscape.rds")
    saveRDS(fit, fit_path)
    saveRDS(alfak2::summarize_alfak2(fit, layer = "global"), summary_path)
    c(
      as.list(task),
      list(
        status = "ok",
        cached = FALSE,
        error_message = NA_character_,
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        fit_path = fit_path,
        landscape_path = summary_path,
        graph_edge_weight = cfg$alfak2_graph_edge_weight,
        anchor_count_reference = if (is.null(anchor_count_reference)) NA_real_ else as.numeric(anchor_count_reference),
        anchor_count_power = cfg$alfak2_anchor_count_power,
        local_convergence = fit$local$diagnostics$convergence,
        local_gradient_norm = fit$local$diagnostics$gradient_norm,
        local_covariance_status = fit$local$diagnostics$covariance_status,
        local_retry_attempted = fit$local$diagnostics$retry_attempted,
        global_factorization_status = fit$global$diagnostics$factorization_status,
        local_nodes = nrow(fit$local$summary),
        global_nodes = nrow(fit$global$summary)
      )
    )
  }, error = function(e) {
    c(
      as.list(task),
      list(
        status = "error",
        cached = FALSE,
        error_message = conditionMessage(e),
        elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
        fit_path = file.path(outdir, "alfak2_fit.rds"),
        landscape_path = file.path(outdir, "landscape.rds"),
        graph_edge_weight = cfg$alfak2_graph_edge_weight,
        anchor_count_reference = if (is.null(anchor_count_reference)) NA_real_ else as.numeric(anchor_count_reference),
        anchor_count_power = cfg$alfak2_anchor_count_power,
        local_convergence = NA_integer_,
        local_gradient_norm = NA_real_,
        local_covariance_status = NA_character_,
        local_retry_attempted = NA,
        global_factorization_status = NA_character_,
        local_nodes = NA_integer_,
        global_nodes = NA_integer_
      )
    )
  })
  saveRDS(res, result_path)
  res
}

run_fit_task <- function(task, cfg, repo_versions) {
  task <- as.list(task)
  if (identical(task$engine, "alfakR")) run_alfakR_fit(task, cfg, repo_versions) else
    run_alfak2_fit(task, cfg, repo_versions)
}

list_to_data_frame <- function(x) {
  if (!length(x)) return(data.frame())
  nms <- unique(unlist(lapply(x, names), use.names = FALSE))
  rows <- lapply(x, function(row) {
    out <- setNames(vector("list", length(nms)), nms)
    for (nm in names(row)) out[[nm]] <- row[[nm]]
    out <- lapply(out, function(value) {
      if (is.null(value)) return(NA)
      if (length(value) == 0L) return(NA)
      if (length(value) > 1L) return(paste(as.character(value), collapse = " || "))
      value
    })
    as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
  })
  do.call(rbind, rows)
}

validate_task_shared_sources <- function(task_tbl) {
  if (is.null(task_tbl) || !nrow(task_tbl)) return(invisible(TRUE))
  key_cols <- c(
    "simulation_id", "lambda_label", "time_start", "time_gap", "time_delta",
    "minobs", intersect(c("pm", "fit_beta_label"), names(task_tbl))
  )
  group_key <- interaction(task_tbl[key_cols], drop = TRUE, lex.order = TRUE)
  groups <- split(task_tbl, group_key)
  bad_input <- names(groups)[vapply(groups, function(df) {
    length(unique(as.character(df$input_rds))) != 1L || length(unique(as.character(df$input_md5))) != 1L
  }, logical(1))]
  bad_truth <- names(groups)[vapply(groups, function(df) {
    length(unique(as.character(df$grf_key))) != 1L || length(unique(as.character(df$grf_rds))) != 1L
  }, logical(1))]
  bad_engines <- names(groups)[vapply(groups, function(df) {
    !all(c("alfak2", "alfakR") %in% as.character(df$engine))
  }, logical(1))]
  if (length(bad_input)) {
    stop("Benchmark tasks do not share the same two-timepoint input within matched conditions.", call. = FALSE)
  }
  if (length(bad_truth)) {
    stop("Benchmark tasks do not share the same GRF truth source within matched conditions.", call. = FALSE)
  }
  if (length(bad_engines)) {
    stop("Each matched condition must contain both alfak2 and alfakR tasks.", call. = FALSE)
  }
  invisible(TRUE)
}

run_fit_task_table <- function(task_tbl, cfg, repo_versions, label) {
  if (is.null(task_tbl) || !nrow(task_tbl)) return(list())
  n_cores <- max(1L, min(as.integer(cfg$n_cores), nrow(task_tbl)))
  message("Running ", label, ": ", nrow(task_tbl), " fit tasks with n_cores=", n_cores, ".")
  task_list <- lapply(seq_len(nrow(task_tbl)), function(i) task_tbl[i, , drop = FALSE])
  if (.Platform$OS.type == "unix" && n_cores > 1L) {
    parallel::mclapply(task_list, run_fit_task, cfg = cfg, repo_versions = repo_versions, mc.cores = n_cores, mc.preschedule = FALSE)
  } else {
    lapply(task_list, run_fit_task, cfg = cfg, repo_versions = repo_versions)
  }
}

build_config <- function(args, repo_dir) {
  methods <- arg_character_vec(args, "methods", c("none", "empirical", "empirical_censored", "empirical_censored_weighted", "empirical_two_step"))
  methods <- sub("^nn_prior_", "", methods)
  methods <- setdiff(methods, "cohort_transition")
  input_policies <- arg_character_vec(args, "alfak2_input_policies", c("full", "minobs_matched"))
  bad_policies <- setdiff(input_policies, c("full", "minobs_matched", "soft_minobs"))
  if (length(bad_policies)) stop("Unsupported alfak2 input policies: ", paste(bad_policies, collapse = ", "), call. = FALSE)
  graph_edge_weight <- as.character(arg_value(args, "alfak2_graph_edge_weight", "mutation"))
  graph_edge_weight <- match.arg(graph_edge_weight, c("mutation", "unit", "normalized"))
  list(
    repo_dir = repo_dir,
    alfak2_repo = normalizePath(arg_value(args, "alfak2_repo", repo_dir), winslash = "/", mustWork = FALSE),
    alfakR_repo = normalizePath(arg_value(args, "alfakR_repo", "/share/lab_crd/lab_crd/taoli/Project/alfakR"), winslash = "/", mustWork = FALSE),
    output_dir = normalize_output_dir(repo_dir, arg_value(args, "output_dir", "benchmark/results/grf_alfak2_vs_alfakR")),
    source_input_dir = {
      x <- arg_value(args, "source_input_dir", NULL)
      if (is.null(x)) NULL else normalize_output_dir(repo_dir, x)
    },
    methods = methods,
    minobs = sort(unique(arg_integer_vec(args, "minobs", c(5L, 10L, 20L)))),
    n_sim = arg_integer(args, "n_sim", 1L),
    lambdas = arg_numeric_vec(args, "lambdas", 0.8),
    time_starts = arg_numeric_vec(args, "time_starts", 0),
    time_gaps = arg_numeric_vec(args, "time_gaps", c(2, 4)),
    nboot = arg_integer(args, "nboot", 5L),
    grid_n = arg_integer(args, "grid_n", 81L),
    n_cores = arg_integer(args, "n_cores", 1L),
    seed = arg_integer(args, "seed", 424242L),
    pm = arg_numeric(args, "pm", 5e-05),
    beta_levels = validate_beta_grid(arg_numeric_vec(args, "beta_levels", c(1e-05, 5e-05, 1e-04, 1e-03, 1e-02))),
    n0 = arg_numeric(args, "n0", 100000),
    nb = arg_numeric(args, "nb", 10000000),
    correct_efflux = arg_logical(args, "correct_efflux", TRUE),
    drop_diploid = arg_logical(args, "drop_diploid", TRUE),
    k_dim = arg_integer(args, "k_dim", 22L),
    n_centroids = arg_integer(args, "n_centroids", 64L),
    time_max = arg_numeric(args, "time_max", 360),
    passage_interval = arg_numeric(args, "passage_interval", 45),
    sample_depth = arg_integer(args, "sample_depth", 2000L),
    abm_pop_size = arg_numeric(args, "abm_pop_size", 50000),
    abm_delta_t = arg_numeric(args, "abm_delta_t", 1),
    abm_max_pop = arg_numeric(args, "abm_max_pop", 2000000),
    abm_culling_survival = arg_numeric(args, "abm_culling_survival", 0.01),
    force_refit = arg_logical(args, "force_refit", FALSE),
    force_sim = arg_logical(args, "force_sim", FALSE),
    reuse_dirty_cache = arg_logical(args, "reuse_dirty_cache", FALSE),
    nn_prior_fit_subset = as.character(arg_value(args, "nn_prior_fit_subset", "hybrid")),
    nn_prior_zero_exposure_quantile = arg_numeric(args, "nn_prior_zero_exposure_quantile", 0.10),
    nn_prior_zero_weight_scale = arg_numeric(args, "nn_prior_zero_weight_scale", 0.50),
    nn_prior_zero_weight_cap_ratio = {
      x <- arg_value(args, "nn_prior_zero_weight_cap_ratio", NA_character_)
      y <- suppressWarnings(as.numeric(x))
      if (is.finite(y)) y else NULL
    },
    nn_prior_zero_birth_fallback_weight = {
      x <- arg_value(args, "nn_prior_zero_birth_fallback_weight", NA_character_)
      y <- suppressWarnings(as.numeric(x))
      if (is.finite(y)) y else NULL
    },
    nn_prior_zero_birth_child_floor = arg_numeric(args, "nn_prior_zero_birth_child_floor", 0.25),
    nn_prior_zero_birth_child_shape = arg_numeric(args, "nn_prior_zero_birth_child_shape", 1),
    nn_prior_zero_birth_replicate_floor = arg_numeric(args, "nn_prior_zero_birth_replicate_floor", 0.50),
    nn_prior_zero_birth_replicate_shape = arg_numeric(args, "nn_prior_zero_birth_replicate_shape", 1),
    nn_prior_two_step_support = as.character(arg_value(args, "nn_prior_two_step_support", "rescue")),
    nn_prior_two_step_support_min = arg_numeric(args, "nn_prior_two_step_support_min", 0.15),
    nn_prior_two_step_cap_floor = arg_numeric(args, "nn_prior_two_step_cap_floor", 0.30),
    alfak2_input_policies = input_policies,
    alfak2_input_depth = as.character(arg_value(args, "alfak2_input_depth", "effective")),
    alfak2_effective_depth = {
      x <- arg_value(args, "alfak2_effective_depth", NA_character_)
      y <- suppressWarnings(as.numeric(x))
      if (is.finite(y)) y else NULL
    },
    alfak2_effective_depth_mode = as.character(arg_value(args, "alfak2_effective_depth_mode", "min")),
    alfak2_observation_model = {
      x <- as.character(arg_value(args, "alfak2_observation_model", ""))
      if (nzchar(x)) x else NULL
    },
    alfak2_dm_concentration = {
      x <- arg_value(args, "alfak2_dm_concentration", NA_character_)
      y <- suppressWarnings(as.numeric(x))
      if (is.finite(y)) y else NULL
    },
    alfak2_min_cn = arg_integer(args, "alfak2_min_cn", 0L),
    alfak2_max_cn = {
      x <- arg_value(args, "alfak2_max_cn", NA_character_)
      y <- suppressWarnings(as.integer(x))
      if (is.finite(y)) y else NA_integer_
    },
    alfak2_local_shell_depth = arg_integer(args, "alfak2_local_shell_depth", 0L),
    alfak2_global_extra_shell = arg_integer(args, "alfak2_global_extra_shell", 1L),
    alfak2_max_nodes = arg_integer(args, "alfak2_max_nodes", 150000L),
    alfak2_lambda_l_grid = arg_numeric_vec(args, "alfak2_lambda_l_grid", 1),
    alfak2_lambda_e_grid = arg_numeric_vec(args, "alfak2_lambda_e_grid", 0.25),
    alfak2_sigma_obs_grid = arg_numeric_vec(args, "alfak2_sigma_obs_grid", 0.05),
    alfak2_graph_edge_weight = graph_edge_weight,
    alfak2_anchor_count_reference = arg_anchor_count_reference(args),
    alfak2_anchor_count_power = arg_numeric(args, "alfak2_anchor_count_power", 1),
    alfak2_legacy_weight = as.character(arg_value(args, "alfak2_legacy_weight", "pi0")),
    alfak2_eval_max = arg_integer(args, "alfak2_eval_max", 500L),
    alfak2_iter_max = arg_integer(args, "alfak2_iter_max", 500L),
    alfak2_retry_max = arg_integer(args, "alfak2_retry_max", 2000L),
    write_node_accuracy = arg_logical(args, "write_node_accuracy", TRUE)
  )
}

build_dirs <- function(output_dir) {
  dirs <- list(
    root = output_dir,
    cache = file.path(output_dir, "cache"),
    fits = file.path(output_dir, "fits"),
    tables = file.path(output_dir, "tables"),
    fit_parts = file.path(output_dir, "tables", "fit_results_parts")
  )
  for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  dirs
}

numeric_in <- function(x, values, tol = 1e-8) {
  x <- as.numeric(x)
  values <- as.numeric(values)
  vapply(x, function(xx) any(abs(xx - values) <= tol), logical(1))
}

resolve_source_cache_path <- function(path, source_dir) {
  path <- as.character(path)
  if (!length(path) || !nzchar(path[[1L]])) return(NA_character_)
  if (file.exists(path)) return(normalizePath(path, winslash = "/", mustWork = TRUE))
  fallback <- file.path(source_dir, "cache", basename(path))
  if (file.exists(fallback)) return(normalizePath(fallback, winslash = "/", mustWork = TRUE))
  path
}

filter_source_input_table <- function(source_tbl, cfg) {
  needed <- c(
    "simulation_id", "lambda", "lambda_label", "time_start", "time_gap", "time_delta",
    "patient_id", "grf_key", "grf_rds", "input_rds", "input_md5", "minobs"
  )
  missing <- setdiff(needed, names(source_tbl))
  if (length(missing)) {
    stop("Source input table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  source_tbl$simulation_id <- as.integer(source_tbl$simulation_id)
  source_tbl$lambda <- as.numeric(source_tbl$lambda)
  source_tbl$time_start <- as.numeric(source_tbl$time_start)
  source_tbl$time_gap <- as.numeric(source_tbl$time_gap)
  source_tbl$time_delta <- as.numeric(source_tbl$time_delta)
  source_tbl$minobs <- as.integer(source_tbl$minobs)
  if ("sim_pm" %in% names(source_tbl)) source_tbl$sim_pm <- as.numeric(source_tbl$sim_pm)
  keep <- source_tbl$simulation_id %in% seq_len(cfg$n_sim) &
    numeric_in(source_tbl$lambda, cfg$lambdas) &
    numeric_in(source_tbl$time_start, cfg$time_starts) &
    numeric_in(source_tbl$time_gap, cfg$time_gaps) &
    source_tbl$minobs %in% cfg$minobs
  if ("sim_pm" %in% names(source_tbl)) {
    keep <- keep & numeric_in(source_tbl$sim_pm, cfg$pm)
  }
  out <- source_tbl[keep, , drop = FALSE]
  out <- out[order(out$simulation_id, out$lambda, out$time_start, out$time_gap, out$minobs), , drop = FALSE]
  if (!nrow(out)) {
    stop(
      "No rows in source input table match requested benchmark grid. ",
      "Check --pm, --n-sim, --lambdas, --time-gaps, --time-starts, and --minobs.",
      call. = FALSE
    )
  }
  out
}

prepare_benchmark_inputs_from_source <- function(cfg, dirs) {
  source_dir <- normalizePath(cfg$source_input_dir, winslash = "/", mustWork = TRUE)
  source_input_path <- file.path(source_dir, "tables", "input_table.tsv")
  if (!file.exists(source_input_path)) {
    stop("Missing source input table: ", source_input_path, call. = FALSE)
  }
  input_tbl <- filter_source_input_table(read_tsv(source_input_path), cfg)
  input_tbl$source_input_dir <- source_dir
  input_tbl$input_source <- "reused"
  input_tbl$grf_rds <- vapply(input_tbl$grf_rds, resolve_source_cache_path, character(1), source_dir = source_dir)
  input_tbl$input_rds <- vapply(input_tbl$input_rds, resolve_source_cache_path, character(1), source_dir = source_dir)
  if ("input_csv" %in% names(input_tbl)) {
    input_tbl$input_csv <- vapply(input_tbl$input_csv, resolve_source_cache_path, character(1), source_dir = source_dir)
  }
  missing_grf <- unique(input_tbl$grf_rds[!file.exists(input_tbl$grf_rds)])
  missing_input <- unique(input_tbl$input_rds[!file.exists(input_tbl$input_rds)])
  if (length(missing_grf) || length(missing_input)) {
    msg <- c(
      if (length(missing_grf)) paste0("missing grf_rds: ", paste(utils::head(missing_grf, 3), collapse = ", ")),
      if (length(missing_input)) paste0("missing input_rds: ", paste(utils::head(missing_input, 3), collapse = ", "))
    )
    stop("Source input files are not available: ", paste(msg, collapse = " | "), call. = FALSE)
  }
  if (!"input_md5" %in% names(input_tbl)) input_tbl$input_md5 <- unname(tools::md5sum(input_tbl$input_rds))

  fit_tasks <- list()
  task_idx <- 0L
  for (i in seq_len(nrow(input_tbl))) {
    row <- input_tbl[i, , drop = FALSE]
    sim_pm <- as.numeric(row_field(row, "sim_pm", cfg$pm))
    for (fit_beta in cfg$beta_levels) {
      fit_beta_label <- pm_to_label(fit_beta)
      for (policy in cfg$alfak2_input_policies) {
        task_idx <- task_idx + 1L
        method <- paste0("alfak2_", cfg$alfak2_input_depth, "_", policy)
        fit_tasks[[task_idx]] <- data.frame(
          engine = "alfak2",
          method = method,
          input_policy = policy,
          nn_prior = NA_character_,
          simulation_id = as.integer(row$simulation_id),
          lambda = as.numeric(row$lambda),
          lambda_label = as.character(row$lambda_label),
          time_start = as.numeric(row$time_start),
          time_gap = as.numeric(row$time_gap),
          time_delta = as.numeric(row$time_delta),
          sim_pm = sim_pm,
          pm = fit_beta,
          fit_beta_label = fit_beta_label,
          patient_id = as.character(row$patient_id),
          grf_key = as.character(row$grf_key),
          grf_rds = as.character(row$grf_rds),
          input_rds = as.character(row$input_rds),
          input_md5 = as.character(row$input_md5),
          minobs = as.integer(row$minobs),
          benchmark_seed = as.integer(cfg$seed + as.integer(row$simulation_id) * 10000L +
                                        round(as.numeric(row$time_gap) * 1000) +
                                        as.integer(row$minobs) * 100L + 900L +
                                        match(policy, cfg$alfak2_input_policies)),
          outdir = file.path(dirs$fits, "alfak2", paste0("lambda_", row$lambda_label),
                             paste0("gap_", path_token(row$time_gap)),
                             paste0("beta_", fit_beta_label), paste0("MINOBS_", row$minobs),
                             method, row$patient_id),
          stringsAsFactors = FALSE
        )
      }
      for (method in cfg$methods) {
        task_idx <- task_idx + 1L
        fit_tasks[[task_idx]] <- data.frame(
          engine = "alfakR",
          method = paste0("alfakR_", method),
          input_policy = "alfakR_minobs_internal",
          nn_prior = method,
          simulation_id = as.integer(row$simulation_id),
          lambda = as.numeric(row$lambda),
          lambda_label = as.character(row$lambda_label),
          time_start = as.numeric(row$time_start),
          time_gap = as.numeric(row$time_gap),
          time_delta = as.numeric(row$time_delta),
          sim_pm = sim_pm,
          pm = fit_beta,
          fit_beta_label = fit_beta_label,
          patient_id = as.character(row$patient_id),
          grf_key = as.character(row$grf_key),
          grf_rds = as.character(row$grf_rds),
          input_rds = as.character(row$input_rds),
          input_md5 = as.character(row$input_md5),
          minobs = as.integer(row$minobs),
          benchmark_seed = as.integer(cfg$seed + as.integer(row$simulation_id) * 10000L +
                                        round(as.numeric(row$time_gap) * 1000) +
                                        as.integer(row$minobs) * 100L +
                                        match(method, cfg$methods)),
          outdir = file.path(dirs$fits, "alfakR", paste0("lambda_", row$lambda_label),
                             paste0("gap_", path_token(row$time_gap)),
                             paste0("beta_", fit_beta_label), paste0("MINOBS_", row$minobs),
                             paste0("nn_prior_", method), row$patient_id),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  task_tbl <- do.call(rbind, fit_tasks)
  task_tbl <- task_tbl[c(which(task_tbl$engine == "alfak2"), which(task_tbl$engine == "alfakR")), , drop = FALSE]
  task_tbl$task_order <- seq_len(nrow(task_tbl))
  validate_task_shared_sources(task_tbl)
  write_tsv(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  saveRDS(cfg, file.path(dirs$root, "benchmark_config.rds"))
  message(
    "Prepared ", nrow(task_tbl), " fit tasks from source inputs: ",
    nrow(input_tbl), " source input rows."
  )
  engine_counts <- table(task_tbl$engine)
  for (engine in names(engine_counts)) {
    message("  ", engine, " tasks: ", unname(engine_counts[[engine]]))
  }
  list(input_table = input_tbl, task_table = task_tbl, grf_lookup = list())
}

prepare_benchmark_inputs <- function(cfg, dirs, repo_versions) {
  if (!is.null(cfg$source_input_dir)) {
    return(prepare_benchmark_inputs_from_source(cfg, dirs))
  }

  input_rows <- list()
  fit_tasks <- list()
  input_idx <- 0L
  task_idx <- 0L
  grf_lookup <- list()
  time_axis_label <- paste0("tmax_", path_token(cfg$time_max), "_pint_", path_token(cfg$passage_interval))
  sim_pm_label <- pm_to_label(cfg$pm)

  for (sim_idx in seq_len(cfg$n_sim)) {
    for (lambda_idx in seq_along(cfg$lambdas)) {
      lambda <- cfg$lambdas[[lambda_idx]]
      lambda_label <- format_grf_label(lambda)
      abm_seed <- cfg$seed + sim_idx * 10000L + lambda_idx * 100L
      grf_key <- paste(sim_idx, lambda_label, paste0("simpm_", sim_pm_label), time_axis_label, sep = "__")
      grf_path <- file.path(dirs$cache, paste0("grf_sim_", grf_key, ".rds"))
      grf_sim <- if (!isTRUE(cfg$force_sim) &&
                     file.exists(grf_path) &&
                     (isTRUE(cfg$reuse_dirty_cache) || !isTRUE(repo_versions$dirty[repo_versions$package == "alfakR"]))) {
        readRDS(grf_path)
      } else {
        message("Simulating GRF ABM: sim=", sim_idx, " lambda=", lambda)
        out <- simulate_nn_prior_grf_abm(
          seed = abm_seed,
          lambda = lambda,
          p = cfg$pm,
          k_dim = cfg$k_dim,
          n_centroids = cfg$n_centroids,
          time_max = cfg$time_max,
          passage_interval = cfg$passage_interval,
          abm_pop_size = cfg$abm_pop_size,
          abm_delta_t = cfg$abm_delta_t,
          abm_max_pop = cfg$abm_max_pop,
          abm_culling_survival = cfg$abm_culling_survival
        )
        saveRDS(out, grf_path)
        out
      }
      grf_lookup[[grf_key]] <- grf_sim

      for (time_start in cfg$time_starts) {
        for (time_gap in cfg$time_gaps) {
          patient_id <- paste0("grf_", sim_idx, "_lambda_", lambda_label, "_simpm_", sim_pm_label, "_", time_axis_label, "_start_", path_token(time_start), "_gap_", path_token(time_gap))
          input_rds <- file.path(dirs$cache, paste0("input_", patient_id, ".rds"))
          input_csv <- file.path(dirs$cache, paste0("input_", patient_id, ".csv"))
          input_seed <- abm_seed + as.integer(round(time_start * 10)) + as.integer(round(time_gap * 1000))
          yi <- if (!file.exists(input_rds) || isTRUE(cfg$force_sim)) {
            out <- build_two_timepoint_yi_from_abm(
              sim_wide = grf_sim$sim_wide,
              time_start = time_start,
              time_gap = time_gap,
              passage_interval = cfg$passage_interval,
              sample_depth = cfg$sample_depth,
              seed = input_seed
            )
            saveRDS(out, input_rds)
            utils::write.csv(data.frame(karyotype = rownames(out$x), out$x, check.names = FALSE), input_csv, row.names = FALSE)
            out
          } else {
            readRDS(input_rds)
          }
          input_md5 <- unname(tools::md5sum(input_rds))
          input_summary <- summarize_input_rows(yi, cfg$minobs, drop_diploid = cfg$drop_diploid)
          for (rr in seq_len(nrow(input_summary))) {
            input_idx <- input_idx + 1L
            input_rows[[input_idx]] <- data.frame(
              simulation_id = sim_idx,
              lambda = lambda,
              lambda_label = lambda_label,
              time_start = time_start,
              time_gap = time_gap,
              time_delta = as.numeric(yi$metadata$time_delta),
              sim_pm = cfg$pm,
              patient_id = patient_id,
              grf_key = grf_key,
              grf_rds = grf_path,
              input_rds = input_rds,
              input_csv = input_csv,
              input_md5 = input_md5,
              input_summary[rr, , drop = FALSE],
              stringsAsFactors = FALSE
            )
          }

          for (minobs in cfg$minobs) {
            for (fit_beta in cfg$beta_levels) {
              fit_beta_label <- pm_to_label(fit_beta)
              for (policy in cfg$alfak2_input_policies) {
                task_idx <- task_idx + 1L
                method <- paste0("alfak2_", cfg$alfak2_input_depth, "_", policy)
                fit_tasks[[task_idx]] <- data.frame(
                  engine = "alfak2",
                  method = method,
                  input_policy = policy,
                  nn_prior = NA_character_,
                  simulation_id = sim_idx,
                  lambda = lambda,
                  lambda_label = lambda_label,
                  time_start = time_start,
                  time_gap = time_gap,
                  time_delta = as.numeric(yi$metadata$time_delta),
                  sim_pm = cfg$pm,
                  pm = fit_beta,
                  fit_beta_label = fit_beta_label,
                  patient_id = patient_id,
                  grf_key = grf_key,
                  grf_rds = grf_path,
                  input_rds = input_rds,
                  input_md5 = input_md5,
                  minobs = as.integer(minobs),
                  benchmark_seed = as.integer(abm_seed + time_gap * 1000L + minobs * 100L + 900L + match(policy, cfg$alfak2_input_policies)),
                  outdir = file.path(dirs$fits, "alfak2", paste0("lambda_", lambda_label), paste0("gap_", path_token(time_gap)), paste0("beta_", fit_beta_label), paste0("MINOBS_", minobs), method, patient_id),
                  stringsAsFactors = FALSE
                )
              }
              for (method in cfg$methods) {
                task_idx <- task_idx + 1L
                fit_tasks[[task_idx]] <- data.frame(
                  engine = "alfakR",
                  method = paste0("alfakR_", method),
                  input_policy = "alfakR_minobs_internal",
                  nn_prior = method,
                  simulation_id = sim_idx,
                  lambda = lambda,
                  lambda_label = lambda_label,
                  time_start = time_start,
                  time_gap = time_gap,
                  time_delta = as.numeric(yi$metadata$time_delta),
                  sim_pm = cfg$pm,
                  pm = fit_beta,
                  fit_beta_label = fit_beta_label,
                  patient_id = patient_id,
                  grf_key = grf_key,
                  grf_rds = grf_path,
                  input_rds = input_rds,
                  input_md5 = input_md5,
                  minobs = as.integer(minobs),
                  benchmark_seed = as.integer(abm_seed + time_gap * 1000L + minobs * 100L + match(method, cfg$methods)),
                  outdir = file.path(dirs$fits, "alfakR", paste0("lambda_", lambda_label), paste0("gap_", path_token(time_gap)), paste0("beta_", fit_beta_label), paste0("MINOBS_", minobs), paste0("nn_prior_", method), patient_id),
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }
      }
    }
  }

  input_tbl <- do.call(rbind, input_rows)
  task_tbl <- do.call(rbind, fit_tasks)
  task_tbl <- task_tbl[c(which(task_tbl$engine == "alfak2"), which(task_tbl$engine == "alfakR")), , drop = FALSE]
  task_tbl$task_order <- seq_len(nrow(task_tbl))
  validate_task_shared_sources(task_tbl)
  write_tsv(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  saveRDS(cfg, file.path(dirs$root, "benchmark_config.rds"))
  message("Prepared ", nrow(task_tbl), " fit tasks.")
  engine_counts <- table(task_tbl$engine)
  for (engine in names(engine_counts)) {
    message("  ", engine, " tasks: ", unname(engine_counts[[engine]]))
  }
  list(input_table = input_tbl, task_table = task_tbl, grf_lookup = grf_lookup)
}

load_prepared_config <- function(cfg, dirs) {
  cfg_path <- file.path(dirs$root, "benchmark_config.rds")
  if (file.exists(cfg_path)) readRDS(cfg_path) else cfg
}

apply_runtime_overrides <- function(cfg, args, repo_dir) {
  if (!is.null(args$alfak2_repo)) {
    cfg$alfak2_repo <- normalizePath(arg_value(args, "alfak2_repo", repo_dir), winslash = "/", mustWork = FALSE)
  }
  if (!is.null(args$alfakR_repo)) {
    cfg$alfakR_repo <- normalizePath(arg_value(args, "alfakR_repo", cfg$alfakR_repo), winslash = "/", mustWork = FALSE)
  }
  if (!is.null(args$output_dir)) {
    cfg$output_dir <- normalize_output_dir(repo_dir, arg_value(args, "output_dir", cfg$output_dir))
  }
  if (!is.null(args$source_input_dir)) {
    cfg$source_input_dir <- normalize_output_dir(repo_dir, arg_value(args, "source_input_dir", cfg$source_input_dir))
  }
  if (!is.null(args$write_node_accuracy)) {
    cfg$write_node_accuracy <- arg_logical(args, "write_node_accuracy", TRUE)
  } else if (is.null(cfg$write_node_accuracy)) {
    cfg$write_node_accuracy <- TRUE
  }
  cfg
}

select_task_row <- function(task_tbl, task_index) {
  task_index <- as.integer(task_index)
  if (!is.finite(task_index) || task_index < 1L) {
    stop("`--task-index` must be a positive 1-based integer.", call. = FALSE)
  }
  if ("task_order" %in% names(task_tbl)) {
    idx <- which(as.integer(task_tbl$task_order) == task_index)
  } else {
    idx <- task_index
  }
  if (!length(idx) || idx[[1]] < 1L || idx[[1]] > nrow(task_tbl)) {
    stop("Task index ", task_index, " is outside the prepared task table.", call. = FALSE)
  }
  task_tbl[idx[[1]], , drop = FALSE]
}

run_one_prepared_task <- function(cfg, dirs, repo_versions, task_index) {
  task_tbl <- read_tsv(file.path(dirs$tables, "task_table.tsv"))
  task <- select_task_row(task_tbl, task_index)
  task <- apply_path_map_df(task, cfg$path_map)
  if (!file.exists(as.character(task$input_rds[[1]]))) {
    stop("Prepared input_rds is missing: ", as.character(task$input_rds[[1]]), call. = FALSE)
  }
  if (!file.exists(as.character(task$grf_rds[[1]]))) {
    stop("Prepared grf_rds is missing: ", as.character(task$grf_rds[[1]]), call. = FALSE)
  }
  actual_md5 <- unname(tools::md5sum(as.character(task$input_rds[[1]])))
  expected_md5 <- as.character(task$input_md5[[1]])
  if (nzchar(expected_md5) && !identical(actual_md5, expected_md5)) {
    stop("Prepared input_rds md5 changed for task ", task_index, ".", call. = FALSE)
  }
  result <- run_fit_task(task, cfg, repo_versions)
  fit_tbl <- list_to_data_frame(list(result))
  part_path <- file.path(dirs$fit_parts, sprintf("task_%06d.tsv", as.integer(task_index)))
  write_tsv(fit_tbl, part_path)
  writeLines("ok", file.path(dirs$fit_parts, sprintf("task_%06d.done", as.integer(task_index))))
  message("Wrote single-task fit result: ", part_path)
  fit_tbl
}

read_fit_result_parts <- function(dirs) {
  part_files <- list.files(dirs$fit_parts, pattern = "^task_[0-9]+\\.tsv$", full.names = TRUE)
  if (!length(part_files)) {
    fallback <- file.path(dirs$tables, "fit_results.tsv")
    if (file.exists(fallback)) return(read_tsv(fallback))
    return(data.frame())
  }
  parts <- lapply(sort(part_files), read_tsv)
  nms <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(x) {
    missing <- setdiff(nms, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, nms, drop = FALSE]
  })
  out <- do.call(rbind, parts)
  if ("task_order" %in% names(out)) {
    out <- out[order(as.integer(out$task_order)), , drop = FALSE]
  }
  out
}

summarize_prepared_results <- function(cfg, dirs) {
  task_tbl <- apply_path_map_df(read_tsv(file.path(dirs$tables, "task_table.tsv")), cfg$path_map)
  input_tbl <- apply_path_map_df(read_tsv(file.path(dirs$tables, "input_table.tsv")), cfg$path_map)
  repo_versions <- if (file.exists(file.path(dirs$tables, "repo_versions.tsv"))) {
    read_tsv(file.path(dirs$tables, "repo_versions.tsv"))
  } else {
    data.frame()
  }
  fit_tbl <- apply_path_map_df(read_fit_result_parts(dirs), cfg$path_map)
  if ("task_order" %in% names(fit_tbl)) {
    fit_tbl <- fit_tbl[order(as.integer(fit_tbl$task_order)), , drop = FALSE]
  }
  write_tsv(fit_tbl, file.path(dirs$tables, "fit_results.tsv"))

  missing <- data.frame()
  if (nrow(task_tbl)) {
    completed <- unique(as.integer(fit_tbl$task_order))
    missing <- task_tbl[!(as.integer(task_tbl$task_order) %in% completed), , drop = FALSE]
  }
  write_tsv(missing, file.path(dirs$tables, "missing_fit_tasks.tsv"))
  if (nrow(missing)) {
    message("Missing fit task results: ", nrow(missing), " / ", nrow(task_tbl))
  }

  write_node_accuracy <- is.null(cfg$write_node_accuracy) || isTRUE(cfg$write_node_accuracy)
  if (write_node_accuracy) {
    node_rows <- list()
    node_idx <- 0L
    for (i in seq_len(nrow(fit_tbl))) {
      nodes <- read_nodes_for_fit_row(fit_tbl[i, , drop = FALSE])
      if (nrow(nodes)) {
        node_idx <- node_idx + 1L
        node_rows[[node_idx]] <- nodes
      }
    }
    node_tbl <- if (length(node_rows)) do.call(rbind, node_rows) else data.frame()
    if (nrow(node_tbl)) {
      node_tbl$estimation_error <- node_tbl$estimated_fitness - node_tbl$true_fitness
    }
    write_tsv(node_tbl, file.path(dirs$tables, "node_accuracy.tsv"))
    summary_tbl <- summarize_accuracy(node_tbl)
  } else {
    message("Skipping node_accuracy.tsv; computing final summary tables in streaming mode.")
    node_tbl <- empty_node_accuracy_table()
    write_tsv(node_tbl, file.path(dirs$tables, "node_accuracy.tsv"))
    summary_tbl <- summarize_fit_results_streaming(fit_tbl)
  }
  write_tsv(summary_tbl, file.path(dirs$tables, "summary_by_lambda_time_minobs_method.tsv"))
  delta_tbl <- make_delta_vs_alfakR(summary_tbl)
  write_tsv(delta_tbl, file.path(dirs$tables, "alfak2_delta_vs_alfakR.tsv"))

  saveRDS(
    list(
      config = cfg,
      repo_versions = repo_versions,
      input_table = input_tbl,
      task_table = task_tbl,
      fit_results = fit_tbl,
      node_accuracy = node_tbl,
      summary = summary_tbl,
      delta = delta_tbl,
      missing_fit_tasks = missing
    ),
    file.path(dirs$root, "grf_alfak2_vs_alfakR_benchmark.rds")
  )

  message("Wrote benchmark outputs under: ", dirs$root)
  message("  ", file.path(dirs$tables, "fit_results.tsv"))
  message("  ", file.path(dirs$tables, "missing_fit_tasks.tsv"))
  message("  ", file.path(dirs$tables, "node_accuracy.tsv"))
  message("  ", file.path(dirs$tables, "summary_by_lambda_time_minobs_method.tsv"))
  message("  ", file.path(dirs$tables, "alfak2_delta_vs_alfakR.tsv"))

  invisible(list(
    input_table = input_tbl,
    task_table = task_tbl,
    fit_results = fit_tbl,
    node_accuracy = node_tbl,
    summary = summary_tbl,
    delta = delta_tbl,
    missing_fit_tasks = missing
  ))
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    return(invisible(NULL))
  }
  repo_dir <- normalizePath(arg_value(args, "repo_dir", resolve_repo_dir()), winslash = "/", mustWork = FALSE)
  cfg <- build_config(args, repo_dir)
  cfg$path_map <- parse_path_map(args)
  dirs <- build_dirs(cfg$output_dir)
  mode <- tolower(as.character(arg_value(args, "mode", "all")))
  valid_modes <- c("all", "prepare", "fit-task", "fit_task", "summarize")
  if (!mode %in% valid_modes) {
    stop("Unsupported --mode=", mode, ". Use one of: all, prepare, fit-task, summarize.", call. = FALSE)
  }
  if (identical(mode, "fit_task")) mode <- "fit-task"
  if (mode %in% c("fit-task", "summarize")) {
    cfg <- load_prepared_config(cfg, dirs)
    cfg <- apply_path_map_list(cfg, parse_path_map(args))
    cfg <- apply_runtime_overrides(cfg, args, repo_dir)
    cfg$path_map <- parse_path_map(args)
    dirs <- build_dirs(cfg$output_dir)
  }
  recompile_dll <- if (!is.null(args$recompile_dll)) {
    arg_logical(args, "recompile_dll", FALSE)
  } else {
    mode %in% c("all", "prepare")
  }

  repo_versions <- rbind(
    repo_state(cfg$alfakR_repo, "alfakR"),
    repo_state(cfg$alfak2_repo, "alfak2")
  )
  if (!identical(mode, "fit-task")) {
    write_tsv(repo_versions, file.path(dirs$tables, "repo_versions.tsv"))
  }

  message("Loading current source trees with pkgload::load_all().")
  message("  alfakR: ", cfg$alfakR_repo)
  message("  alfak2: ", cfg$alfak2_repo)
  load_current_repos(cfg$alfakR_repo, cfg$alfak2_repo, recompile_dll = recompile_dll)

  if (any(repo_versions$dirty) && !isTRUE(cfg$reuse_dirty_cache)) {
    message("At least one repo is dirty; fit caches will not be reused unless --reuse-dirty-cache=true.")
  }

  if (identical(mode, "prepare")) {
    return(invisible(prepare_benchmark_inputs(cfg, dirs, repo_versions)))
  }
  if (identical(mode, "fit-task")) {
    slurm_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", NA_character_)))
    task_index <- arg_integer(args, "task_index", if (is.finite(slurm_idx)) slurm_idx else 1L)
    return(invisible(run_one_prepared_task(cfg, dirs, repo_versions, task_index)))
  }
  if (identical(mode, "summarize")) {
    return(invisible(summarize_prepared_results(cfg, dirs)))
  }

  input_rows <- list()
  fit_tasks <- list()
  input_idx <- 0L
  task_idx <- 0L
  grf_lookup <- list()
  time_axis_label <- paste0("tmax_", path_token(cfg$time_max), "_pint_", path_token(cfg$passage_interval))
  sim_pm_label <- pm_to_label(cfg$pm)

  for (sim_idx in seq_len(cfg$n_sim)) {
    for (lambda_idx in seq_along(cfg$lambdas)) {
      lambda <- cfg$lambdas[[lambda_idx]]
      lambda_label <- format_grf_label(lambda)
      abm_seed <- cfg$seed + sim_idx * 10000L + lambda_idx * 100L
      grf_key <- paste(sim_idx, lambda_label, paste0("simpm_", sim_pm_label), time_axis_label, sep = "__")
      grf_path <- file.path(dirs$cache, paste0("grf_sim_", grf_key, ".rds"))
      grf_sim <- if (!isTRUE(cfg$force_sim) &&
                     file.exists(grf_path) &&
                     (isTRUE(cfg$reuse_dirty_cache) || !isTRUE(repo_versions$dirty[repo_versions$package == "alfakR"]))) {
        readRDS(grf_path)
      } else {
        message("Simulating GRF ABM: sim=", sim_idx, " lambda=", lambda)
        out <- simulate_nn_prior_grf_abm(
          seed = abm_seed,
          lambda = lambda,
          p = cfg$pm,
          k_dim = cfg$k_dim,
          n_centroids = cfg$n_centroids,
          time_max = cfg$time_max,
          passage_interval = cfg$passage_interval,
          abm_pop_size = cfg$abm_pop_size,
          abm_delta_t = cfg$abm_delta_t,
          abm_max_pop = cfg$abm_max_pop,
          abm_culling_survival = cfg$abm_culling_survival
        )
        saveRDS(out, grf_path)
        out
      }
      grf_lookup[[grf_key]] <- grf_sim

      for (time_start in cfg$time_starts) {
        for (time_gap in cfg$time_gaps) {
          patient_id <- paste0("grf_", sim_idx, "_lambda_", lambda_label, "_simpm_", sim_pm_label, "_", time_axis_label, "_start_", path_token(time_start), "_gap_", path_token(time_gap))
          input_rds <- file.path(dirs$cache, paste0("input_", patient_id, ".rds"))
          input_csv <- file.path(dirs$cache, paste0("input_", patient_id, ".csv"))
          input_seed <- abm_seed + as.integer(round(time_start * 10)) + as.integer(round(time_gap * 1000))
          yi <- if (!file.exists(input_rds) || isTRUE(cfg$force_sim)) {
            out <- build_two_timepoint_yi_from_abm(
              sim_wide = grf_sim$sim_wide,
              time_start = time_start,
              time_gap = time_gap,
              passage_interval = cfg$passage_interval,
              sample_depth = cfg$sample_depth,
              seed = input_seed
            )
            saveRDS(out, input_rds)
            utils::write.csv(data.frame(karyotype = rownames(out$x), out$x, check.names = FALSE), input_csv, row.names = FALSE)
            out
          } else {
            readRDS(input_rds)
          }
          input_md5 <- unname(tools::md5sum(input_rds))
          input_summary <- summarize_input_rows(yi, cfg$minobs, drop_diploid = cfg$drop_diploid)
          for (rr in seq_len(nrow(input_summary))) {
            input_idx <- input_idx + 1L
            input_rows[[input_idx]] <- data.frame(
              simulation_id = sim_idx,
              lambda = lambda,
              lambda_label = lambda_label,
              time_start = time_start,
              time_gap = time_gap,
              time_delta = as.numeric(yi$metadata$time_delta),
              sim_pm = cfg$pm,
              patient_id = patient_id,
              grf_key = grf_key,
              grf_rds = grf_path,
              input_rds = input_rds,
              input_csv = input_csv,
              input_md5 = input_md5,
              input_summary[rr, , drop = FALSE],
              stringsAsFactors = FALSE
            )
          }

          for (minobs in cfg$minobs) {
            for (fit_beta in cfg$beta_levels) {
              fit_beta_label <- pm_to_label(fit_beta)
              for (policy in cfg$alfak2_input_policies) {
                task_idx <- task_idx + 1L
                method <- paste0("alfak2_", cfg$alfak2_input_depth, "_", policy)
                fit_tasks[[task_idx]] <- data.frame(
                  engine = "alfak2",
                  method = method,
                  input_policy = policy,
                  nn_prior = NA_character_,
                  simulation_id = sim_idx,
                  lambda = lambda,
                  lambda_label = lambda_label,
                  time_start = time_start,
                  time_gap = time_gap,
                  time_delta = as.numeric(yi$metadata$time_delta),
                  sim_pm = cfg$pm,
                  pm = fit_beta,
                  fit_beta_label = fit_beta_label,
                  patient_id = patient_id,
                  grf_key = grf_key,
                  grf_rds = grf_path,
                  input_rds = input_rds,
                  input_md5 = input_md5,
                  minobs = as.integer(minobs),
                  benchmark_seed = as.integer(abm_seed + time_gap * 1000L + minobs * 100L + 900L + match(policy, cfg$alfak2_input_policies)),
                  outdir = file.path(dirs$fits, "alfak2", paste0("lambda_", lambda_label), paste0("gap_", path_token(time_gap)), paste0("beta_", fit_beta_label), paste0("MINOBS_", minobs), method, patient_id),
                  stringsAsFactors = FALSE
                )
              }
              for (method in cfg$methods) {
                task_idx <- task_idx + 1L
                fit_tasks[[task_idx]] <- data.frame(
                  engine = "alfakR",
                  method = paste0("alfakR_", method),
                  input_policy = "alfakR_minobs_internal",
                  nn_prior = method,
                  simulation_id = sim_idx,
                  lambda = lambda,
                  lambda_label = lambda_label,
                  time_start = time_start,
                  time_gap = time_gap,
                  time_delta = as.numeric(yi$metadata$time_delta),
                  sim_pm = cfg$pm,
                  pm = fit_beta,
                  fit_beta_label = fit_beta_label,
                  patient_id = patient_id,
                  grf_key = grf_key,
                  grf_rds = grf_path,
                  input_rds = input_rds,
                  input_md5 = input_md5,
                  minobs = as.integer(minobs),
                  benchmark_seed = as.integer(abm_seed + time_gap * 1000L + minobs * 100L + match(method, cfg$methods)),
                  outdir = file.path(dirs$fits, "alfakR", paste0("lambda_", lambda_label), paste0("gap_", path_token(time_gap)), paste0("beta_", fit_beta_label), paste0("MINOBS_", minobs), paste0("nn_prior_", method), patient_id),
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }
      }
    }
  }

  input_tbl <- do.call(rbind, input_rows)
  task_tbl <- do.call(rbind, fit_tasks)
  task_tbl <- task_tbl[c(which(task_tbl$engine == "alfak2"), which(task_tbl$engine == "alfakR")), , drop = FALSE]
  task_tbl$task_order <- seq_len(nrow(task_tbl))
  validate_task_shared_sources(task_tbl)
  write_tsv(input_tbl, file.path(dirs$tables, "input_table.tsv"))
  write_tsv(task_tbl, file.path(dirs$tables, "task_table.tsv"))
  saveRDS(cfg, file.path(dirs$root, "benchmark_config.rds"))

  message("Prepared ", nrow(task_tbl), " fit tasks.")
  fit_rows <- c(
    run_fit_task_table(task_tbl[task_tbl$engine == "alfak2", , drop = FALSE], cfg, repo_versions, "alfak2 phase"),
    run_fit_task_table(task_tbl[task_tbl$engine == "alfakR", , drop = FALSE], cfg, repo_versions, "alfakR phase")
  )
  fit_tbl <- list_to_data_frame(fit_rows)
  write_tsv(fit_tbl, file.path(dirs$tables, "fit_results.tsv"))

  node_rows <- list()
  node_idx <- 0L
  for (i in seq_len(nrow(fit_tbl))) {
    fr <- fit_tbl[i, , drop = FALSE]
    if (!identical(as.character(fr$status[[1]]), "ok")) next
    grf_key <- as.character(fr$grf_key[[1]])
    grf_sim <- grf_lookup[[grf_key]]
    if (is.null(grf_sim)) {
      grf_path <- as.character(fr$grf_rds[[1]])
      if (!nzchar(grf_path) || is.na(grf_path)) {
        grf_path <- file.path(dirs$cache, paste0("grf_sim_", grf_key, ".rds"))
      }
      if (file.exists(grf_path)) grf_sim <- readRDS(grf_path)
    }
    if (is.null(grf_sim)) next
    if (identical(as.character(fr$engine[[1]]), "alfakR")) {
      landscape <- read_rds_if_exists(as.character(fr$landscape_path[[1]]))
      nodes <- coerce_alfakR_landscape_nodes(
        landscape = landscape,
        method = as.character(fr$method[[1]]),
        nn_prior = as.character(fr$nn_prior[[1]]),
        fit_row = fr,
        centroids = grf_sim$centroids,
        lambda = as.numeric(fr$lambda[[1]])
      )
    } else {
      fit <- read_rds_if_exists(as.character(fr$fit_path[[1]]))
      if (is.null(fit)) next
      nodes <- coerce_alfak2_nodes(
        fit = fit,
        method = as.character(fr$method[[1]]),
        input_policy = as.character(fr$input_policy[[1]]),
        fit_row = fr,
        centroids = grf_sim$centroids,
        lambda = as.numeric(fr$lambda[[1]]),
        use_legacy_scale = TRUE
      )
    }
    if (nrow(nodes)) {
      node_idx <- node_idx + 1L
      node_rows[[node_idx]] <- nodes
    }
  }
  node_tbl <- if (length(node_rows)) do.call(rbind, node_rows) else data.frame()
  if (nrow(node_tbl)) {
    node_tbl$estimation_error <- node_tbl$estimated_fitness - node_tbl$true_fitness
  }
  write_tsv(node_tbl, file.path(dirs$tables, "node_accuracy.tsv"))

  summary_tbl <- summarize_accuracy(node_tbl)
  write_tsv(summary_tbl, file.path(dirs$tables, "summary_by_lambda_time_minobs_method.tsv"))
  delta_tbl <- make_delta_vs_alfakR(summary_tbl)
  write_tsv(delta_tbl, file.path(dirs$tables, "alfak2_delta_vs_alfakR.tsv"))

  saveRDS(
    list(
      config = cfg,
      repo_versions = repo_versions,
      input_table = input_tbl,
      task_table = task_tbl,
      fit_results = fit_tbl,
      node_accuracy = node_tbl,
      summary = summary_tbl,
      delta = delta_tbl
    ),
    file.path(dirs$root, "grf_alfak2_vs_alfakR_benchmark.rds")
  )

  message("Wrote benchmark outputs under: ", dirs$root)
  message("  ", file.path(dirs$tables, "repo_versions.tsv"))
  message("  ", file.path(dirs$tables, "input_table.tsv"))
  message("  ", file.path(dirs$tables, "task_table.tsv"))
  message("  ", file.path(dirs$tables, "fit_results.tsv"))
  message("  ", file.path(dirs$tables, "node_accuracy.tsv"))
  message("  ", file.path(dirs$tables, "summary_by_lambda_time_minobs_method.tsv"))
  message("  ", file.path(dirs$tables, "alfak2_delta_vs_alfakR.tsv"))

  invisible(list(
    input_table = input_tbl,
    task_table = task_tbl,
    fit_results = fit_tbl,
    node_accuracy = node_tbl,
    summary = summary_tbl,
    delta = delta_tbl
  ))
}

if (sys.nframe() == 0L) {
  main()
}
