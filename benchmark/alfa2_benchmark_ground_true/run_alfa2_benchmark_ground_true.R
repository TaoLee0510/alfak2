#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Run the ALFA-K-original ground-truth benchmark for alfak2 and alfakR.\n\n",
    "Usage:\n",
    "  Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=prepare\n",
    "  Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=ground-truth --ground-truth-index=1\n",
    "  Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=fit-task --task-index=1\n",
    "  Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=summarize\n\n",
    "Modes:\n",
    "  prepare       Write ground_truth_index.tsv, method_index.tsv, and run_index.tsv.\n",
    "  ground-truth  Generate one or all ground-truth RDS files using the ALFA-K original ABM process.\n",
    "  fit-task      Run one run_index.tsv row. Missing ground truth is generated on demand.\n",
    "  summarize     Combine completed run_cache/*.rds metrics into summary tables.\n",
    "  all           Prepare, generate all ground truth, run all tasks sequentially, and summarize.\n\n",
    "Core options:\n",
    "  --output-dir=benchmark/results/alfa2_benchmark_ground_true\n",
    "  --alfak2-repo=/share/lab_crd/lab_crd/taoli/Project/alfak2\n",
    "  --alfakR-repo=/share/lab_crd/lab_crd/taoli/Project/alfakR\n",
    "  --sample-depths=1000,200\n",
    "  --wavelengths=0.2,0.4,0.8,1.6\n",
    "  --ground-truth-times=0,180\n",
    "  --passage-times=0,180\n",
    "  --ground-truth-reps=1:5\n",
    "  --fit-repeats=1:5\n",
    "  --soft-minobs=5,10,20\n",
    "  --ntp=2\n",
    "  --alfakR-dt=1\n",
    "  --nboot=45\n",
    "  --forward-prediction-reps=5\n",
    "  --force=false\n",
    sep = ""
  )
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

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

arg_integer <- function(args, name, default) {
  value <- suppressWarnings(as.integer(arg_value(args, name, default)))
  if (!is.finite(value)) default else value
}

arg_numeric <- function(args, name, default) {
  value <- suppressWarnings(as.numeric(arg_value(args, name, default)))
  if (!is.finite(value)) default else value
}

parse_numeric_vec <- function(value, default) {
  if (is.null(value)) return(default)
  parts <- trimws(strsplit(as.character(value), ",", fixed = TRUE)[[1L]])
  out <- suppressWarnings(as.numeric(parts))
  out <- out[is.finite(out)]
  if (!length(out)) default else out
}

parse_integer_vec <- function(value, default) {
  if (is.null(value)) return(default)
  value <- trimws(as.character(value))
  if (grepl("^[0-9]+:[0-9]+$", value)) {
    parts <- as.integer(strsplit(value, ":", fixed = TRUE)[[1L]])
    return(seq(parts[[1L]], parts[[2L]]))
  }
  out <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  out <- out[is.finite(out)]
  if (!length(out)) default else out
}

arg_numeric_vec <- function(args, name, default) parse_numeric_vec(arg_value(args, name, NULL), default)
arg_integer_vec <- function(args, name, default) parse_integer_vec(arg_value(args, name, NULL), default)

resolve_repo_dir <- function(start = getwd()) {
  candidates <- unique(normalizePath(
    file.path(start, c(".", "..", "../..", "../../..")),
    winslash = "/",
    mustWork = FALSE
  ))
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) && dir.exists(file.path(cand, "benchmark"))) {
      return(cand)
    }
  }
  stop("Could not locate the alfak2 repository root.", call. = FALSE)
}

normalize_output_dir <- function(repo_dir, output_dir) {
  if (grepl("^/", output_dir)) {
    normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(repo_dir, output_dir), winslash = "/", mustWork = FALSE)
  }
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

rbind_fill <- function(xs) {
  xs <- Filter(function(x) !is.null(x) && is.data.frame(x) && nrow(x), xs)
  if (!length(xs)) return(data.frame())
  cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(cols, names(x))
    for (m in missing) x[[m]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

label_number <- function(x) gsub("\\.", "p", as.character(x))
pad2 <- function(x) sprintf("%02d", as.integer(x))

safe_read_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) NULL)
}

numeric_vec_equal <- function(x, y, tol = 1e-8) {
  x <- suppressWarnings(as.numeric(x))
  y <- suppressWarnings(as.numeric(y))
  length(x) == length(y) &&
    all(is.finite(x)) &&
    all(is.finite(y)) &&
    all(abs(x - y) <= tol)
}

required_benchmark_times <- function() c(0, 180)

validate_exact_benchmark_times <- function(x, label) {
  required <- required_benchmark_times()
  if (!numeric_vec_equal(x, required)) {
    stop(
      label, " must be exactly ", paste(required, collapse = ","),
      " for alfa2_benchmark_ground_true.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

default_alfak2_methods <- function() {
  c(
    "graph_gaussian_baseline",
    "edge_effect_empirical_bayes",
    "edge_effect_interaction_path_ensemble",
    "kronecker_or_graph_trend_filtering",
    "local_NNGP_or_GPnn",
    "delta_tree_ensemble",
    "tabpfn_nearfield_feature_model",
    "truncated_nearfield_gmrf",
    "local_polynomial_stencil"
  )
}

make_dirs <- function(output_dir) {
  dirs <- list(
    root = output_dir,
    ground_truth = file.path(output_dir, "ground_truth"),
    runs = file.path(output_dir, "runs"),
    cache = file.path(output_dir, "run_cache"),
    tables = file.path(output_dir, "tables"),
    logs = file.path(output_dir, "logs")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

build_ground_truth_index <- function(output_dir,
                                     sample_depths = c(1000L, 200L),
                                     wavelengths = c(0.2, 0.4, 0.8, 1.6),
                                     ground_truth_reps = 1:5,
                                     seed_base = 720000L,
                                     pmis = 5e-05) {
  grid <- expand.grid(
    sample_depth_index = seq_along(sample_depths),
    wavelength_index = seq_along(wavelengths),
    ground_truth_repeat = ground_truth_reps,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$sample_depth <- as.integer(sample_depths[grid$sample_depth_index])
  grid$wavelength <- as.numeric(wavelengths[grid$wavelength_index])
  grid$pmis <- pmis
  grid$ground_truth_id <- sprintf(
    "depth%s_w%s_gt%s",
    grid$sample_depth,
    label_number(grid$wavelength),
    pad2(grid$ground_truth_repeat)
  )
  grid$ground_truth_rds <- file.path(
    output_dir,
    "ground_truth",
    paste0("depth_", grid$sample_depth),
    paste0(
      "w_", label_number(grid$wavelength),
      "_m_", label_number(grid$pmis),
      "_rep_", pad2(grid$ground_truth_repeat),
      ".Rds"
    )
  )
  grid$landscape_seed <- as.integer(seed_base +
    100000L * as.integer(grid$sample_depth_index) +
    1000L * as.integer(grid$wavelength_index) +
    as.integer(grid$ground_truth_repeat))
  grid <- grid[order(grid$sample_depth, grid$wavelength, grid$ground_truth_repeat), ]
  grid$ground_truth_index <- seq_len(nrow(grid))
  grid[, c(
    "ground_truth_index", "ground_truth_id", "sample_depth", "sample_depth_index",
    "wavelength", "wavelength_index", "ground_truth_repeat", "pmis",
    "landscape_seed", "ground_truth_rds"
  )]
}

build_method_index <- function(soft_minobs = c(5L, 10L, 20L),
                               alfak2_methods = default_alfak2_methods(),
                               alfakR_priors = c("none", "empirical", "empirical_censored", "empirical_censored_weighted")) {
  alfak2_inputs <- data.frame(
    input_mode = c("full", rep("soft_minobs", length(soft_minobs))),
    soft_minobs = c(NA_integer_, as.integer(soft_minobs)),
    stringsAsFactors = FALSE
  )
  alfak2 <- merge(
    alfak2_inputs,
    data.frame(extrapolation_method = alfak2_methods, stringsAsFactors = FALSE),
    all = TRUE
  )
  alfak2$package <- "alfak2"
  alfak2$minobs <- NA_integer_
  alfak2$NN_prior <- NA_character_
  alfak2$NN_prior_display <- NA_character_

  alfakR <- expand.grid(
    minobs = as.integer(soft_minobs),
    NN_prior = alfakR_priors,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  alfakR$package <- "alfakR"
  alfakR$input_mode <- NA_character_
  alfakR$soft_minobs <- NA_integer_
  alfakR$extrapolation_method <- NA_character_
  alfakR$NN_prior_display <- ifelse(alfakR$NN_prior == "none", "None", alfakR$NN_prior)

  cols <- c("package", "input_mode", "soft_minobs", "extrapolation_method", "minobs", "NN_prior", "NN_prior_display")
  out <- rbind(alfak2[, cols], alfakR[, cols])
  out$method_pair_id <- ifelse(
    out$package == "alfak2",
    ifelse(
      out$input_mode == "full",
      paste("alfak2", "full", out$extrapolation_method, sep = ":"),
      paste("alfak2", paste0("soft_minobs", out$soft_minobs), out$extrapolation_method, sep = ":")
    ),
    paste("alfakR", paste0("minobs", out$minobs), out$NN_prior_display, sep = ":")
  )
  out$method_index <- seq_len(nrow(out))
  out[, c("method_index", "method_pair_id", cols)]
}

build_run_index <- function(ground_truth_index, method_index, fit_repeats = 1:5, seed_base = 910000L) {
  repeat_grid <- data.frame(fit_repeat = as.integer(fit_repeats), stringsAsFactors = FALSE)
  gt_rep <- merge(ground_truth_index, repeat_grid, all = TRUE)
  out <- merge(gt_rep, method_index, all = TRUE)
  out <- out[order(
    out$sample_depth, out$wavelength, out$ground_truth_repeat, out$fit_repeat,
    out$package, out$input_mode, out$soft_minobs, out$extrapolation_method, out$minobs, out$NN_prior
  ), ]
  out$task_id <- seq_len(nrow(out))
  out$run_id <- sprintf("alfa2_gt_%06d", out$task_id)
  out$fit_seed <- as.integer(seed_base +
    100000L * as.integer(out$sample_depth_index) +
    10000L * as.integer(out$wavelength_index) +
    100L * as.integer(out$ground_truth_repeat) +
    10L * as.integer(out$fit_repeat) +
    as.integer(out$method_index))
  out[, c("task_id", "run_id", setdiff(names(out), c("task_id", "run_id")))]
}

write_parameter_pair_counts <- function(method_index, run_index, tables_dir) {
  method_counts <- aggregate(
    method_pair_id ~ package,
    method_index,
    function(x) length(unique(x))
  )
  names(method_counts)[names(method_counts) == "method_pair_id"] <- "n_method_parameter_pairs"

  depth_counts <- aggregate(
    method_pair_id ~ sample_depth + package,
    unique(run_index[, c("sample_depth", "package", "method_pair_id")]),
    function(x) length(unique(x))
  )
  names(depth_counts)[names(depth_counts) == "method_pair_id"] <- "n_method_parameter_pairs_per_sample_depth"
  write_tsv(method_counts, file.path(tables_dir, "method_parameter_pair_counts.tsv"))
  write_tsv(depth_counts, file.path(tables_dir, "method_parameter_pair_counts_by_depth.tsv"))
}

prepare_indices <- function(cfg, dirs) {
  gt <- build_ground_truth_index(
    output_dir = cfg$output_dir,
    sample_depths = cfg$sample_depths,
    wavelengths = cfg$wavelengths,
    ground_truth_reps = cfg$ground_truth_reps,
    seed_base = cfg$ground_truth_seed_base,
    pmis = cfg$pmis
  )
  methods <- build_method_index(
    soft_minobs = cfg$soft_minobs,
    alfak2_methods = default_alfak2_methods(),
    alfakR_priors = cfg$alfakR_priors
  )
  runs <- build_run_index(
    ground_truth_index = gt,
    method_index = methods,
    fit_repeats = cfg$fit_repeats,
    seed_base = cfg$fit_seed_base
  )
  write_tsv(gt, file.path(dirs$tables, "ground_truth_index.tsv"))
  write_tsv(methods, file.path(dirs$tables, "method_index.tsv"))
  write_tsv(runs, file.path(dirs$tables, "run_index.tsv"))
  write_parameter_pair_counts(methods, runs, dirs$tables)
  write_tsv(data.frame(
    key = c(
      "sample_depths", "wavelengths", "ground_truth_reps", "fit_repeats",
      "ground_truth_times", "passage_times", "soft_minobs", "alfakR_priors",
      "ntp", "alfakR_dt", "nboot", "forward_prediction_reps", "pmis",
      "n0", "nb", "alfak2_local_shell_depth", "alfak2_global_extra_shell",
      "alfak2_max_nodes"
    ),
    value = c(
      paste(cfg$sample_depths, collapse = ","),
      paste(cfg$wavelengths, collapse = ","),
      paste(cfg$ground_truth_reps, collapse = ","),
      paste(cfg$fit_repeats, collapse = ","),
      paste(cfg$ground_truth_times, collapse = ","),
      paste(cfg$passage_times, collapse = ","),
      paste(cfg$soft_minobs, collapse = ","),
      paste(cfg$alfakR_priors, collapse = ","),
      cfg$ntp,
      cfg$alfakR_dt,
      cfg$nboot,
      cfg$forward_prediction_reps,
      cfg$pmis,
      cfg$n0,
      cfg$nb,
      cfg$alfak2_local_shell_depth,
      cfg$alfak2_global_extra_shell,
      cfg$alfak2_max_nodes
    ),
    stringsAsFactors = FALSE
  ), file.path(dirs$tables, "benchmark_config.tsv"))
  list(ground_truth_index = gt, method_index = methods, run_index = runs)
}

load_indices <- function(dirs) {
  paths <- file.path(dirs$tables, c("ground_truth_index.tsv", "method_index.tsv", "run_index.tsv"))
  if (!all(file.exists(paths))) {
    stop("Missing index TSVs. Run --mode=prepare first.", call. = FALSE)
  }
  list(
    ground_truth_index = read_tsv(paths[[1L]]),
    method_index = read_tsv(paths[[2L]]),
    run_index = read_tsv(paths[[3L]])
  )
}

validate_source_repo <- function(path, package) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  desc <- file.path(path, "DESCRIPTION")
  if (!file.exists(desc)) {
    stop("Missing DESCRIPTION for ", package, " source repo: ", path, call. = FALSE)
  }
  dcf <- read.dcf(desc)
  pkg <- as.character(dcf[1L, "Package"])
  if (!identical(pkg, package)) {
    stop("Expected package `", package, "` at ", path, " but DESCRIPTION says `", pkg, "`.", call. = FALSE)
  }
  path
}

loaded_namespace_path <- function(package) {
  normalizePath(getNamespaceInfo(package, "path"), winslash = "/", mustWork = TRUE)
}

load_repositories <- function(repo_dir, alfakR_repo) {
  repo_dir <- validate_source_repo(repo_dir, "alfak2")
  alfakR_repo <- validate_source_repo(alfakR_repo, "alfakR")
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop(
      "Package `pkgload` is required so the benchmark loads alfak2/alfakR from source repos, ",
      "not from the active R module library.",
      call. = FALSE
    )
  }

  pkgload::load_all(repo_dir, quiet = TRUE)
  pkgload::load_all(alfakR_repo, quiet = TRUE)

  loaded_alfak2 <- loaded_namespace_path("alfak2")
  loaded_alfakR <- loaded_namespace_path("alfakR")
  if (!identical(loaded_alfak2, repo_dir)) {
    stop("alfak2 was loaded from ", loaded_alfak2, " instead of ", repo_dir, call. = FALSE)
  }
  if (!identical(loaded_alfakR, alfakR_repo)) {
    stop("alfakR was loaded from ", loaded_alfakR, " instead of ", alfakR_repo, call. = FALSE)
  }

  message("Loaded alfak2 source repo: ", loaded_alfak2)
  message("Loaded alfakR source repo: ", loaded_alfakR)
  invisible(TRUE)
}

ground_truth_has_times <- function(path, expected_times) {
  xi <- safe_read_rds(path)
  if (is.null(xi) || is.null(xi$abm_output) || is.null(xi$abm_output$x)) return(FALSE)
  actual_times <- suppressWarnings(as.numeric(colnames(xi$abm_output$x)))
  numeric_vec_equal(actual_times, expected_times)
}

generate_ground_truth <- function(row, ground_truth_times = c(0, 180), force = FALSE) {
  out_path <- as.character(row$ground_truth_rds[[1L]])
  ground_truth_times <- as.numeric(ground_truth_times)
  validate_exact_benchmark_times(ground_truth_times, "ground_truth_times")
  if (!force && ground_truth_has_times(out_path, ground_truth_times)) return(out_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(out_path, ".lock")
  lock_start <- Sys.time()
  while (!dir.create(lock_dir, showWarnings = FALSE)) {
    if (!force && ground_truth_has_times(out_path, ground_truth_times)) return(out_path)
    waited <- as.numeric(difftime(Sys.time(), lock_start, units = "mins"))
    if (is.finite(waited) && waited > 180) {
      stop("Timed out waiting for ground-truth lock: ", lock_dir, call. = FALSE)
    }
    Sys.sleep(10)
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  if (!force && ground_truth_has_times(out_path, ground_truth_times)) return(out_path)

  set.seed(as.integer(row$landscape_seed[[1L]]))

  gen_randscape <- function(founder, Nwaves, scalef = NULL, wavelength = 0.8) {
    if (is.null(scalef)) scalef <- 1 / (pi * sqrt(Nwaves))
    f0 <- 0
    while (f0 <= 0.4) {
      pk <- lapply(seq_len(Nwaves), function(i) {
        sample((-10):20, length(founder), replace = TRUE)
      })

      d <- sapply(pk, function(ci) {
        sqrt(sum((founder - ci)^2))
      })
      f0 <- sum(sin(d / wavelength) * scalef)
    }
    do.call(rbind, pk)
  }

  resample_sim <- function(sim, n_samples, output_times) {
    sim_times <- suppressWarnings(as.numeric(sim$time))
    if (!length(sim_times) || any(!is.finite(sim_times))) {
      stop("ABM output has non-finite simulation times.", call. = FALSE)
    }
    output_times <- as.numeric(output_times)
    time_idx <- vapply(output_times, function(target) {
      exact <- which(abs(sim_times - target) <= 1e-8)
      if (length(exact)) exact[[1L]] else which.min(abs(sim_times - target))
    }, integer(1L))
    if (length(unique(time_idx)) != length(output_times)) {
      stop("Could not map requested benchmark times to distinct ABM snapshots.", call. = FALSE)
    }
    mat <- t(as.matrix(sim[time_idx, -1, drop = FALSE]))
    colnames(mat) <- as.character(output_times)
    mat <- apply(mat, 2, function(p) rmultinom(1, n_samples, prob = p / sum(p)))
    keep <- rowSums(mat) > 0
    rownames(mat) <- colnames(sim)[-1]
    mat <- mat[keep, , drop = FALSE]
    mat <- mat[order(rowSums(mat), decreasing = TRUE), ]
    attr(mat, "actual_abm_times") <- sim_times[time_idx]
    mat
  }

  founder <- rep(2, 22)
  Nwaves <- 10
  wavelength <- as.numeric(row$wavelength[[1L]])
  l <- gen_randscape(founder, Nwaves, wavelength = wavelength)
  times <- ground_truth_times
  x0 <- c(1)
  names(x0) <- paste(founder, collapse = ".")
  pmis <- as.numeric(row$pmis[[1L]])

  sim <- alfakR::run_abm_simulation_grf(
    centroids = l,
    lambda = wavelength,
    p = pmis,
    times = times,
    x0 = x0,
    abm_pop_size = 5e4,
    abm_max_pop = 2e6,
    abm_delta_t = 0.1,
    abm_culling_survival = 0.01,
    abm_record_interval = -1,
    abm_seed = 42,
    normalize_freq = FALSE
  )

  sim_rs <- resample_sim(sim, as.integer(row$sample_depth[[1L]]), output_times = times)
  yi <- list(x = data.frame(sim_rs, check.names = FALSE), dt = 1)
  out <- list(
    abm_output = yi,
    true_landscape = l,
    ground_truth_times = times,
    actual_abm_times = as.numeric(attr(sim_rs, "actual_abm_times"))
  )
  tmp_path <- paste0(out_path, ".tmp.", Sys.getpid(), ".Rds")
  saveRDS(out, tmp_path)
  if (!file.rename(tmp_path, out_path)) {
    unlink(tmp_path, force = TRUE)
    stop("Failed to move temporary ground truth into place: ", out_path, call. = FALSE)
  }
  out_path
}

select_passage_counts <- function(xi, passage_times = c(0, 180), ntp = length(passage_times)) {
  yi <- xi$abm_output
  pass_times_all <- suppressWarnings(as.numeric(colnames(yi$x)))
  if (!length(pass_times_all) || any(!is.finite(pass_times_all))) {
    stop("Ground-truth count matrix has non-numeric timepoint columns.", call. = FALSE)
  }
  pass_times <- as.numeric(passage_times)
  validate_exact_benchmark_times(pass_times, "passage_times")
  if (length(pass_times) != as.integer(ntp)) {
    stop(
      "Requested passage_times length is ", length(pass_times),
      " but ntp is ", as.integer(ntp), ".",
      call. = FALSE
    )
  }
  keep <- match(pass_times, pass_times_all)
  if (anyNA(keep)) {
    stop(
      "Could not select required passage times ",
      paste(pass_times, collapse = ","),
      " from ground truth with times ",
      paste(pass_times_all, collapse = ","),
      ".",
      call. = FALSE
    )
  }
  counts <- as.matrix(yi$x[, keep, drop = FALSE])
  colnames(counts) <- as.character(pass_times)
  storage.mode(counts) <- "integer"
  list(
    counts = counts,
    passage_times = pass_times,
    dt = if (length(pass_times) >= 2L) diff(range(pass_times)) else yi$dt %||% 1
  )
}

cache_matches_current_config <- function(x, cfg) {
  if (is.null(x) || is.null(x$benchmark_config)) return(FALSE)
  numeric_vec_equal(x$benchmark_config$ground_truth_times, cfg$ground_truth_times) &&
    numeric_vec_equal(x$benchmark_config$passage_times, cfg$passage_times) &&
    numeric_vec_equal(x$benchmark_config$alfakR_dt, cfg$alfakR_dt) &&
    numeric_vec_equal(x$benchmark_config$forward_prediction_reps, cfg$forward_prediction_reps)
}

supported_second_layer_shells <- function(candidates) {
  choices <- tryCatch(
    eval(formals(alfak2:::second_layer_metric_values)$shell),
    error = function(e) character()
  )
  if (!length(choices)) return(candidates)
  intersect(candidates, choices)
}

validate_benchmark_time_config <- function(cfg) {
  validate_exact_benchmark_times(cfg$ground_truth_times, "ground_truth_times")
  validate_exact_benchmark_times(cfg$passage_times, "passage_times")
  if (!identical(as.integer(cfg$ntp), 2L)) {
    stop("ntp must be exactly 2 for alfa2_benchmark_ground_true.", call. = FALSE)
  }
  if (!numeric_vec_equal(cfg$alfakR_dt, 1)) {
    stop("alfakR_dt must be exactly 1 for alfa2_benchmark_ground_true.", call. = FALSE)
  }
  forward_prediction_reps <- as.integer(cfg$forward_prediction_reps)
  if (!is.finite(forward_prediction_reps) || forward_prediction_reps < 0L) {
    stop("forward_prediction_reps must be non-negative.", call. = FALSE)
  }
  invisible(TRUE)
}

parse_karyotype_matrix <- function(labels) {
  mats <- strsplit(as.character(labels), ".", fixed = TRUE)
  n <- length(mats[[1L]])
  out <- do.call(rbind, lapply(mats, function(x) as.numeric(x)))
  if (ncol(out) != n) stop("Inconsistent karyotype dimensions.", call. = FALSE)
  out
}

get_true_fitness <- function(k, wavelength, true_landscape) {
  Nwaves <- nrow(true_landscape)
  scalef <- 1 / (pi * sqrt(Nwaves))
  d <- apply(true_landscape, 1, function(ci) {
    sqrt(sum((k - ci)^2))
  })
  sum(sin(d / wavelength) * scalef)
}

true_fitness_for_labels <- function(labels, wavelength, true_landscape) {
  k <- parse_karyotype_matrix(labels)
  as.numeric(apply(k, 1, get_true_fitness, wavelength = wavelength, true_landscape = true_landscape))
}

prepare_alfak2_counts <- function(counts, input_mode, soft_minobs = NA_integer_) {
  counts <- as.matrix(counts)
  storage.mode(counts) <- "integer"
  if (identical(input_mode, "full")) return(counts)
  if (!identical(input_mode, "soft_minobs")) {
    stop("Unsupported alfak2 input mode: ", input_mode, call. = FALSE)
  }
  minobs <- as.integer(soft_minobs)
  if (!is.finite(minobs) || minobs <= 0L) {
    stop("soft_minobs must be a positive integer for soft_minobs input.", call. = FALSE)
  }
  row_totals <- rowSums(counts, na.rm = TRUE)
  weights <- pmin(1, pmax(0, row_totals) / as.numeric(minobs))
  weights[!is.finite(weights)] <- 0
  weights <- cbind(t0 = weights, t1 = weights)
  rownames(weights) <- rownames(counts)
  attr(counts, "observation_weights") <- weights
  attr(counts, "soft_minobs") <- list(minobs = minobs, rule = "row_total_over_minobs")
  counts
}

make_eval_graph <- function(counts,
                            wavelength,
                            true_landscape,
                            dt,
                            beta,
                            local_shell_depth = 0L,
                            global_extra_shell = 2L,
                            max_nodes = 150000L) {
  k <- parse_karyotype_matrix(rownames(counts))
  max_cn <- max(k, na.rm = TRUE) + local_shell_depth + global_extra_shell
  data <- alfak2::prepare_alfak2_data(counts, dt = dt, beta = beta)
  graph <- alfak2::build_karyotype_graph(
    data,
    transition_kernel = "exact",
    shell_depth = local_shell_depth + global_extra_shell,
    min_cn = 0,
    max_cn = as.integer(max_cn),
    max_nodes = max_nodes
  )
  labels <- as.character(graph$labels)
  truth <- true_fitness_for_labels(labels, wavelength, true_landscape)
  observed_weight <- rowSums(data$counts)
  observed_weight <- observed_weight[match(labels, names(observed_weight))]
  nodes <- data.frame(
    node_id = seq_along(labels),
    karyotype = labels,
    support_distance = as.integer(graph$support_distance),
    support_tier = as.character(graph$support_tier),
    truth = truth,
    eval_weight = as.numeric(observed_weight),
    stringsAsFactors = FALSE
  )
  edge_keep <- as.integer(graph$support_distance[graph$edge_to]) ==
    as.integer(graph$support_distance[graph$edge_from]) + 1L &
    as.integer(graph$support_distance[graph$edge_to]) <= 2L
  edges <- data.frame(
    from = as.integer(graph$edge_from[edge_keep]),
    to = as.integer(graph$edge_to[edge_keep]),
    parent_karyotype = labels[graph$edge_from[edge_keep]],
    child_karyotype = labels[graph$edge_to[edge_keep]],
    parent_distance = as.integer(graph$support_distance[graph$edge_from[edge_keep]]),
    child_distance = as.integer(graph$support_distance[graph$edge_to[edge_keep]]),
    stringsAsFactors = FALSE
  )
  list(graph = graph, nodes = nodes, edges = edges)
}

attach_predictions <- function(eval_graph, predictions) {
  nodes <- eval_graph$nodes
  pred_idx <- match(nodes$karyotype, as.character(predictions$karyotype))
  nodes$pred <- NA_real_
  nodes$pred_sd <- NA_real_
  nodes$pred_anchor_calibrated <- NA_real_
  nodes$prediction_status <- "missing"
  ok <- !is.na(pred_idx)
  if (any(ok)) {
    nodes$pred[ok] <- as.numeric(predictions$fitness_mean[pred_idx[ok]])
    nodes$pred_anchor_calibrated[ok] <- nodes$pred[ok]
    if ("fitness_sd" %in% names(predictions)) {
      nodes$pred_sd[ok] <- as.numeric(predictions$fitness_sd[pred_idx[ok]])
    }
    if ("prediction_status" %in% names(predictions)) {
      nodes$prediction_status[ok] <- as.character(predictions$prediction_status[pred_idx[ok]])
    }
  }
  edges <- eval_graph$edges
  if (nrow(edges)) {
    from_idx <- match(edges$parent_karyotype, nodes$karyotype)
    to_idx <- match(edges$child_karyotype, nodes$karyotype)
    edges$pred_gradient <- nodes$pred[to_idx] - nodes$pred[from_idx]
    edges$truth_gradient <- nodes$truth[to_idx] - nodes$truth[from_idx]
  } else {
    edges$pred_gradient <- numeric(0)
    edges$truth_gradient <- numeric(0)
  }
  list(nodes = nodes, edges = edges)
}

full_landscape_eval <- function(predictions, wavelength, true_landscape) {
  if (is.null(predictions) || !nrow(predictions)) {
    return(list(nodes = data.frame(), edges = data.frame()))
  }
  predictions <- predictions[!is.na(predictions$karyotype) & nzchar(as.character(predictions$karyotype)), , drop = FALSE]
  predictions <- predictions[!duplicated(as.character(predictions$karyotype)), , drop = FALSE]
  if (!nrow(predictions)) return(list(nodes = data.frame(), edges = data.frame()))
  truth <- true_fitness_for_labels(as.character(predictions$karyotype), wavelength, true_landscape)
  nodes <- data.frame(
    node_id = seq_len(nrow(predictions)),
    karyotype = as.character(predictions$karyotype),
    support_distance = 0L,
    support_tier = "full_lscape",
    truth = truth,
    eval_weight = NA_real_,
    pred = as.numeric(predictions$fitness_mean),
    pred_sd = if ("fitness_sd" %in% names(predictions)) as.numeric(predictions$fitness_sd) else NA_real_,
    pred_anchor_calibrated = as.numeric(predictions$fitness_mean),
    prediction_status = if ("prediction_status" %in% names(predictions)) as.character(predictions$prediction_status) else "predicted",
    stringsAsFactors = FALSE
  )
  list(nodes = nodes, edges = data.frame(
    from = integer(), to = integer(), parent_karyotype = character(), child_karyotype = character(),
    parent_distance = integer(), child_distance = integer(), pred_gradient = numeric(), truth_gradient = numeric(),
    stringsAsFactors = FALSE
  ))
}

coerce_alfak2_predictions <- function(fit) {
  x <- alfak2::summarize_alfak2(fit, layer = "global")
  data.frame(
    karyotype = as.character(x$karyotype),
    fitness_mean = as.numeric(x$fitness_mean),
    fitness_sd = if ("fitness_sd" %in% names(x)) as.numeric(x$fitness_sd) else NA_real_,
    prediction_status = if ("prediction_status" %in% names(x)) as.character(x$prediction_status) else "predicted",
    stringsAsFactors = FALSE
  )
}

coerce_alfakR_predictions <- function(landscape_path) {
  x <- readRDS(landscape_path)
  if (!is.data.frame(x) || !"k" %in% names(x)) {
    stop("alfakR landscape output has an unexpected schema.", call. = FALSE)
  }
  data.frame(
    karyotype = as.character(x$k),
    fitness_mean = if ("mean" %in% names(x)) as.numeric(x$mean) else if ("median" %in% names(x)) as.numeric(x$median) else NA_real_,
    fitness_sd = if ("sd" %in% names(x)) as.numeric(x$sd) else NA_real_,
    prediction_status = if ("fq" %in% names(x) && "nn" %in% names(x)) {
      ifelse(x$fq %in% TRUE, "fq", ifelse(x$nn %in% TRUE, "nn", "other"))
    } else {
      "predicted"
    },
    stringsAsFactors = FALSE
  )
}

alfak_original_safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = method))
}

alfak_original_R2R <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (length(obs) < 2L) return(NA_real_)
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  denom <- sum((obs - mean(obs))^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  1 - sum((pred - obs)^2) / denom
}

alfak_original_subset_metrics <- function(truth, pred, keep) {
  keep <- keep & is.finite(truth) & is.finite(pred)
  list(
    r = alfak_original_safe_cor(pred[keep], truth[keep], method = "pearson"),
    rho = alfak_original_safe_cor(pred[keep], truth[keep], method = "spearman"),
    R = alfak_original_R2R(truth[keep], pred[keep])
  )
}

alfak_original_landscape_summary <- function(nodes, result) {
  if (!is.data.frame(nodes) || !nrow(nodes)) return(data.frame())
  if (!all(c("truth", "pred", "support_distance") %in% names(nodes))) return(data.frame())
  truth <- as.numeric(nodes$truth)
  pred <- as.numeric(nodes$pred)
  d <- suppressWarnings(as.integer(nodes$support_distance))
  finite <- is.finite(truth) & is.finite(pred)
  fq <- finite & d == 0L
  nn <- finite & d == 1L
  d2 <- finite & !(d %in% c(0L, 1L))
  all_m <- alfak_original_subset_metrics(truth, pred, finite)
  fq_m <- alfak_original_subset_metrics(truth, pred, fq)
  nn_m <- alfak_original_subset_metrics(truth, pred, nn)
  d2_m <- alfak_original_subset_metrics(truth, pred, d2)
  data.frame(
    r = all_m$r,
    rfq = fq_m$r,
    rnn = nn_m$r,
    rd2 = d2_m$r,
    rho = all_m$rho,
    rhofq = fq_m$rho,
    rhonn = nn_m$rho,
    rhod2 = d2_m$rho,
    R = all_m$R,
    Rfq = fq_m$R,
    Rnn = nn_m$R,
    Rd2 = d2_m$R,
    Rxv = suppressWarnings(as.numeric(result$cross_validation_R2 %||% NA_real_))[1L],
    nfq = sum(fq, na.rm = TRUE),
    nall = sum(finite, na.rm = TRUE),
    nnn = sum(nn, na.rm = TRUE),
    nd2 = sum(d2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

prediction_lscape <- function(predictions) {
  if (!is.data.frame(predictions) || !nrow(predictions)) return(data.frame())
  out <- data.frame(
    k = as.character(predictions$karyotype),
    mean = as.numeric(predictions$fitness_mean),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$mean) & !is.na(out$k) & nzchar(out$k), , drop = FALSE]
  out <- out[!duplicated(out$k), , drop = FALSE]
  out
}

match_time_columns <- function(times_needed, available_times) {
  vapply(times_needed, function(target) {
    exact <- which(abs(available_times - target) <= 1e-8)
    if (length(exact)) exact[[1L]] else which.min(abs(available_times - target))
  }, integer(1L))
}

normalize_count_vector <- function(x) {
  nm <- names(x)
  x <- suppressWarnings(as.numeric(x))
  names(x) <- nm
  x[!is.finite(x)] <- 0
  total <- sum(x)
  if (!is.finite(total) || total <= 0) {
    out <- rep(0, length(x))
    names(out) <- nm
    return(out)
  }
  x / total
}

safe_centroid_angle <- function(a, b) {
  ma <- sqrt(sum(a^2))
  mb <- sqrt(sum(b^2))
  if (!is.finite(ma) || !is.finite(mb) || ma == 0 || mb == 0) return(NaN)
  cosv <- sum(a * b) / (ma * mb)
  cosv <- max(-1, min(1, cosv))
  180 * acos(cosv) / pi
}

safe_cosine <- function(a, b) {
  ma <- sqrt(sum(a^2))
  mb <- sqrt(sum(b^2))
  if (!is.finite(ma) || !is.finite(mb) || ma == 0 || mb == 0) return(NaN)
  sum(a * b) / (ma * mb)
}

safe_overlap <- function(a, b) {
  denom <- min(sum(a), sum(b))
  if (!is.finite(denom) || denom <= 0) return(NaN)
  sum(pmin(a, b)) / denom
}

safe_wasserstein <- function(k, a, b) {
  if (!requireNamespace("transport", quietly = TRUE)) return(NaN)
  tryCatch({
    transport::wasserstein(
      transport::wpp(k, mass = a),
      transport::wpp(k, mass = b)
    )
  }, error = function(e) NaN)
}

average_prediction_outputs <- function(preds, lscape_k, measure_tps) {
  mats <- lapply(preds, function(y) {
    y_times <- suppressWarnings(as.numeric(y$time))
    idx <- match_time_columns(measure_tps - measure_tps[[1L]], y_times)
    if (anyNA(idx)) stop("Could not match ABM prediction output times.", call. = FALSE)
    mat <- as.matrix(y[idx, setdiff(colnames(y), "time"), drop = FALSE])
    out <- matrix(0, nrow = length(measure_tps), ncol = length(lscape_k),
                  dimnames = list(as.character(measure_tps), lscape_k))
    common <- intersect(colnames(mat), lscape_k)
    if (length(common)) out[, common] <- mat[, common, drop = FALSE]
    out
  })
  Reduce(`+`, mats) / length(mats)
}

alfak_original_forward_metrics <- function(predictions, xi, cfg, run_dir) {
  reps <- as.integer(cfg$forward_prediction_reps %||% 0L)
  if (!is.finite(reps) || reps <= 0L) return(data.frame())
  lscape <- prediction_lscape(predictions)
  if (!nrow(lscape)) return(data.frame())
  all_times <- suppressWarnings(as.numeric(colnames(xi$abm_output$x)))
  measure_tps <- as.numeric(cfg$passage_times)
  if (length(measure_tps) != 2L || !numeric_vec_equal(measure_tps, c(0, 180))) {
    stop("ALFA-K aligned forward metrics require passage_times exactly 0,180.", call. = FALSE)
  }
  idx <- match_time_columns(measure_tps, all_times)
  if (anyNA(idx)) {
    stop("Could not match passage_times to ground-truth count columns for forward metrics.", call. = FALSE)
  }
  start_time <- measure_tps[[1L]]
  report_times <- measure_tps[measure_tps > start_time] - start_time
  counts <- as.matrix(xi$abm_output$x[, idx, drop = FALSE])
  colnames(counts) <- as.character(measure_tps)

  x0_counts <- counts[, 1L]
  names(x0_counts) <- rownames(counts)
  x0_counts <- x0_counts[names(x0_counts) %in% lscape$k]
  x0_counts <- normalize_count_vector(x0_counts)
  x0 <- rep(0, nrow(lscape))
  names(x0) <- lscape$k
  if (length(x0_counts)) x0[names(x0_counts)] <- x0_counts

  pred_dir <- file.path(run_dir, "abm_preds_aligned")
  dir.create(pred_dir, recursive = TRUE, showWarnings = FALSE)
  pred_outputs <- vector("list", reps)
  for (rep_idx in seq_len(reps)) {
    pred_path <- file.path(pred_dir, paste0("rep_", pad2(rep_idx), ".Rds"))
    if (file.exists(pred_path)) {
      pred_outputs[[rep_idx]] <- readRDS(pred_path)
    } else {
      y <- alfakR::predict_evo(
        lscape = lscape,
        p = cfg$pmis,
        times = c(0, report_times),
        x0 = x0,
        prediction_type = "ABM",
        abm_pop_size = 2e5,
        abm_max_pop = 2e7,
        abm_delta_t = 0.1,
        abm_culling_survival = 0.01,
        abm_record_interval = 10,
        abm_seed = rep_idx
      )
      saveRDS(y, pred_path)
      pred_outputs[[rep_idx]] <- y
    }
  }

  pred_mat <- average_prediction_outputs(pred_outputs, lscape$k, measure_tps)
  truth_mat <- matrix(0, nrow = length(measure_tps), ncol = nrow(lscape),
                      dimnames = list(as.character(measure_tps), lscape$k))
  common_truth <- intersect(rownames(counts), lscape$k)
  for (i in seq_along(measure_tps)) {
    if (length(common_truth)) {
      truth_mat[i, common_truth] <- normalize_count_vector(counts[common_truth, i])
    }
  }

  y0 <- truth_mat[1L, ]
  k <- parse_karyotype_matrix(lscape$k)
  k0 <- colSums(as.numeric(y0) * k)
  rows <- vector("list", length(measure_tps))
  for (i in seq_along(measure_tps)) {
    y_true <- truth_mat[i, ]
    y_pred <- pred_mat[i, ]
    ky <- colSums(as.numeric(y_true) * k)
    kp <- colSums(as.numeric(y_pred) * k)
    rows[[i]] <- data.frame(
      passage = i - 1L,
      eval_time = measure_tps[[i]],
      overlap_prediction = safe_overlap(y_pred, y_true),
      overlap_baseline = safe_overlap(y0, y_true),
      cosine_prediction = safe_cosine(y_pred, y_true),
      cosine_baseline = safe_cosine(y0, y_true),
      euclidean_prediction = sqrt(sum((kp - ky)^2)),
      euclidean_baseline = sqrt(sum((k0 - ky)^2)),
      wasserstein_prediction = safe_wasserstein(k, y_pred, y_true),
      wasserstein_baseline = safe_wasserstein(k, y0, y_true),
      angle = safe_centroid_angle(ky - k0, kp - k0),
      stringsAsFactors = FALSE
    )
  }
  wide <- rbind_fill(rows)
  metrics <- c("overlap", "cosine", "euclidean", "wasserstein")
  out <- lapply(metrics, function(metric) {
    prediction <- wide[[paste0(metric, "_prediction")]]
    baseline <- wide[[paste0(metric, "_baseline")]]
    win <- if (metric %in% c("overlap", "cosine")) prediction > baseline else prediction < baseline
    data.frame(
      passage = wide$passage,
      eval_time = wide$eval_time,
      metric = metric,
      prediction = prediction,
      baseline = baseline,
      angle = wide$angle,
      win = win,
      stringsAsFactors = FALSE
    )
  })
  rbind_fill(out)
}

run_alfak2_one <- function(row, selected, cfg, run_dir) {
  counts <- prepare_alfak2_counts(
    selected$counts,
    input_mode = as.character(row$input_mode[[1L]]),
    soft_minobs = suppressWarnings(as.integer(row$soft_minobs[[1L]]))
  )
  k <- parse_karyotype_matrix(rownames(selected$counts))
  max_cn <- max(k, na.rm = TRUE) + cfg$alfak2_local_shell_depth + cfg$alfak2_global_extra_shell
  started <- proc.time()[["elapsed"]]
  set.seed(as.integer(row$fit_seed[[1L]]))
  fit <- alfak2::fit_alfak2(
    counts,
    dt = selected$dt,
    beta = cfg$pmis,
    transition_kernel = "exact",
    local_shell_depth = cfg$alfak2_local_shell_depth,
    global_extra_shell = cfg$alfak2_global_extra_shell,
    min_cn = 0,
    max_cn = as.integer(max_cn),
    max_nodes = cfg$alfak2_max_nodes,
    lambda_l_grid = cfg$alfak2_lambda_l_grid,
    lambda_e_grid = cfg$alfak2_lambda_e_grid,
    sigma_obs_grid = cfg$alfak2_sigma_obs_grid,
    graph_edge_weight = "normalized",
    extrapolation_method = as.character(row$extrapolation_method[[1L]]),
    max_prediction_distance = 2,
    anchor_count_reference = if (identical(as.character(row$input_mode[[1L]]), "soft_minobs")) as.numeric(row$soft_minobs[[1L]]) else NULL,
    input_depth = "raw",
    observation_weight_mode = "fractional_count",
    alfakR_scale = TRUE,
    n0 = cfg$n0,
    nb = cfg$nb,
    correct_efflux = FALSE,
    control = list(eval.max = cfg$alfak2_eval_max, iter.max = cfg$alfak2_iter_max),
    retry_control = list(eval.max = cfg$alfak2_retry_eval_max, iter.max = cfg$alfak2_retry_iter_max)
  )
  fit_path <- file.path(run_dir, "fit.rds")
  saveRDS(fit, fit_path)
  list(
    status = "ok",
    failure_status = "ok",
    error_message = NA_character_,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    predictions = coerce_alfak2_predictions(fit),
    fit_path = fit_path,
    dependency_status = fit$global$diagnostics$dependency_status %||% "ok",
    cross_validation_R2 = NA_real_
  )
}

run_alfakR_one <- function(row, selected, cfg, run_dir) {
  started <- proc.time()[["elapsed"]]
  set.seed(as.integer(row$fit_seed[[1L]]))
  yi <- list(x = data.frame(selected$counts, check.names = FALSE), dt = cfg$alfakR_dt)
  suppressMessages(
    alfakR::alfak(
      yi = yi,
      outdir = run_dir,
      passage_times = selected$passage_times,
      minobs = as.integer(row$minobs[[1L]]),
      nboot = cfg$nboot,
      n0 = cfg$n0,
      nb = cfg$nb,
      pm = cfg$pmis,
      nn_prior = as.character(row$NN_prior[[1L]])
    )
  )
  landscape_path <- file.path(run_dir, "landscape.Rds")
  xval <- safe_read_rds(file.path(run_dir, "xval.Rds"))
  list(
    status = "ok",
    failure_status = "ok",
    error_message = NA_character_,
    runtime_seconds = proc.time()[["elapsed"]] - started,
    predictions = coerce_alfakR_predictions(landscape_path),
    fit_path = landscape_path,
    dependency_status = "ok",
    cross_validation_R2 = suppressWarnings(as.numeric(xval))[1L] %||% NA_real_
  )
}

add_row_metadata <- function(metrics, row, result) {
  if (!nrow(metrics)) return(metrics)
  metrics$failure_status <- result$failure_status
  metrics$fit_status <- result$status
  metrics$error_message <- result$error_message
  metrics$output_path <- result$fit_path
  metrics$dependency_status <- result$dependency_status %||% NA_character_
  cbind(as.data.frame(row, stringsAsFactors = FALSE), metrics)
}

run_one_task <- function(row, cfg, dirs, force = FALSE) {
  cache_path <- file.path(dirs$cache, paste0(as.character(row$run_id[[1L]]), ".rds"))
  if (!force && file.exists(cache_path)) {
    cached <- safe_read_rds(cache_path)
    if (cache_matches_current_config(cached, cfg)) return(cached)
  }

  gt_path <- as.character(row$ground_truth_rds[[1L]])
  generate_ground_truth(row, ground_truth_times = cfg$ground_truth_times, force = FALSE)
  xi <- readRDS(gt_path)
  selected <- select_passage_counts(xi, passage_times = cfg$passage_times, ntp = cfg$ntp)

  run_dir <- file.path(dirs$runs, as.character(row$run_id[[1L]]))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  eval_graph <- make_eval_graph(
    counts = selected$counts,
    wavelength = as.numeric(row$wavelength[[1L]]),
    true_landscape = xi$true_landscape,
    dt = selected$dt,
    beta = cfg$pmis,
    local_shell_depth = cfg$alfak2_local_shell_depth,
    global_extra_shell = cfg$alfak2_global_extra_shell,
    max_nodes = cfg$alfak2_max_nodes
  )

  result <- tryCatch({
    if (identical(as.character(row$package[[1L]]), "alfak2")) {
      run_alfak2_one(row, selected, cfg, run_dir)
    } else {
      run_alfakR_one(row, selected, cfg, run_dir)
    }
  }, error = function(e) {
    list(
      status = "failed",
      failure_status = "failed",
      error_message = conditionMessage(e),
      runtime_seconds = NA_real_,
      predictions = data.frame(karyotype = character(), fitness_mean = numeric(), fitness_sd = numeric()),
      fit_path = NA_character_,
      dependency_status = "failed",
      cross_validation_R2 = NA_real_
    )
  })

  attached <- attach_predictions(eval_graph, result$predictions)
  eval_shells <- supported_second_layer_shells(c("all_eval", "all_nearfield", "d0", "d1", "d2"))
  metrics_eval <- if (length(eval_shells)) {
    alfak2:::second_layer_metric_table(
      attached$nodes,
      attached$edges,
      shells = eval_shells,
      runtime_seconds = result$runtime_seconds,
      failure_status = result$failure_status
    )
  } else {
    data.frame()
  }
  full_eval <- full_landscape_eval(result$predictions, as.numeric(row$wavelength[[1L]]), xi$true_landscape)
  full_shells <- supported_second_layer_shells("full_lscape")
  metrics_full <- if (length(full_shells)) {
    alfak2:::second_layer_metric_table(
      full_eval$nodes,
      full_eval$edges,
      shells = full_shells,
      runtime_seconds = result$runtime_seconds,
      failure_status = result$failure_status
    )
  } else {
    data.frame()
  }
  metrics <- rbind_fill(list(metrics_eval, metrics_full))
  metrics <- add_row_metadata(metrics, row, result)
  alfak_landscape <- alfak_original_landscape_summary(attached$nodes, result)
  alfak_landscape <- add_row_metadata(alfak_landscape, row, result)
  alfak_forward <- tryCatch({
    alfak_original_forward_metrics(result$predictions, xi, cfg, run_dir)
  }, error = function(e) {
    data.frame(
      passage = integer(),
      eval_time = numeric(),
      metric = character(),
      prediction = numeric(),
      baseline = numeric(),
      angle = numeric(),
      win = logical(),
      stringsAsFactors = FALSE
    )
  })
  alfak_forward <- add_row_metadata(alfak_forward, row, result)
  out <- list(
    row = row,
    result = result,
    metrics = metrics,
    alfak_original_landscape = alfak_landscape,
    alfak_original_forward = alfak_forward,
    selected_metadata = list(
      passage_times = selected$passage_times,
      dt = selected$dt,
      n_karyotypes = nrow(selected$counts)
    ),
    benchmark_config = list(
      ground_truth_times = cfg$ground_truth_times,
      passage_times = cfg$passage_times,
      alfakR_dt = cfg$alfakR_dt,
      forward_prediction_reps = cfg$forward_prediction_reps
    )
  )
  saveRDS(out, cache_path)
  out
}

generate_ground_truth_mode <- function(indices, cfg, force = FALSE) {
  gt <- indices$ground_truth_index
  idx <- arg_integer(cfg$args, "ground_truth_index", NA_integer_)
  if (is.finite(idx)) {
    row <- gt[gt$ground_truth_index == idx, , drop = FALSE]
    if (!nrow(row)) stop("No ground-truth row for --ground-truth-index=", idx, call. = FALSE)
    generate_ground_truth(row, ground_truth_times = cfg$ground_truth_times, force = force)
    return(invisible(TRUE))
  }
  for (i in seq_len(nrow(gt))) {
    message(sprintf("Generating ground truth %d/%d: %s", i, nrow(gt), gt$ground_truth_id[[i]]))
    generate_ground_truth(gt[i, , drop = FALSE], ground_truth_times = cfg$ground_truth_times, force = force)
  }
  invisible(TRUE)
}

fit_task_mode <- function(indices, cfg, dirs, force = FALSE) {
  task_index <- arg_integer(cfg$args, "task_index", NA_integer_)
  if (!is.finite(task_index)) stop("--task-index is required for --mode=fit-task.", call. = FALSE)
  row <- indices$run_index[indices$run_index$task_id == task_index, , drop = FALSE]
  if (!nrow(row)) stop("No run-index row for --task-index=", task_index, call. = FALSE)
  out <- run_one_task(row, cfg, dirs, force = force)
  message(sprintf("Finished task %d with status: %s", task_index, out$result$status))
  invisible(out)
}

method_label <- function(x) {
  ifelse(
    x$package == "alfak2",
    ifelse(
      x$input_mode == "full",
      paste("alfak2", "full", x$extrapolation_method, sep = ":"),
      paste("alfak2", paste0("soft_minobs", x$soft_minobs), x$extrapolation_method, sep = ":")
    ),
    paste("alfakR", paste0("minobs", x$minobs), ifelse(x$NN_prior == "none", "None", x$NN_prior), sep = ":")
  )
}

read_cache_metrics <- function(cache_paths) {
  metrics <- vector("list", length(cache_paths))
  for (i in seq_along(cache_paths)) {
    x <- safe_read_rds(cache_paths[[i]])
    if (!is.null(x) && is.data.frame(x$metrics) && nrow(x$metrics)) {
      metrics[[i]] <- x$metrics
    }
    if (i %% 500L == 0L) {
      message(sprintf("Read metrics from %d/%d run caches", i, length(cache_paths)))
    }
  }
  rbind_fill(metrics)
}

read_cache_table <- function(cache_paths, field) {
  out <- vector("list", length(cache_paths))
  for (i in seq_along(cache_paths)) {
    x <- safe_read_rds(cache_paths[[i]])
    if (!is.null(x) && is.data.frame(x[[field]]) && nrow(x[[field]])) {
      out[[i]] <- x[[field]]
    }
    if (i %% 500L == 0L) {
      message(sprintf("Read %s from %d/%d run caches", field, i, length(cache_paths)))
    }
  }
  rbind_fill(out)
}

landscape_wide_to_long <- function(x) {
  metric_cols <- intersect(
    c("r", "rfq", "rnn", "rd2", "rho", "rhofq", "rhonn", "rhod2", "R", "Rfq", "Rnn", "Rd2", "Rxv", "nfq"),
    names(x)
  )
  meta_cols <- setdiff(names(x), metric_cols)
  rows <- lapply(metric_cols, function(metric) {
    y <- x
    y$metric <- metric
    y$value <- suppressWarnings(as.numeric(y[[metric]]))
    y[, c(meta_cols, "metric", "value"), drop = FALSE]
  })
  rbind_fill(rows)
}

forward_wide_to_long <- function(x) {
  if (!nrow(x)) return(data.frame())
  prediction <- x
  prediction$value_type <- "prediction"
  prediction$value <- suppressWarnings(as.numeric(prediction$prediction))
  baseline <- x
  baseline$value_type <- "baseline"
  baseline$value <- suppressWarnings(as.numeric(baseline$baseline))
  rbind_fill(list(prediction, baseline))
}

make_group_key <- function(x, group_cols) {
  key_data <- x[, group_cols, drop = FALSE]
  for (nm in names(key_data)) {
    value <- key_data[[nm]]
    value <- as.character(value)
    value[is.na(value) | !nzchar(value)] <- "<NA>"
    key_data[[nm]] <- value
  }
  interaction(key_data, drop = TRUE, lex.order = TRUE)
}

metric_value_column <- function(metrics) {
  candidates <- c("value", "metric_value")
  hit <- candidates[candidates %in% names(metrics)]
  if (!length(hit)) {
    stop("Metrics table has no numeric metric value column. Expected one of: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[[1L]]
}

summarize_metric_groups <- function(metrics, group_cols, value_col) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(metrics)
    dt[, .metric_value := suppressWarnings(as.numeric(get(value_col)))]
    dt[, .finite_metric := is.finite(.metric_value)]
    out <- dt[, {
      vals <- .metric_value[.finite_metric]
      list(
        n_runs = data.table::uniqueN(run_id),
        n_success = data.table::uniqueN(run_id[fit_status == "ok"]),
        n_values = length(vals),
        mean = if (length(vals)) mean(vals) else NA_real_,
        sd = if (length(vals) > 1L) stats::sd(vals) else NA_real_,
        median = if (length(vals)) stats::median(vals) else NA_real_,
        q25 = if (length(vals)) stats::quantile(vals, 0.25, names = FALSE) else NA_real_,
        q75 = if (length(vals)) stats::quantile(vals, 0.75, names = FALSE) else NA_real_
      )
    }, by = group_cols]
    drop_cols <- intersect(c(".metric_value", ".finite_metric"), names(out))
    if (length(drop_cols)) out[, (drop_cols) := NULL]
    return(as.data.frame(out))
  }

  metric_values <- suppressWarnings(as.numeric(metrics[[value_col]]))
  finite_metric <- is.finite(metric_values)
  key <- make_group_key(metrics, group_cols)
  rows <- lapply(levels(key), function(k) {
    idx_all <- which(key == k)
    idx <- idx_all[finite_metric[idx_all]]
    all_x <- metrics[idx_all, , drop = FALSE]
    vals <- metric_values[idx]
    data.frame(
      all_x[1, group_cols, drop = FALSE],
      n_runs = length(unique(all_x$run_id)),
      n_success = length(unique(all_x$run_id[all_x$fit_status == "ok"])),
      n_values = length(vals),
      mean = if (length(vals)) mean(vals) else NA_real_,
      sd = if (length(vals) > 1L) stats::sd(vals) else NA_real_,
      median = if (length(vals)) stats::median(vals) else NA_real_,
      q25 = if (length(vals)) stats::quantile(vals, 0.25, names = FALSE) else NA_real_,
      q75 = if (length(vals)) stats::quantile(vals, 0.75, names = FALSE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  rbind_fill(rows)
}

summarize_status_groups <- function(metrics, status_cols) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::as.data.table(metrics)
    out <- dt[, list(n_runs = data.table::uniqueN(run_id)), by = status_cols]
    return(as.data.frame(out))
  }
  status_key <- make_group_key(metrics, status_cols)
  do.call(rbind, lapply(levels(status_key), function(k) {
    x <- metrics[status_key == k, , drop = FALSE]
    data.frame(x[1, status_cols, drop = FALSE], n_runs = length(unique(x$run_id)), stringsAsFactors = FALSE)
  }))
}

summarize_mode <- function(indices, dirs) {
  cache_paths <- list.files(dirs$cache, pattern = "[.]rds$", full.names = TRUE)
  metrics <- read_cache_metrics(cache_paths)
  write_tsv(metrics, file.path(dirs$tables, "metrics_long.tsv"))
  if (!nrow(metrics)) {
    warning("No completed metrics found.")
    return(invisible(metrics))
  }
  metrics$method_label <- method_label(metrics)
  value_col <- metric_value_column(metrics)
  group_cols <- c(
    "sample_depth", "wavelength", "package", "input_mode", "soft_minobs",
    "extrapolation_method", "minobs", "NN_prior", "method_label",
    "shell", "prediction_scale", "metric"
  )
  summary <- summarize_metric_groups(metrics, group_cols, value_col)
  write_tsv(summary, file.path(dirs$tables, "summary_by_depth_wavelength_method_metric.tsv"))

  status_cols <- c("sample_depth", "wavelength", "package", "method_label", "fit_status", "failure_status")
  status <- summarize_status_groups(metrics, status_cols)
  write_tsv(status, file.path(dirs$tables, "fit_status_counts.tsv"))

  alfak_landscape <- read_cache_table(cache_paths, "alfak_original_landscape")
  write_tsv(alfak_landscape, file.path(dirs$tables, "alfak_original_landscape_by_run.tsv"))
  if (nrow(alfak_landscape)) {
    alfak_landscape$method_label <- method_label(alfak_landscape)
    alfak_landscape_long <- landscape_wide_to_long(alfak_landscape)
    write_tsv(alfak_landscape_long, file.path(dirs$tables, "alfak_original_landscape_long.tsv"))
    landscape_group_cols <- c(
      "sample_depth", "wavelength", "package", "input_mode", "soft_minobs",
      "extrapolation_method", "minobs", "NN_prior", "method_label", "metric"
    )
    landscape_summary <- summarize_metric_groups(alfak_landscape_long, landscape_group_cols, "value")
    write_tsv(landscape_summary, file.path(dirs$tables, "alfak_original_landscape_summary.tsv"))
  }

  alfak_forward <- read_cache_table(cache_paths, "alfak_original_forward")
  write_tsv(alfak_forward, file.path(dirs$tables, "alfak_original_forward_prediction_metrics.tsv"))
  if (nrow(alfak_forward)) {
    alfak_forward$method_label <- method_label(alfak_forward)
    alfak_forward_long <- forward_wide_to_long(alfak_forward)
    write_tsv(alfak_forward_long, file.path(dirs$tables, "alfak_original_forward_prediction_long.tsv"))
    forward_group_cols <- c(
      "sample_depth", "wavelength", "package", "input_mode", "soft_minobs",
      "extrapolation_method", "minobs", "NN_prior", "method_label",
      "metric", "passage", "eval_time", "value_type"
    )
    forward_summary <- summarize_metric_groups(alfak_forward_long, forward_group_cols, "value")
    write_tsv(forward_summary, file.path(dirs$tables, "alfak_original_forward_prediction_summary.tsv"))
    alfak_forward_win <- alfak_forward
    alfak_forward_win$value <- as.numeric(alfak_forward_win$win)
    win_group_cols <- c(
      "sample_depth", "wavelength", "package", "input_mode", "soft_minobs",
      "extrapolation_method", "minobs", "NN_prior", "method_label",
      "metric", "passage", "eval_time"
    )
    win_summary <- summarize_metric_groups(alfak_forward_win, win_group_cols, "value")
    write_tsv(win_summary, file.path(dirs$tables, "alfak_original_forward_win_summary.tsv"))
  }
  invisible(metrics)
}

build_config <- function(args, repo_dir) {
  output_dir <- normalize_output_dir(repo_dir, arg_value(args, "output_dir", "benchmark/results/alfa2_benchmark_ground_true"))
  cfg <- list(
    args = args,
    repo_dir = repo_dir,
    output_dir = output_dir,
    alfakR_repo = normalizePath(arg_value(args, "alfakR_repo", file.path(dirname(repo_dir), "alfakR")),
                                winslash = "/", mustWork = FALSE),
    sample_depths = as.integer(arg_integer_vec(args, "sample_depths", c(1000L, 200L))),
    wavelengths = as.numeric(arg_numeric_vec(args, "wavelengths", c(0.2, 0.4, 0.8, 1.6))),
    ground_truth_times = as.numeric(arg_numeric_vec(args, "ground_truth_times", c(0, 180))),
    passage_times = as.numeric(arg_numeric_vec(args, "passage_times", c(0, 180))),
    ground_truth_reps = as.integer(arg_integer_vec(args, "ground_truth_reps", 1:5)),
    fit_repeats = as.integer(arg_integer_vec(args, "fit_repeats", 1:5)),
    soft_minobs = as.integer(arg_integer_vec(args, "soft_minobs", c(5L, 10L, 20L))),
    alfakR_priors = c("none", "empirical", "empirical_censored", "empirical_censored_weighted"),
    ntp = arg_integer(args, "ntp", 2L),
    alfakR_dt = arg_numeric(args, "alfakR_dt", 1),
    nboot = arg_integer(args, "nboot", 45L),
    forward_prediction_reps = arg_integer(args, "forward_prediction_reps", 5L),
    pmis = arg_numeric(args, "pmis", 5e-05),
    n0 = arg_numeric(args, "n0", 2e5),
    nb = arg_numeric(args, "nb", 2e7),
    ground_truth_seed_base = arg_integer(args, "ground_truth_seed_base", 720000L),
    fit_seed_base = arg_integer(args, "fit_seed_base", 910000L),
    alfak2_local_shell_depth = arg_integer(args, "alfak2_local_shell_depth", 0L),
    alfak2_global_extra_shell = arg_integer(args, "alfak2_global_extra_shell", 2L),
    alfak2_max_nodes = arg_integer(args, "alfak2_max_nodes", 150000L),
    alfak2_eval_max = arg_integer(args, "alfak2_eval_max", 500L),
    alfak2_iter_max = arg_integer(args, "alfak2_iter_max", 500L),
    alfak2_retry_eval_max = arg_integer(args, "alfak2_retry_eval_max", 2000L),
    alfak2_retry_iter_max = arg_integer(args, "alfak2_retry_iter_max", 2000L),
    alfak2_lambda_l_grid = as.numeric(arg_numeric_vec(args, "alfak2_lambda_l_grid", 0.2)),
    alfak2_lambda_e_grid = as.numeric(arg_numeric_vec(args, "alfak2_lambda_e_grid", 0.01)),
    alfak2_sigma_obs_grid = as.numeric(arg_numeric_vec(args, "alfak2_sigma_obs_grid", 0.05))
  )
  validate_benchmark_time_config(cfg)
  cfg
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (isTRUE(args$help)) {
    usage()
    quit(save = "no", status = 0)
  }
  repo_dir <- normalizePath(arg_value(args, "alfak2_repo", resolve_repo_dir()), winslash = "/", mustWork = TRUE)
  cfg <- build_config(args, repo_dir)
  dirs <- make_dirs(cfg$output_dir)
  mode <- arg_value(args, "mode", "prepare")
  force <- arg_logical(args, "force", FALSE)

  if (!identical(mode, "prepare")) {
    load_repositories(cfg$repo_dir, cfg$alfakR_repo)
  }

  if (identical(mode, "prepare")) {
    prepare_indices(cfg, dirs)
    message("Prepared benchmark indices under: ", dirs$tables)
  } else if (identical(mode, "ground-truth")) {
    indices <- load_indices(dirs)
    generate_ground_truth_mode(indices, cfg, force = force)
  } else if (identical(mode, "fit-task")) {
    indices <- load_indices(dirs)
    fit_task_mode(indices, cfg, dirs, force = force)
  } else if (identical(mode, "summarize")) {
    indices <- load_indices(dirs)
    summarize_mode(indices, dirs)
  } else if (identical(mode, "all")) {
    indices <- prepare_indices(cfg, dirs)
    generate_ground_truth_mode(indices, cfg, force = force)
    for (i in seq_len(nrow(indices$run_index))) {
      message(sprintf("Running task %d/%d", i, nrow(indices$run_index)))
      run_one_task(indices$run_index[i, , drop = FALSE], cfg, dirs, force = force)
    }
    summarize_mode(indices, dirs)
  } else {
    usage()
    stop("Unknown --mode: ", mode, call. = FALSE)
  }
}

if (sys.nframe() == 0L) main()
