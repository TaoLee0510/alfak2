# Unobserved karyotypes predicted too fit at intermediate ruggedness
# (README.md in this directory).
#
# The datasets are those of two earlier studies, linked, not copied or
# regenerated: the paper-sized young populations of ../lattice_residual
# ("rugged") and the development-sized populations of ../short_interval
# ("diagnosis"), wavelengths 0.4 and 0.8 only. Fits and scores of this study
# go to its own results/<source>.

ub_dir <- if (exists("ub_script_dir")) ub_script_dir else getwd()
lr_script_dir <- normalizePath(file.path(ub_dir, "..", "lattice_residual"))
source(file.path(lr_script_dir, "setup.R"))

ub_default_output <- file.path(si_repo, "validation", "sparse_joint", "unobserved_bias", "results")

ub_sources <- list(
  rugged = file.path(lr_default_output, "rugged"),
  diagnosis = si_default_output,
  # The development benchmark's development landscapes, as copied by
  # ../short_interval/devbench.R (used by laplace_mean.R only).
  dev = file.path(si_default_output, "devbench")
)

ub_config <- function(source = c("rugged", "diagnosis", "dev")) {
  source <- match.arg(source)
  config <- if (source == "rugged") lr_config() else si_config()
  config$source <- source
  config$source_dir <- ub_sources[[source]]
  config$output_dir <- file.path(ub_default_output, source)
  config$library_root <- lr_default_output
  config$wavelengths <- c("w0.4", "w0.8")
  sparse <- function(...) {
    list(type = "sparse_joint", build = "af362af", birth_rate = 1, min_cn = 1L, max_cn = 8L,
         arguments = list(...))
  }
  config$methods <- list(
    zero = list(type = "zero"),
    # The package as committed (af362af), at its defaults.
    base = sparse(),
    # Tighter spread of first-sample abundances around balance: their
    # differences must then come from fitness.
    bsd05 = sparse(balance_sd = 0.5),
    bsd025 = sparse(balance_sd = 0.25),
    # History = true age of the population (probe; only a simulation knows it).
    age = sparse(history = "age"),
    # Without the balance prior.
    free = sparse(initial_state = "free")
  )
  ub_add_residual_probes(config)
}

# Link the source's datasets and scored nodes, and write an index of the
# wavelengths studied here.
ub_link <- function(config) {
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  for (part in c("datasets", "nodes")) {
    link <- file.path(config$output_dir, part)
    if (!file.exists(link)) file.symlink(file.path(config$source_dir, part), link)
  }
  index <- utils::read.csv(file.path(config$source_dir, "dataset_index.csv"), stringsAsFactors = FALSE)
  index <- index[index$wavelength_class %in% config$wavelengths, ]
  utils::write.csv(index, si_index_path(config), row.names = FALSE)
  cat("Linked", nrow(index), "datasets of", config$source_dir, "\n")
}

# Where the error sits: each scored karyotype by distance to the data, shell
# (distance to the nearest major) and whether it is in the model's domain,
# with the true and predicted relative fitness, the reported standard
# deviation and whether the 95% interval covers the truth.
ub_node_table <- function(config, dataset_id, method) {
  dataset <- readRDS(stt_dataset_path(config, dataset_id))
  nodes <- stt_dataset_nodes_si(config, dataset)
  fit <- readRDS(si_fit_path(config, dataset_id, method))
  if (!identical(fit$status, "ok") || is.null(fit$domain)) return(NULL)
  truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  score <- stt_score_predictions(nodes, truth, fit$prediction, config$scoring$reference_min_count)$nodes
  p <- fit$prediction[match(score$karyotype, fit$prediction$karyotype), ]
  shell <- pmin(p$shell, 3)
  group <- paste0("d", score$distance, ifelse(score$karyotype %in% fit$domain$karyotype, "_in", "_out"), "_s", shell)
  data.frame(dataset_id = dataset_id, method = method, group = group, truth = score$truth_relative,
             pred = score$pred_relative, sd = p$fitness_sd,
             covered = abs(score$pred_relative - score$truth_relative) <= stats::qnorm(0.975) * p$fitness_sd,
             stringsAsFactors = FALSE)
}

# Per dataset and method: the development benchmark's scores by distance, the
# node groups, and the fitted scales, baseline and history.
ub_stage_score <- function(config, workers) {
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  methods <- intersect(names(config$methods), list.files(file.path(config$output_dir, "fits")))
  out <- si_parallel(index$dataset_id, function(dataset_id) {
    dataset <- readRDS(stt_dataset_path(config, dataset_id))
    nodes <- stt_dataset_nodes_si(config, dataset)
    truth <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
    lapply(methods, function(method) {
      path <- si_fit_path(config, dataset_id, method)
      if (!file.exists(path)) return(NULL)
      fit <- readRDS(path)
      if (!identical(fit$status, "ok")) return(NULL)
      m <- stt_score_predictions(nodes, truth, fit$prediction, config$scoring$reference_min_count)$metrics
      m <- m[m$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse", "rmse_zero", "bias", "spearman", "coverage95")]
      d <- fit$diagnostics
      info <- data.frame(dataset_id = dataset_id, method = method,
                         history = if (is.null(d$history)) NA else d$history,
                         fitted_birth_rate = if (is.null(d$fitted_birth_rate)) NA else d$fitted_birth_rate,
                         scale_residual = if (is.null(fit$scales)) NA else fit$scales[["residual"]],
                         scale_dosage = if (is.null(fit$scales)) NA else fit$scales[["dosage"]],
                         fit_seconds = if (is.null(fit$fit_seconds)) NA else fit$fit_seconds)
      list(distance = cbind(info, m), nodes = if (method == "zero") NULL else ub_node_table(config, dataset_id, method))
    })
  }, workers)
  out <- unlist(out, recursive = FALSE)
  key <- index[c("dataset_id", "wavelength_class", "t0", "interval")]
  distance <- merge(key, si_bind(lapply(out, `[[`, "distance")))
  node_rows <- si_bind(lapply(out, `[[`, "nodes"))
  groups <- stats::aggregate(cbind(truth, pred, sd, covered) ~ dataset_id + method + group, node_rows, mean)
  groups$n <- stats::aggregate(truth ~ dataset_id + method + group, node_rows, length)$truth
  groups <- merge(key, groups)
  utils::write.csv(distance, file.path(config$output_dir, "metrics_distance.csv"), row.names = FALSE)
  utils::write.csv(groups, file.path(config$output_dir, "metrics_groups.csv"), row.names = FALSE)
  cat("Scored", nrow(distance), "rows in", config$output_dir, "\n")
}

# Probe: the residual scale held at `scale` (log birth rate per karyotype)
# instead of estimated. The package has no such option; inside a worker the
# problem builder is wrapped so that the residual scale's prior is a spike at
# `scale` (sd 0.01 on the log scale), and restored afterwards. Probe only.
ub_with_residual_scale <- function(scale, expr) {
  ns <- asNamespace("alfak2")
  original <- get("sparse_joint_problem", envir = ns)
  wrapped <- function(data, settings) {
    problem <- original(data, settings)
    problem$tmb_data$log_scale_prior_mean[4L] <- log(scale)
    problem$tmb_data$log_scale_prior_sd[4L] <- 0.01
    problem$parameters$log_scale[4L] <- log(scale)
    problem
  }
  utils::assignInNamespace("sparse_joint_problem", wrapped, ns = "alfak2")
  on.exit(utils::assignInNamespace("sparse_joint_problem", original, ns = "alfak2"))
  force(expr)
}

ub_add_residual_probes <- function(config, scales = c(0.03, 0.06, 0.12)) {
  for (s in scales) {
    name <- sprintf("res%03d", round(s * 100))
    config$methods[[name]] <- list(type = "sparse_joint", build = "af362af", birth_rate = 1, min_cn = 1L,
                                   max_cn = 8L, arguments = list(), residual_scale = s)
  }
  config
}

# si_run_method() of the short-interval study, with the residual-scale probe.
si_run_method_base <- si_run_method
si_run_method <- function(config, dataset, nodes, method) {
  if (is.null(method$residual_scale)) return(si_run_method_base(config, dataset, nodes, method))
  ub_with_residual_scale(method$residual_scale, si_run_method_base(config, dataset, nodes, method))
}
