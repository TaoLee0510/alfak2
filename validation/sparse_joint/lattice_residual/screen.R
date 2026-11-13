#!/usr/bin/env Rscript
# Screen for change C before implementing it: would a residual correlated
# along the lattice predict unobserved karyotypes better?
#
# Each dataset is refitted with the committed package (af362af) at the history
# its `base` fit chose. A residual is fitness minus the fitted trend, both
# relative to the count-weighted mean of the major karyotypes. A residual
# correlated along the lattice moves an unobserved domain karyotype towards
# the residuals of its observed neighbours, and a karyotype outside the domain
# towards those of its domain neighbours. This script applies such shifts,
# rho times the neighbours' mean, to the base predictions for rho in 0-1 and
# scores them with the development benchmark's scorer. Sources of the shift:
#
#   major             fitted residuals of the neighbouring majors
#   observed          fitted residuals of all observed neighbours
#   observed_true     true residuals of all observed neighbours (bound)
#   around_major_true the true mean residual of the unobserved neighbours of
#                     the neighbouring majors: how far the karyotypes around a
#                     major fall below the trend (bound for a random shell
#                     offset per major)
#   domain            fitted residuals of all domain neighbours, observed or
#                     not: the zero counts of unobserved ones pooled
#
# The trend is held fixed, so this is a screen, not change C. It also reports
# the lag-1 autocorrelation of the true residual and of true fitness.
#
#   Rscript validation/sparse_joint/lattice_residual/screen.R [--workers=6] [--set=rugged]
#
# Output: results/<set>/screen.csv, screen_lag1.rds and screen_evidence.rds
# (fitted against true residuals of unobserved domain karyotypes).

if (!exists("lr_script_dir")) lr_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(lr_script_dir, "setup.R"))

# Labels of the one-step neighbours of each karyotype, one row per karyotype
# and one column per chromosome and direction (NA where out of bounds).
lr_neighbour_labels <- function(labels, min_cn = 1L, max_cn = 8L) {
  ns <- asNamespace("alfak2")
  cn <- ns$parse_karyotypes(labels)
  out <- matrix(NA_character_, nrow(cn), 2L * ncol(cn))
  for (chr in seq_len(ncol(cn))) for (direction in c(-1L, 1L)) {
    moved <- cn
    moved[, chr] <- moved[, chr] + direction
    ok <- moved[, chr] >= min_cn & moved[, chr] <= max_cn
    column <- 2L * (chr - 1L) + (direction + 3L) %/% 2L
    out[ok, column] <- ns$format_karyotypes(moved[ok, , drop = FALSE])
  }
  out
}

# Mean of `value` over each row's neighbours that are in names(value), and
# their number.
lr_neighbour_mean <- function(neighbours, value) {
  hit <- matrix(value[neighbours], nrow(neighbours))
  hit[is.na(neighbours)] <- NA
  n <- rowSums(!is.na(hit))
  list(mean = ifelse(n > 0, rowSums(hit, na.rm = TRUE) / pmax(n, 1), 0), n = n)
}

lr_screen_dataset <- function(config, dataset_id, rhos) {
  dataset <- readRDS(stt_dataset_path(config, dataset_id))
  nodes <- stt_dataset_nodes_si(config, dataset)
  base <- readRDS(si_fit_path(config, dataset_id, "base"))
  if (!identical(base$status, "ok")) return(NULL)
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  fit <- alfak2::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = 1,
                                  min_cn = 1L, max_cn = 8L, history = base$diagnostics$history)
  s <- fit$summary
  trend <- fit$trend
  ns <- asNamespace("alfak2")
  # Trend, as relative fitness, for any karyotypes.
  trend_relative <- function(labels) {
    if (anyDuplicated(labels)) {
      u <- unique(labels)
      return(trend_relative(u)[match(labels, u)])
    }
    cn <- ns$parse_karyotypes(labels)
    design <- ns$sparse_joint_design(cn, trend$centre, trend$major, trend$shell_power)$design
    self <- 2 * ns$sparse_joint_self_probability(cn, dataset$beta) - 1
    self * exp(as.numeric(design %*% trend$coefficients)) - trend$reference$fitness
  }
  weight <- (s$count_t0 + s$count_t1)[s$major]
  reference <- function(rate) rate - sum(weight * stt_true_growth_rate(s$karyotype[s$major], dataset$truth)) / sum(weight)
  true_relative <- function(labels) reference(stt_true_growth_rate(labels, dataset$truth))
  domain_trend <- trend_relative(s$karyotype)
  fitted_residual <- stats::setNames(s$relative_fitness - domain_trend, s$karyotype)
  true_residual <- stats::setNames(true_relative(s$karyotype) - domain_trend, s$karyotype)
  observed <- s$count_t0 + s$count_t1 > 0
  domain_neighbours <- lr_neighbour_labels(s$karyotype)
  outside <- setdiff(nodes$karyotype, s$karyotype)
  outside_neighbours <- lr_neighbour_labels(outside)
  truth_nodes <- stt_true_growth_rate(nodes$karyotype, dataset$truth)
  rows <- list()
  # Sources of the shift of an unobserved domain karyotype: the residuals of
  # its major neighbours (fitted), or of all its observed neighbours (fitted,
  # or true: what knowing them exactly would give).
  # Mean true residual of the unobserved domain neighbours of each major:
  # how far below the trend the karyotypes around that major are.
  unobserved_residual <- true_residual[!observed]
  around_major <- lr_neighbour_mean(lr_neighbour_labels(s$karyotype[s$major]), unobserved_residual)$mean
  names(around_major) <- s$karyotype[s$major]
  for (source in c("major", "observed", "observed_true", "around_major_true", "domain")) {
    if (source == "domain") {
      # All domain neighbours, the unobserved ones with the residual their
      # zero counts gave them: pooling that evidence along the lattice.
      from <- lr_neighbour_mean(domain_neighbours, fitted_residual)
      residual <- fitted_residual
    } else if (source == "around_major_true") {
      from <- lr_neighbour_mean(domain_neighbours, around_major)
      residual <- fitted_residual
    } else {
      known <- if (source == "major") s$major else observed
      residual <- if (source == "observed_true") true_residual else fitted_residual
      from <- lr_neighbour_mean(domain_neighbours, residual[known])
    }
    for (rho in rhos) {
      domain_shift <- stats::setNames(ifelse(observed, 0, rho * from$mean), s$karyotype)
      residual_c <- stats::setNames(ifelse(observed, residual, fitted_residual + domain_shift), s$karyotype)
      outside_shift <- stats::setNames(rho * lr_neighbour_mean(outside_neighbours, residual_c)$mean, outside)
      shift <- c(domain_shift, outside_shift)
      prediction <- base$prediction
      add <- shift[prediction$karyotype]
      prediction$fitness <- prediction$fitness + ifelse(is.na(add), 0, add)
      score <- stt_score_predictions(nodes, truth_nodes, prediction, config$scoring$reference_min_count)
      m <- score$metrics[score$metrics$stratum %in% c("d0", "d1", "d2"), ]
      rows[[length(rows) + 1L]] <- data.frame(dataset_id = dataset_id, source = source, rho = rho,
                                              stratum = m$stratum, rmse = m$rmse, rmse_zero = m$rmse_zero,
                                              spearman = m$spearman, stringsAsFactors = FALSE)
    }
  }
  # Smoothness of the true residual on the lattice: each scored karyotype one
  # or two steps from the data against one of its one-step neighbours, drawn
  # at random.
  scored <- nodes$karyotype[nodes$distance > 0]
  scored_neighbours <- lr_neighbour_labels(scored)
  set.seed(dataset$seed)
  partner <- apply(scored_neighbours, 1L, function(x) { x <- x[!is.na(x)]; x[sample.int(length(x), 1L)] })
  lag1 <- data.frame(dataset_id = dataset_id, distance = nodes$distance[nodes$distance > 0],
                     residual = true_relative(scored) - trend_relative(scored),
                     partner = true_relative(partner) - trend_relative(partner),
                     fitness = true_relative(scored), partner_fitness = true_relative(partner))
  evidence <- data.frame(dataset_id = dataset_id, fitted = fitted_residual[!observed],
                         true = true_residual[!observed])
  list(scores = do.call(rbind, rows), lag1 = lag1, evidence = evidence)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  value <- function(name, default) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit)) sub("^[^=]*=", "", hit[[1L]]) else default
  }
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  config <- lr_config(set = value("set", "rugged"))
  .libPaths(c(si_library("af362af", config$library_root), .libPaths()))
  loadNamespace("alfak2")
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  index <- index[file.exists(si_fit_path(config, index$dataset_id, "base")), ]
  out <- si_parallel(index$dataset_id, function(id) lr_screen_dataset(config, id, c(0, 0.25, 0.5, 0.75, 1)),
                     as.integer(value("workers", "6")))
  scores <- merge(index[c("dataset_id", "wavelength_class", "interval")], si_bind(lapply(out, `[[`, "scores")))
  lag1 <- merge(index[c("dataset_id", "wavelength_class", "interval")], si_bind(lapply(out, `[[`, "lag1")))
  utils::write.csv(scores, file.path(config$output_dir, "screen.csv"), row.names = FALSE)
  saveRDS(lag1, file.path(config$output_dir, "screen_lag1.rds"))
  evidence <- merge(index[c("dataset_id", "wavelength_class", "interval")], si_bind(lapply(out, `[[`, "evidence")))
  saveRDS(evidence, file.path(config$output_dir, "screen_evidence.rds"))
  scores$ratio <- scores$rmse / scores$rmse_zero
  options(width = 200)
  print(stats::aggregate(ratio ~ wavelength_class + stratum + source + rho, scores[scores$stratum != "d0", ],
                         function(x) round(stats::median(x), 3)), row.names = FALSE)
  # Lag-1 autocorrelation within datasets, of the residual and of fitness.
  centre <- function(x, g) x - stats::ave(x, g)
  for (w in sort(unique(lag1$wavelength_class))) for (d in 1:2) {
    x <- lag1[lag1$wavelength_class == w & lag1$distance == d, ]
    cat(sprintf("%s d%d: lag-1 correlation of the residual %.3f, of fitness %.3f\n", w, d,
                stats::cor(centre(x$residual, x$dataset_id), centre(x$partner, x$dataset_id)),
                stats::cor(centre(x$fitness, x$dataset_id), centre(x$partner_fitness, x$dataset_id))))
  }
}
