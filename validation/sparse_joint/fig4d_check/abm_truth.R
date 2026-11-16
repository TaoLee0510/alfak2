#!/usr/bin/env Rscript
# Fig. 4d check, step 2: is alfak2's delta f between a frequent karyotype and
# its one-chromosome neighbours right where the truth is known?
#
# Reads the ABM and missegregation fits of the paper benchmark round 05a8855
# (read only) and the true Gaussian-random-field landscapes of the paper's
# simulations. For each fit (minobs 10), each frequent karyotype (landscape
# fq) and each of its 44 neighbours: delta f estimated (landscape mean) and
# true (the paper's getf), and where alfak2's value of the neighbour comes
# from (observed_D, unseen_D, observed_out, trend; see decompose.R).
#
# Runs inside the benchmark image on HPC (run_hpc.sh):
#   Rscript abm_truth.R OUT_DIR [SOURCE]   SOURCE: abm or misseg (default both)

args <- commandArgs(TRUE)
out_dir <- args[1]
sources <- if (length(args) > 1) args[2] else c("abm", "misseg")
data_dir <- Sys.getenv("PAPER_DATA", "/share/lab_crd/M010_ALFAK_2023/ALFA-K/data")
fits_root <- "/share/lab_crd/taoli/Project/alfak2-paper-benchmark/results_05a8855/fits"
raw <- c(abm = file.path(data_dir, "raw", "ABM"), misseg = file.path(data_dir, "raw", "ABM_varyMisseg"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

grf_fitness <- function(k, wavelength, centroids) {
  kmat <- do.call(rbind, lapply(strsplit(k, ".", fixed = TRUE), as.numeric))
  scalef <- 1 / (pi * sqrt(nrow(centroids)))
  vapply(seq_len(nrow(kmat)), function(i) {
    d <- sqrt(colSums((t(centroids) - kmat[i, ])^2))
    sum(sin(d / wavelength) * scalef)
  }, numeric(1))
}
neighbours <- function(k) {
  v <- as.integer(strsplit(k, ".", fixed = TRUE)[[1]])
  unlist(lapply(seq_along(v), function(i) vapply(c(-1L, 1L), function(d) {
    w <- v
    w[i] <- w[i] + d
    paste(w, collapse = ".")
  }, character(1))))
}
wavelength <- function(sim_id) as.numeric(sub("p", ".", strsplit(sim_id, "_", fixed = TRUE)[[1]][2], fixed = TRUE))

one_sim <- function(src, sim_id) {
  sim <- readRDS(file.path(raw[[src]], paste0(sim_id, ".Rds")))
  w <- wavelength(sim_id)
  rows <- list()
  for (design in c("w2_ends", "w4_ends", "w8_ends", "w4_all", "w8_all")) {
    a2dir <- file.path(fits_root, src, sim_id, design, "alfak2", "minobs_10")
    fit <- if (file.exists(file.path(a2dir, "alfak2_fit.Rds"))) readRDS(file.path(a2dir, "alfak2_fit.Rds")) else NULL
    for (method in c("alfak2", "alfakR")) {
      lp <- file.path(fits_root, src, sim_id, design, method, "minobs_10", "landscape.Rds")
      if (!file.exists(lp)) next
      ls <- readRDS(lp)
      lut <- stats::setNames(ls$mean, ls$k)
      focal <- ls$k[ls$fq]
      if (!length(focal)) next
      d <- do.call(rbind, lapply(focal, function(fk) data.frame(k_id = fk, nn = neighbours(fk), stringsAsFactors = FALSE)))
      d <- d[!grepl("(^|\\.)0(\\.|$)", d$nn), ]
      d$est <- unname(lut[d$k_id] - lut[d$nn])
      keys <- unique(c(d$k_id, d$nn))
      truth <- stats::setNames(grf_fitness(keys, w, sim$true_landscape), keys)
      d$true <- unname(truth[d$k_id] - truth[d$nn])
      if (!is.null(fit)) {
        counts <- fit$data$counts
        lab <- if (!is.null(rownames(counts))) rownames(counts) else fit$data$labels
        seen <- lab[rowSums(counts) > 0]
        in_d <- d$nn %in% fit$summary$karyotype
        obs <- d$nn %in% seen
        d$source <- ifelse(in_d & obs, "observed_D", ifelse(in_d, "unseen_D", ifelse(obs, "observed_out", "trend")))
        d$history <- fit$settings$history
        d$dt_total <- sum(fit$data$dt)
        d$cells_mean <- mean(colSums(counts))
      } else {
        d$source <- NA_character_
        d$history <- d$dt_total <- d$cells_mean <- NA_real_
      }
      d$method <- method
      d$design <- design
      d$sim_id <- sim_id
      d$source_set <- src
      d$w <- w
      rows[[length(rows) + 1L]] <- d
    }
  }
  do.call(rbind, rows)
}

for (src in sources) {
  sims <- list.files(file.path(fits_root, src))
  res <- lapply(sims, function(s) tryCatch(one_sim(src, s), error = function(e) {
    message(src, " ", s, ": ", conditionMessage(e))
    NULL
  }))
  x <- do.call(rbind, res)
  x <- x[is.finite(x$est) & is.finite(x$true), ]
  saveRDS(x, file.path(out_dir, paste0("truth_rows_", src, ".rds")))
  cat(src, ":", length(sims), "simulations,", nrow(x), "rows\n")
}
