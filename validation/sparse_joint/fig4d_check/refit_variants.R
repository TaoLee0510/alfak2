#!/usr/bin/env Rscript
# Fig. 4d check, step 3: refit one Salehi trajectory (first and last time
# point, as the benchmark's "alfak2 2 time points" arm) under variants that
# isolate what drives alfak2's delta f of unseen neighbours:
#   baseline    the benchmark's call (beta 5e-5, min_cn 0)
#   depth320    every time point drawn down to 320 cells without replacement
#               (the PDX median), seed fixed per trajectory
#   beta_quarter  beta / 4: a quarter of the missegregation influx
#   history_T   history fixed at the interval (no longer candidates)
#   free        initial_state = "free" (no balance prior)
# For each variant: f of every frequent karyotype (at least min_obs cells over
# the time points, as alfakR's fq) and of its 44 one-chromosome neighbours from
# predict(); delta f = f(frequent) - f(neighbour), with the neighbour's source.
#
# Runs inside the benchmark image (run_hpc.sh), one trajectory per array task:
#   Rscript refit_variants.R TASKS.tsv OUT_DIR   (SLURM_ARRAY_TASK_ID = row)

suppressPackageStartupMessages(library(alfak2))
args <- commandArgs(TRUE)
tasks <- utils::read.delim(args[1], stringsAsFactors = FALSE)
out_dir <- args[2]
row <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
task <- tasks[row, ]
data_dir <- Sys.getenv("PAPER_DATA", "/share/lab_crd/M010_ALFAK_2023/ALFA-K/data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0(task$fi, ".rds"))
if (file.exists(out_file)) quit(save = "no")

yi <- readRDS(file.path(data_dir, "processed", "salehi", "alfak_inputs", paste0(task$fi, ".Rds")))
x <- as.matrix(yi$x)
if (ncol(x) > 2L) x <- x[, c(1L, ncol(x)), drop = FALSE]
storage.mode(x) <- "double"
times <- as.numeric(colnames(yi$x))[c(1L, ncol(yi$x))] * yi$dt
dt <- diff(times)

downsample <- function(x, n, seed) {
  set.seed(seed)
  for (j in seq_len(ncol(x))) {
    if (sum(x[, j]) <= n) next
    cells <- rep(seq_len(nrow(x)), x[, j])
    keep <- sample(cells, n)
    x[, j] <- tabulate(keep, nrow(x))
  }
  x[rowSums(x) > 0, , drop = FALSE]
}
neighbours <- function(k) {
  v <- as.integer(strsplit(k, ".", fixed = TRUE)[[1]])
  unlist(lapply(seq_along(v), function(i) vapply(c(-1L, 1L), function(d) {
    w <- v
    w[i] <- w[i] + d
    paste(w, collapse = ".")
  }, character(1))))
}

variants <- list(
  baseline = list(),
  depth320 = list(counts = downsample(x, 320, 7000 + row)),
  beta_quarter = list(beta = 5e-5 / 4),
  history_T = list(history = dt),
  free = list(initial_state = "free"))

rows <- lapply(names(variants), function(v) {
  a <- variants[[v]]
  counts <- if (is.null(a$counts)) x else a$counts
  call <- list(counts, dt = dt, beta = if (is.null(a$beta)) 5e-5 else a$beta, min_cn = 0)
  if (!is.null(a$history)) call$history <- a$history
  if (!is.null(a$initial_state)) call$initial_state <- a$initial_state
  start <- proc.time()[["elapsed"]]
  fit <- tryCatch(suppressWarnings(do.call(fit_sparse_joint, call)), error = function(e) e)
  if (inherits(fit, "error")) {
    return(data.frame(variant = v, k_id = NA, nn = NA, deltaf = NA, source = NA, history = NA,
                      error = conditionMessage(fit), stringsAsFactors = FALSE))
  }
  # Frequent karyotypes as in the original fit (the full counts), so that the
  # same karyotypes are compared across variants.
  fq <- rownames(x)[rowSums(x) >= as.numeric(task$min_obs)]
  fq <- fq[!grepl("(^|\\.)-", fq)]
  d <- do.call(rbind, lapply(fq, function(fk) data.frame(k_id = fk, nn = neighbours(fk), stringsAsFactors = FALSE)))
  d <- d[!grepl("-", d$nn, fixed = TRUE), ]
  keys <- unique(c(d$k_id, d$nn))
  p <- predict(fit, keys)
  f <- stats::setNames(p$fitness_mean, p$karyotype)
  d$deltaf <- unname(f[d$k_id] - f[d$nn])
  seen <- rownames(counts)[rowSums(counts) > 0]
  in_d <- d$nn %in% fit$summary$karyotype
  obs <- d$nn %in% seen
  d$source <- ifelse(in_d & obs, "observed_D", ifelse(in_d, "unseen_D", ifelse(obs, "observed_out", "trend")))
  d$history <- fit$diagnostics$history
  d$cells_mean <- mean(colSums(counts))
  d$seconds <- proc.time()[["elapsed"]] - start
  d$error <- NA_character_
  d$variant <- v
  d
})
cols <- unique(unlist(lapply(rows, names)))
res <- do.call(rbind, lapply(rows, function(r) {
  for (c in setdiff(cols, names(r))) r[[c]] <- NA
  r[cols]
}))
res$fi <- task$fi
res$context <- task$context
res$treated <- task$treated
res$dt <- dt
saveRDS(res, out_file)
cat(task$fi, "done:", paste(unique(res$variant), collapse = ","), "\n")
