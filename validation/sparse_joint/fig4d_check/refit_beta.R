#!/usr/bin/env Rscript
# Fig. 4d check, step 3b: do the data prefer a missegregation rate? One
# Salehi trajectory (first and last time point, as the "alfak2 2 time points"
# arm) fitted at beta = 5e-5 x {1/16, 1/4, 1, 4}, with min_cn 0 as the
# benchmark. beta is fixed in each fit, so the objectives (negative log
# marginal likelihood of the same counts) compare the rates. Also records
# delta f between each frequent karyotype and its 44 one-chromosome
# neighbours, with the neighbour's source.
#
# Runs inside the benchmark image (run_hpc.sh), one trajectory per array task:
#   Rscript refit_beta.R TASKS.tsv OUT_DIR [MULTIPLES]   (SLURM_ARRAY_TASK_ID = row;
#   MULTIPLES comma-separated, default 0.0625,0.25,1,4)

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
neighbours <- function(k) {
  v <- as.integer(strsplit(k, ".", fixed = TRUE)[[1]])
  unlist(lapply(seq_along(v), function(i) vapply(c(-1L, 1L), function(d) {
    w <- v
    w[i] <- w[i] + d
    paste(w, collapse = ".")
  }, character(1))))
}
fq <- rownames(x)[rowSums(x) >= as.numeric(task$min_obs)]
pairs <- do.call(rbind, lapply(fq, function(fk) data.frame(k_id = fk, nn = neighbours(fk), stringsAsFactors = FALSE)))
pairs <- pairs[!grepl("-", pairs$nn, fixed = TRUE), ]
keys <- unique(c(pairs$k_id, pairs$nn))

multiples <- if (length(args) > 2) as.numeric(strsplit(args[3], ",", fixed = TRUE)[[1]]) else c(1 / 16, 1 / 4, 1, 4)
rows <- lapply(multiples, function(m) {
  beta <- 5e-5 * m
  fit <- tryCatch(suppressWarnings(fit_sparse_joint(x, dt = dt, beta = beta, min_cn = 0)), error = function(e) e)
  if (inherits(fit, "error")) {
    return(data.frame(beta_multiple = m, objective = NA, history = NA, k_id = NA, nn = NA, deltaf = NA,
                      source = NA, error = conditionMessage(fit), stringsAsFactors = FALSE))
  }
  p <- predict(fit, keys)
  f <- stats::setNames(p$fitness_mean, p$karyotype)
  d <- pairs
  d$deltaf <- unname(f[d$k_id] - f[d$nn])
  in_d <- d$nn %in% fit$summary$karyotype
  obs <- d$nn %in% rownames(x)[rowSums(x) > 0]
  d$source <- ifelse(in_d & obs, "observed_D", ifelse(in_d, "unseen_D", ifelse(obs, "observed_out", "trend")))
  d$beta_multiple <- m
  d$objective <- fit$diagnostics$objective
  d$history <- fit$diagnostics$history
  d$error <- NA_character_
  d[c("beta_multiple", "objective", "history", "k_id", "nn", "deltaf", "source", "error")]
})
res <- do.call(rbind, rows)
res$fi <- task$fi
res$context <- task$context
res$treated <- task$treated
res$cells_mean <- mean(colSums(x))
res$dt <- dt
saveRDS(res, out_file)
cat(task$fi, "done\n")
