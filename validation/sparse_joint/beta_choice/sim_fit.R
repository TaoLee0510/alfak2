#!/usr/bin/env Rscript
# Choosing the missegregation rate by marginal likelihood, step 0: one fit of
# one simulated input at one rate (one Slurm array task).
#
# Inputs are the paper's ABM and missegregation simulations, built as the
# paper benchmark builds them (its R/abm.R, read only): the first and last
# time point of the window of `window` passages before day 120. The fit is
# alfak2::fit_sparse_joint(x, dt, beta = beta_multiple * 5e-5) with the
# benchmark's other settings (defaults). Records the objective (negative log
# marginal likelihood, comparable across rates for the same counts), the
# chosen history and convergence, and the landscape's accuracy against the
# true fitness on the benchmark's karyotype set (minobs 10: frequent
# karyotypes, their one-step neighbours and theirs), as the benchmark's
# landscape_scores(): r, rfq, rnn, rd2 and rescaled R^2. Also the delta f
# rows between each frequent karyotype and its neighbours.
#
#   Rscript sim_fit.R TASKS.tsv OUT_DIR    (SLURM_ARRAY_TASK_ID = row)

suppressPackageStartupMessages(library(alfak2))
args <- commandArgs(TRUE)
tasks <- utils::read.delim(args[1], stringsAsFactors = FALSE)
out_dir <- args[2]
row <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
task <- tasks[row, ]
# ABM wavelength 0.8 and missegregation rate 5e-5 share simulation ids: the
# source is part of the name.
out_file <- file.path(out_dir, sprintf("%s_%s_w%d_b%s.rds", task$source, task$sim_id, task$window,
                                      format(task$beta_multiple)))
if (file.exists(out_file)) quit(save = "no")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bench <- "/share/lab_crd/taoli/Project/alfak2-paper-benchmark"
data_dir <- Sys.getenv("PAPER_DATA", "/share/lab_crd/M010_ALFAK_2023/ALFA-K/data")
source(file.path(bench, "R", "abm.R"))
raw <- if (task$source == "abm") file.path(data_dir, "raw", "ABM") else file.path(data_dir, "raw", "ABM_varyMisseg")
sim <- read_abm_sim(task$sim_id, raw)
info <- abm_sim_info(task$sim_id)
yi <- abm_input(sim, list(window = task$window, ends = TRUE), last_time = 120)
x <- as.matrix(yi$x)
storage.mode(x) <- "double"
times <- as.numeric(colnames(yi$x)) * yi$dt
beta <- 5e-5 * task$beta_multiple

r2r <- function(obs, pred) {
  obs <- obs - mean(obs)
  pred <- pred - mean(pred)
  1 - sum((pred - obs)^2) / sum((obs - mean(obs))^2)
}

start <- proc.time()[["elapsed"]]
fit <- tryCatch(suppressWarnings(fit_sparse_joint(x, dt = diff(times), beta = beta)), error = function(e) e)
seconds <- proc.time()[["elapsed"]] - start
base <- data.frame(source = task$source, sim_id = task$sim_id, window = task$window, w = info$w,
                   true_beta = info$pmis, beta_multiple = task$beta_multiple, beta = beta, seconds = seconds,
                   stringsAsFactors = FALSE)
if (inherits(fit, "error")) {
  base$error <- conditionMessage(fit)
  saveRDS(list(summary = base, deltaf = NULL), out_file)
  quit(save = "no")
}

fq <- rownames(x)[rowSums(x) >= 10]
nn <- alfak2:::alfakR_neighbours(fq)
k <- unique(c(fq, nn, alfak2:::alfakR_neighbours(nn)))
p <- predict(fit, k)
m <- p$fitness_mean
f <- grf_fitness(k, info$w, sim$true_landscape)
isfq <- k %in% fq
isnn <- k %in% nn & !isfq
d2 <- !(isfq | isnn)
safe_cor <- function(a, b) if (length(a) > 2 && stats::sd(a) > 0 && stats::sd(b) > 0) stats::cor(a, b) else NA_real_

base$error <- NA_character_
base$objective <- fit$diagnostics$objective
base$history <- fit$diagnostics$history
base$convergence <- fit$diagnostics$convergence
base$nfq <- length(fq)
base$r <- safe_cor(m, f)
base$rfq <- safe_cor(m[isfq], f[isfq])
base$rnn <- safe_cor(m[isnn], f[isnn])
base$rd2 <- safe_cor(m[d2], f[d2])
base$R <- r2r(f, m)
base$Rfq <- if (sum(isfq) > 1) r2r(f[isfq], m[isfq]) else NA_real_
base$cells_mean <- mean(colSums(x))
base$dt <- diff(times)

lut_m <- stats::setNames(m, k)
lut_f <- stats::setNames(f, k)
seen <- rownames(x)[rowSums(x) > 0]
dd <- do.call(rbind, lapply(fq, function(fk) {
  v <- as.integer(strsplit(fk, ".", fixed = TRUE)[[1]])
  nb <- unlist(lapply(seq_along(v), function(i) vapply(c(-1L, 1L), function(s) {
    u <- v
    u[i] <- u[i] + s
    paste(u, collapse = ".")
  }, character(1))))
  nb <- nb[nb %in% k]
  if (!length(nb)) return(NULL)
  data.frame(k_id = fk, nn = nb, est = unname(lut_m[fk] - lut_m[nb]), true = unname(lut_f[fk] - lut_f[nb]),
             in_domain = nb %in% fit$summary$karyotype, observed = nb %in% seen, stringsAsFactors = FALSE)
}))
saveRDS(list(summary = base, deltaf = dd), out_file)
cat(task$sim_id, task$window, task$beta_multiple, "done in", round(seconds), "s\n")
