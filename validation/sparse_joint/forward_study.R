# Accuracy of the closed-form forward model on realistic benchmark populations.
#
# For development landscapes and three intervals, the true population at the
# first timepoint is propagated with the true growth rates on the domain of
# karyotypes alive at either timepoint: exactly (uniformisation), with the
# first- and second-order closed forms, and with four first-order substeps.
# The last two columns compare the exact mean dynamics, and no change at all,
# with the population the simulator actually produced.
#
# Usage, from the repository root, with alfak2 installed:
#   Rscript validation/sparse_joint/forward_study.R
suppressMessages(library(alfak2))
for (f in c("landscape.R", "config.R", "generate.R")) source(file.path("benchmark/sparse_two_timepoint", f))
config <- stt_config(); t0 <- config$t0
rows <- list()
for (landscape in c("dev_L04", "dev_L05", "dev_L07", "dev_L08", "dev_L10", "dev_L11")) for (interval in c(180, 360, 1080)) {
  truth <- readRDS(stt_truth_path(config, landscape, 1L))
  a <- truth$snapshots[[as.character(t0)]]; b <- truth$snapshots[[as.character(t0 + interval)]]
  labels <- union(names(a), names(b)); K <- stt_parse_karyotypes(labels); n <- length(labels)
  kernel <- alfak2:::sparse_division_kernel(K, truth$beta)
  fitness <- stt_grf_fitness(K, truth$centroids, truth$wavelength)
  log_birth <- log(pmax(stt_effective_growth_rate(fitness, truth$abm_delta_t), 1e-6))
  n0 <- numeric(n); n0[match(names(a), labels)] <- a; log_n0 <- log(pmax(n0, 1e-300)) - log(sum(n0))
  n1 <- numeric(n); n1[match(names(b), labels)] <- b; p_abm <- n1 / sum(n1)
  norm <- function(x) { x <- exp(x - max(x)); x / sum(x) }
  exact <- norm(alfak2:::exact_forward_reference(log_n0, log_birth, kernel, interval))
  o1 <- norm(alfak2:::closed_form_forward_r(log_n0, log_birth, kernel, interval, 1L))
  o2 <- norm(alfak2:::closed_form_forward_r(log_n0, log_birth, kernel, interval, 2L))
  stepped <- function(S) { x <- log_n0; for (s in seq_len(S)) x <- alfak2:::closed_form_forward_r(x, log_birth, kernel, interval / S, 1L); norm(x) }
  s4 <- stepped(4L)
  # error on states that matter for sampling: exact frequency >= 1e-4
  err <- function(p) { use <- exact >= 1e-4; e <- abs(log(pmax(p[use], 1e-300) / exact[use])); c(tv = 0.5 * sum(abs(p - exact)), q50 = median(e), q90 = unname(quantile(e, 0.9)), max = max(e)) }
  e1 <- err(o1); e2 <- err(o2); e4 <- err(s4)
  rows[[length(rows) + 1L]] <- data.frame(landscape = landscape, class = truth$wavelength_class, days = interval, states = n, edges = length(kernel$from),
    relevant = sum(exact >= 1e-4), absent_at_t0 = sum(exact >= 1e-4 & n0 == 0),
    tv_o1 = e1[["tv"]], q90_o1 = e1[["q90"]], max_o1 = e1[["max"]], tv_o2 = e2[["tv"]], q90_o2 = e2[["q90"]], max_o2 = e2[["max"]],
    tv_s4 = e4[["tv"]], q90_s4 = e4[["q90"]], max_s4 = e4[["max"]], tv_exact_vs_abm = 0.5 * sum(abs(exact - p_abm)), tv_t0_vs_abm = 0.5 * sum(abs(norm(log_n0) - p_abm)))
}
out <- do.call(rbind, rows); num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], function(x) signif(x, 3))
options(width = 250); print(out, row.names = FALSE)
