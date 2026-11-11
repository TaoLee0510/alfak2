# Checks of the sparse joint model on real two-timepoint samples, where there
# is no truth to score against.
#
#   Rscript validation/sparse_joint/real_cohort.R --input-dir=DIR --output-dir=DIR
#       [--birth-rate-range=0.02,1] [--beta=5e-5] [--split-min-cells=200] [--workers=2] [--skip-cohort]
#
# Inputs are ALFA-K input files (.Rds lists whose `x` is a count matrix with
# karyotype row names and the sampling days as column names). They are read
# in place; everything written goes to --output-dir, which should lie outside
# the repository when the samples are patient data.
#
# What is checked:
#   1. every sample fits, at each of three birth rates across the range;
#   2. rankings agree across those rates;
#   3. the counts the fit expects resemble the observed counts;
#   4. two halves of the cells of a sample give the same answer, within the
#      reported uncertainty;
#   5. the cohort fit, and how far it moves each sample.
# Set STT_R_LIB to the library holding the alfak2 build under test.

library_dir <- Sys.getenv("STT_R_LIB", unset = "")
if (nzchar(library_dir)) .libPaths(c(library_dir, .libPaths()))
suppressPackageStartupMessages(library(alfak2))
args <- commandArgs(TRUE)
option <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  if (grepl("=", hit[[1L]], fixed = TRUE)) sub("^[^=]*=", "", hit[[1L]]) else TRUE
}
input_dir <- option("input-dir")
output_dir <- option("output-dir")
if (is.null(input_dir) || is.null(output_dir)) stop("--input-dir and --output-dir are required.", call. = FALSE)
rate_range <- as.numeric(strsplit(option("birth-rate-range", "0.02,1"), ",", fixed = TRUE)[[1L]])
beta <- as.numeric(option("beta", "5e-5"))
split_min_cells <- as.numeric(option("split-min-cells", "200"))
workers <- as.integer(option("workers", "2"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
central_rate <- exp(mean(log(rate_range)))

files <- list.files(input_dir, pattern = "[.][Rr]ds$", full.names = TRUE)
samples <- lapply(files, function(file) {
  input <- readRDS(file)
  counts <- as.matrix(if (is.list(input) && !is.null(input$x)) input$x else input)
  storage.mode(counts) <- "integer"
  list(counts = counts[rowSums(counts) > 0, , drop = FALSE], dt = diff(as.numeric(colnames(counts))))
})
names(samples) <- sub("[.][Rr]ds$", "", basename(files))
cat("Samples:", length(samples), " birth rates:", signif(rate_range[1L], 3), "to", signif(rate_range[2L], 3),
    "per day, centre", signif(central_rate, 3), "\n")
copy_numbers <- unlist(lapply(samples, function(s) range(alfak2:::parse_karyotypes(rownames(s$counts)))))
min_cn <- min(1L, min(copy_numbers))
max_cn <- max(copy_numbers) + 1L

run <- function(items, fun) {
  parallel::mclapply(items, function(item) tryCatch(fun(item), error = function(e) conditionMessage(e)),
                     mc.cores = workers, mc.preschedule = FALSE)
}
fit_one <- function(counts, dt, ...) {
  fit_sparse_joint(counts, dt = dt, beta = beta, min_cn = min_cn, max_cn = max_cn, ...)
}

# 1-3. Single-sample fits across the range of birth rates.
solo_path <- file.path(output_dir, "solo_fits.rds")
solo <- if (file.exists(solo_path)) readRDS(solo_path) else {
  out <- run(names(samples), function(name) {
    sample <- samples[[name]]
    started <- proc.time()[["elapsed"]]
    fit <- fit_sparse_joint_rates(sample$counts, dt = sample$dt, beta = beta, min_cn = min_cn, max_cn = max_cn,
                                  birth_rate_range = rate_range)
    fit$seconds <- proc.time()[["elapsed"]] - started
    fit
  })
  names(out) <- names(samples)
  saveRDS(out, solo_path)
  out
}
solo_table <- do.call(rbind, lapply(names(samples), function(name) {
  sample <- samples[[name]]
  fit <- solo[[name]]
  row <- data.frame(sample = name, interval = sample$dt, cells_t0 = sum(sample$counts[, 1L]),
                    cells_t1 = sum(sample$counts[, 2L]), karyotypes = nrow(sample$counts),
                    status = if (is.character(fit)) fit else "ok", stringsAsFactors = FALSE)
  if (is.character(fit)) return(cbind(row, states = NA, majors = NA, dropped = NA, converged = NA,
                                       seconds = NA, rank_agreement = NA, dispersion_t1 = NA,
                                       deficit_ratio = NA, per_division_ratio = NA))
  centre <- fit$fits[[fit$central]]
  summary <- centre$summary
  prediction <- stats::predict(fit, summary$karyotype)
  seen <- summary$count_t1 > 0 | summary$expected_t1 >= 1
  minor <- !summary$major
  deficit <- vapply(fit$fits, function(f) mean(f$summary$relative_fitness[minor]), numeric(1))
  cbind(row,
        states = centre$diagnostics$n_states, majors = centre$diagnostics$n_major,
        dropped = centre$diagnostics$dropped_karyotypes,
        converged = sum(vapply(fit$fits, function(f) f$diagnostics$convergence == 0L, logical(1))),
        seconds = fit$seconds,
        # lowest Spearman correlation between the predictions at two rates
        rank_agreement = attr(prediction, "rank_agreement"),
        # Pearson dispersion of the second sample's counts around the fit
        dispersion_t1 = sum(((summary$count_t1 - summary$expected_t1)^2 / pmax(summary$expected_t1, 1e-8))[seen]) /
          sum(seen),
        # mean deficit of non-major karyotypes: highest over lowest rate, per day and per division
        deficit_ratio = deficit[length(deficit)] / deficit[1L],
        per_division_ratio = (deficit[length(deficit)] / fit$rates[length(deficit)]) / (deficit[1L] / fit$rates[1L]))
}))
utils::write.csv(solo_table, file.path(output_dir, "solo_summary.csv"), row.names = FALSE)

# 4. Split-half reproducibility at the central rate.
split_half <- function(counts, seed) {
  set.seed(seed)
  halves <- lapply(seq_len(2L), function(column) {
    cells <- rep(seq_len(nrow(counts)), counts[, column])
    first <- sample(length(cells), length(cells) %/% 2L)
    cbind(tabulate(cells[first], nrow(counts)), tabulate(cells[-first], nrow(counts)))
  })
  lapply(seq_len(2L), function(half) {
    out <- cbind(halves[[1L]][, half], halves[[2L]][, half])
    dimnames(out) <- dimnames(counts)
    storage.mode(out) <- "integer"
    out[rowSums(out) > 0, , drop = FALSE]
  })
}
eligible <- names(samples)[vapply(samples, function(s) all(colSums(s$counts) >= split_min_cells), logical(1))]
split_path <- file.path(output_dir, "split_half.rds")
split_rows <- if (file.exists(split_path)) readRDS(split_path) else {
  out <- run(eligible, function(name) {
    sample <- samples[[name]]
    halves <- split_half(sample$counts, seed = 20261110L + match(name, names(samples)))
    fits <- lapply(halves, function(half) fit_one(half, sample$dt, birth_rate = central_rate))
    targets <- rownames(sample$counts)
    predictions <- lapply(fits, function(fit) stats::predict(fit, targets))
    total <- rowSums(sample$counts)
    # The halves have their own references: compare after removing the mean
    # difference over the karyotypes with at least 10 cells.
    difference <- predictions[[1L]]$relative_fitness - predictions[[2L]]$relative_fitness
    anchor <- total >= 10
    difference <- difference - mean(difference[anchor])
    compare <- function(use, label) {
      if (sum(use) < 3L) return(NULL)
      spread <- function(column) sqrt(predictions[[1L]][[column]][use]^2 + predictions[[2L]][[column]][use]^2)
      data.frame(sample = name, karyotypes = label, n = sum(use),
                 spearman = stats::cor(predictions[[1L]]$relative_fitness[use], predictions[[2L]]$relative_fitness[use],
                                       method = "spearman"),
                 rms_difference = sqrt(mean(difference[use]^2)),
                 within_reported = mean(abs(difference[use]) <= stats::qnorm(0.975) * spread("relative_sd")),
                 within_model = mean(abs(difference[use]) <= stats::qnorm(0.975) * spread("relative_sd_model")),
                 stringsAsFactors = FALSE)
    }
    rbind(compare(total >= 10, "10+ cells"), compare(total >= 2 & total < 10, "2-9 cells"),
          compare(total == 1, "1 cell"))
  })
  failed <- vapply(out, is.character, logical(1))
  if (any(failed)) cat("Split-half failures:", paste(eligible[failed], unlist(out[failed]), collapse = "; "), "\n")
  out <- do.call(rbind, out[!failed])
  saveRDS(out, split_path)
  out
}
utils::write.csv(split_rows, file.path(output_dir, "split_half.csv"), row.names = FALSE)

# 5. Cohort fit at the central rate.
cohort_table <- NULL
if (!isTRUE(option("skip-cohort", FALSE))) {
  cohort_path <- file.path(output_dir, "cohort_fit.rds")
  ok <- names(samples)[solo_table$status == "ok"]
  cohort <- if (file.exists(cohort_path)) readRDS(cohort_path) else {
    started <- proc.time()[["elapsed"]]
    fit <- fit_sparse_joint_cohort(lapply(samples[ok], `[[`, "counts"),
                                   dt = vapply(samples[ok], `[[`, numeric(1), "dt"), beta = beta,
                                   birth_rate = central_rate, min_cn = min_cn, max_cn = max_cn)
    fit$seconds <- proc.time()[["elapsed"]] - started
    saveRDS(fit, cohort_path)
    fit
  }
  cohort_table <- do.call(rbind, lapply(ok, function(name) {
    own <- cohort$samples[[name]]
    alone <- solo[[name]]$fits[[solo[[name]]$central]]$summary
    at <- match(own$karyotype, alone$karyotype)
    use <- !is.na(at) & (own$count_t0 + own$count_t1) > 0
    data.frame(sample = name, states = nrow(own), observed = sum(use),
               spearman_with_alone = if (sum(use) >= 3L) stats::cor(own$relative_fitness[use],
                                                                      alone$relative_fitness[at[use]],
                                                                      method = "spearman") else NA_real_,
               rms_shift = if (any(use)) sqrt(mean((own$relative_fitness[use] - alone$relative_fitness[at[use]])^2)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  utils::write.csv(cohort_table, file.path(output_dir, "cohort_against_alone.csv"), row.names = FALSE)
  recurrent <- cohort$landscape[order(-cohort$landscape$n_observed, -cohort$landscape$cells), ]
  utils::write.csv(recurrent, file.path(output_dir, "cohort_landscape.csv"), row.names = FALSE)
}

options(width = 200)
cat("\n== single-sample fits ==\n")
print(format(solo_table[-1L], digits = 3), row.names = FALSE)
ok <- solo_table$status == "ok"
cat(sprintf("\nfitted %d of %d; median seconds %.0f (max %.0f); median states %d (max %d)\n", sum(ok), nrow(solo_table),
            stats::median(solo_table$seconds[ok]), max(solo_table$seconds[ok]),
            as.integer(stats::median(solo_table$states[ok])), max(solo_table$states[ok])))
quartiles <- function(x) paste(signif(stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE), 2), collapse = " / ")
cat("rank agreement across birth rates (quartiles):", quartiles(solo_table$rank_agreement), "\n")
cat("deficit of non-major karyotypes, highest over lowest rate: per day", quartiles(solo_table$deficit_ratio),
    "; per division", quartiles(solo_table$per_division_ratio), "\n")
cat("dispersion of second-sample counts (1 = Poisson):", quartiles(solo_table$dispersion_t1), "\n")
if (!is.null(split_rows) && nrow(split_rows)) {
  cat("\n== split halves:", length(unique(split_rows$sample)), "samples with at least", split_min_cells,
      "cells at both timepoints ==\n")
  columns <- c("n", "spearman", "rms_difference", "within_reported", "within_model")
  print(format(stats::aggregate(split_rows[columns], split_rows["karyotypes"], stats::median), digits = 3), row.names = FALSE)
}
if (!is.null(cohort_table)) {
  d <- cohort$diagnostics
  cat(sprintf("\n== cohort fit: %d samples, %d karyotypes (%d shared), %d states, %.0f s, convergence %d ==\n",
              d$n_samples, d$n_nodes, d$shared_nodes, d$n_states, cohort$seconds, d$convergence))
  cat("scales:", paste(names(cohort$trend$scales), signif(cohort$trend$scales, 3), sep = "=", collapse = "  "), "\n")
  cat("agreement of a sample's cohort estimate with its own fit, Spearman (quartiles):",
      quartiles(cohort_table$spearman_with_alone), "\n")
  cat("\nmost recurrent karyotypes of the shared landscape (chromosome:copies where not 2):\n")
  top <- utils::head(recurrent, 12L)
  describe <- function(label) {
    cn <- as.integer(strsplit(label, ".", fixed = TRUE)[[1L]])
    changed <- which(cn != 2L)
    if (!length(changed)) return("all 2")
    paste(sprintf("%d:%d", changed, cn[changed]), collapse = " ")
  }
  print(data.frame(karyotype = vapply(top$karyotype, describe, character(1)), observed_in = top$n_observed,
                   major_in = top$n_major, cells = top$cells,
                   relative_fitness = signif(top$relative_fitness, 3), relative_sd = signif(top$relative_sd, 2)),
        row.names = FALSE)
}
