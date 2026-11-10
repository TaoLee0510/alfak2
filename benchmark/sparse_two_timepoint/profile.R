# Observable summary of a two-timepoint count matrix.
#
# These statistics need no truth, so the same function describes simulated
# datasets and real inputs. They are used to check that the benchmark's
# sampling window produces data that look like the real two-timepoint samples
# in diversity and in turnover between the timepoints.

stt_count_profile <- function(counts) {
  counts <- as.matrix(counts)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  total <- rowSums(counts)
  p0 <- counts[, 1L] / sum(counts[, 1L])
  p1 <- counts[, 2L] / sum(counts[, 2L])
  shared <- counts[, 1L] > 0 & counts[, 2L] > 0
  # Spread of log frequency change among well-sampled shared karyotypes: the
  # directly visible amount of selection over the interval.
  well_sampled <- shared & total >= 20
  c(
    depth_t0 = sum(counts[, 1L]),
    depth_t1 = sum(counts[, 2L]),
    n_karyotypes = nrow(counts),
    shared_fraction = mean(shared),
    shared_mass_t0 = sum(p0[shared]),
    shared_mass_t1 = sum(p1[shared]),
    total_variation = 0.5 * sum(abs(p1 - p0)),
    top_fraction = max(total) / sum(total),
    singleton_fraction = mean(total == 1),
    n_count_ge5 = sum(total >= 5),
    log_ratio_sd = if (sum(well_sampled) >= 3L) stats::sd(log(p1[well_sampled] / p0[well_sampled])) else NA_real_
  )
}

# Profiles of every ALFA-K style input (`list(x = counts, dt = ...)`) in a
# directory. The files are read in place and never copied.
stt_profile_directory <- function(directory) {
  files <- list.files(directory, pattern = "\\.[Rr]ds$", full.names = TRUE)
  rows <- lapply(files, function(file) {
    input <- readRDS(file)
    counts <- as.matrix(if (is.list(input) && !is.null(input$x)) input$x else input)
    days <- suppressWarnings(as.numeric(colnames(counts)))
    data.frame(sample = sub("\\.[Rr]ds$", "", basename(file)),
               interval = if (length(days) == 2L && all(is.finite(days))) diff(days) else NA_real_,
               t(stt_count_profile(counts)), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

stt_profile_quantiles <- function(profiles, columns = c("n_karyotypes", "shared_fraction", "total_variation",
                                                         "top_fraction", "singleton_fraction", "log_ratio_sd")) {
  out <- t(vapply(columns, function(column) {
    stats::quantile(profiles[[column]], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE, names = FALSE)
  }, numeric(5)))
  colnames(out) <- c("q10", "q25", "median", "q75", "q90")
  data.frame(statistic = columns, out, row.names = NULL, stringsAsFactors = FALSE)
}
