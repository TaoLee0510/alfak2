# Truth landscape of the sparse two-timepoint benchmark.
#
# The generator is the Gaussian-random-field landscape of the original ALFA-K
# agent-based simulations: fitness is a sum of sine waves of the Euclidean
# distance to random centroids, and is the per-day division rate of a karyotype.

stt_parse_karyotypes <- function(labels) {
  parts <- strsplit(as.character(labels), ".", fixed = TRUE)
  n_chr <- lengths(parts)
  if (length(unique(n_chr)) != 1L) stop("Karyotype labels have inconsistent dimensions.", call. = FALSE)
  matrix(as.numeric(unlist(parts, use.names = FALSE)), ncol = n_chr[[1L]], byrow = TRUE,
         dimnames = list(as.character(labels), NULL))
}

stt_format_karyotypes <- function(karyotypes) {
  karyotypes <- as.matrix(karyotypes)
  do.call(paste, c(as.data.frame(karyotypes), sep = "."))
}

# Centroids are redrawn until the founder divides faster than
# `min_founder_fitness`, as in the original simulations.
stt_generate_centroids <- function(founder, n_waves = 10L, wavelength, min_founder_fitness = 0.4) {
  repeat {
    centroids <- matrix(sample((-10):20, n_waves * length(founder), replace = TRUE),
                        nrow = n_waves, byrow = TRUE)
    if (stt_grf_fitness(matrix(founder, nrow = 1L), centroids, wavelength) > min_founder_fitness) {
      return(centroids)
    }
  }
}

# Division rate per day for each row of `karyotypes` (matrix or labels).
stt_grf_fitness <- function(karyotypes, centroids, wavelength) {
  if (is.character(karyotypes)) karyotypes <- stt_parse_karyotypes(karyotypes)
  karyotypes <- as.matrix(karyotypes)
  squared <- outer(rowSums(karyotypes^2), rowSums(centroids^2), "+") - 2 * karyotypes %*% t(centroids)
  rowSums(sin(sqrt(pmax(squared, 0)) / wavelength)) / (pi * sqrt(nrow(centroids)))
}

# The ABM is a tau-leap: each step a cell divides with probability
# fitness * delta_t, and non-positive fitness means no division. The mean
# per-day growth rate that the simulated counts actually follow is therefore
# log(1 + fitness * delta_t) / delta_t, not the nominal fitness. Scoring uses
# this rate so that no method is charged for the simulator's step bias.
stt_effective_growth_rate <- function(fitness, abm_delta_t) {
  log1p(pmax(fitness, 0) * abm_delta_t) / abm_delta_t
}

stt_true_growth_rate <- function(labels, truth) {
  stt_effective_growth_rate(
    stt_grf_fitness(labels, truth$centroids, truth$wavelength),
    truth$abm_delta_t
  )
}
