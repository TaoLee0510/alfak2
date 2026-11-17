# Checks of the benchmark simulator against quantities known in closed form.
#
# The daughter distribution is compared with the single-chromosome
# missegregation kernel of alfak2 (`alfak2_pij_cpp`), which is an independent
# implementation of the same biology.

stt_validate_abm <- function(config) {
  checks <- list()
  record <- function(name, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = isTRUE(passed), detail = detail,
                                                 stringsAsFactors = FALSE)
  }
  set.seed(90210)
  n_chr <- 4L
  wavelength <- 1.5
  centroids <- stt_generate_centroids(rep(2, n_chr), 10L, wavelength)
  parent <- c(2L, 3L, 1L, 2L)
  parent_label <- paste(parent, collapse = ".")
  fitness <- stt_grf_fitness(matrix(parent, nrow = 1L), centroids, wavelength)
  one_step <- function(beta, min_cn, max_cn, n_cells, seed, delta_t) {
    local_config <- config
    local_config$abm <- utils::modifyList(config$abm, list(delta_t = delta_t, max_pop = 0,
                                                           min_cn = min_cn, max_cn = max_cn))
    stt_run_abm(stats::setNames(n_cells, parent_label), centroids, wavelength, local_config,
                record_days = delta_t, seed = seed, beta = beta)[[1L]]
  }
  kernel <- function(child, beta) {
    prod(mapply(function(from, to) alfak2:::alfak2_pij_cpp(from, to, beta), parent, child))
  }

  # 1. Daughter distribution without viability losses.
  beta <- 0.01
  n_cells <- 2e7
  delta_t <- 0.5 / fitness
  after <- one_step(beta, min_cn = 0L, max_cn = 20L, n_cells, seed = 11, delta_t)
  n_divisions <- sum(after) - n_cells
  children <- setdiff(names(after), parent_label)
  expected <- 2 * n_divisions * vapply(children, function(label) kernel(stt_parse_karyotypes(label)[1L, ], beta), numeric(1))
  large <- expected >= 50
  z <- (after[children][large] - expected[large]) / sqrt(expected[large])
  record("division probability is fitness * delta_t",
         abs(n_divisions / n_cells - 0.5) < 5 * sqrt(0.25 / n_cells),
         sprintf("dividing fraction %.5f, expected 0.5", n_divisions / n_cells))
  record("daughter karyotypes follow the missegregation kernel",
         max(abs(z)) < 4.5 && abs(mean(z)) < 4.5 / sqrt(length(z)),
         sprintf("%d daughter types with expected count >= 50; max |z| = %.2f, mean z = %.2f",
                 sum(large), max(abs(z)), mean(z)))
  expected_parent <- (n_cells - n_divisions) + 2 * n_divisions * kernel(parent, beta)
  record("faithful daughters and non-dividing parents are conserved",
         abs(after[[parent_label]] - expected_parent) < 5 * sqrt(2 * n_divisions * (1 - kernel(parent, beta))),
         sprintf("parent count %.0f, expected %.0f", after[[parent_label]], expected_parent))

  # 2. Viability bounds: daughters outside [min_cn, max_cn] die.
  bounded <- one_step(beta, min_cn = 1L, max_cn = 3L, n_cells, seed = 12, delta_t)
  bounded_karyotypes <- stt_parse_karyotypes(names(bounded))
  record("no karyotype outside the viability bounds survives",
         all(bounded_karyotypes >= 1 & bounded_karyotypes <= 3),
         sprintf("copy-number range %d..%d under bounds 1..3", min(bounded_karyotypes), max(bounded_karyotypes)))
  viable_probability <- sum(vapply(names(bounded), function(label) kernel(stt_parse_karyotypes(label)[1L, ], beta), numeric(1)))
  # Total after one step = non-dividers + viable daughters, so the dividing
  # count solves n_cells - D + 2 * D * viable_probability = total (in mean).
  implied_divisions <- (sum(bounded) - n_cells) / (2 * viable_probability - 1)
  record("mass lost equals the kernel probability outside the bounds",
         abs(implied_divisions / n_cells - 0.5) < 2e-3,
         sprintf("retained single-daughter probability %.5f; implied dividing fraction %.5f",
                 viable_probability, implied_divisions / n_cells))

  # 3. Mean growth factor per step is exactly 1 + fitness * delta_t.
  local_config <- config
  local_config$abm <- utils::modifyList(config$abm, list(delta_t = 0.1, max_pop = 0, min_cn = 0L, max_cn = 20L))
  grown <- stt_run_abm(stats::setNames(1e6, parent_label), centroids, wavelength, local_config,
                       record_days = 5, seed = 13, beta = 0)[[1L]]
  expected_growth <- (1 + fitness * 0.1)^50
  record("growth without missegregation matches (1 + fitness * delta_t)^steps",
         length(grown) == 1L && abs(sum(grown) / 1e6 / expected_growth - 1) < 5e-3,
         sprintf("growth factor %.4f, expected %.4f", sum(grown) / 1e6, expected_growth))

  # 4. Selection between two clones follows the effective growth rates.
  other <- c(3L, 3L, 1L, 2L)
  other_label <- paste(other, collapse = ".")
  other_fitness <- stt_grf_fitness(matrix(other, nrow = 1L), centroids, wavelength)
  raced <- stt_run_abm(stats::setNames(c(5e5, 5e5), c(parent_label, other_label)), centroids, wavelength,
                       local_config, record_days = 10, seed = 14, beta = 0)[[1L]]
  observed <- log(raced[[other_label]] / raced[[parent_label]]) / 10
  expected_rate <- diff(stt_effective_growth_rate(c(fitness, other_fitness), 0.1))
  record("log frequency ratio changes at the effective growth-rate difference",
         abs(observed - expected_rate) < 1e-3,
         sprintf("observed %.5f per day, expected %.5f", observed, expected_rate))

  # 5. Passaging keeps the population at or below the cap.
  local_config$abm <- utils::modifyList(config$abm, list(delta_t = 0.1, max_pop = 2e5, culling_survival = 0.01,
                                                         min_cn = 1L, max_cn = 8L))
  passaged <- stt_run_abm(stats::setNames(5e4, parent_label), centroids, wavelength, local_config,
                          record_days = seq(10, 200, by = 10), seed = 15, beta = 5e-5)
  totals <- vapply(passaged, sum, numeric(1))
  record("population stays within the passaging cap",
         all(totals <= 2e5) && min(totals) >= 1e3,
         sprintf("population range %.0f..%.0f with cap 2e5 and survival 0.01", min(totals), max(totals)))

  # 6. Same seed, same trajectory.
  again <- stt_run_abm(stats::setNames(5e4, parent_label), centroids, wavelength, local_config,
                       record_days = seq(10, 200, by = 10), seed = 15, beta = 5e-5)
  record("simulation is reproducible for a fixed seed", identical(passaged, again), "identical snapshots")

  do.call(rbind, checks)
}
