# Design of the sparse two-timepoint benchmark.
#
# The simulated data are meant to look like real two-timepoint tumour samples.
# The reference is a set of 20 GBM inputs: 22 chromosomes, copy numbers 1-4,
# 8-3889 cells per timepoint with often strongly unbalanced depths, intervals
# of 180-3600 days, typically ~65 karyotypes of which a third are seen at both
# timepoints, a largest clone of ~13%, and a total-variation distance of ~0.8
# between the two samples.
#
# Three choices follow from matching those observables (README, "Design"):
# - 22 chromosomes. With 8 chromosomes the simulated population is trapped on a
#   local peak within ~700 days and samples become 98% one karyotype.
# - First sample at day 1440. Earlier, clones replace each other completely
#   within a year; later, the smoother landscapes also freeze.
# - Wavelengths 0.4, 0.8, 1.2 and 1.6. The 0.4 class is the frozen,
#   near-monoclonal regime that about one real sample in ten resembles; the
#   others give ongoing adaptation with real-like diversity and turnover.

stt_config <- function(output_dir = file.path("benchmark", "results", "sparse_two_timepoint")) {
  wavelength_classes <- data.frame(
    wavelength_class = c("w0.4", "w0.8", "w1.2", "w1.6"),
    wavelength = c(0.4, 0.8, 1.2, 1.6),
    stringsAsFactors = FALSE
  )
  landscapes_per_class <- 3L
  make_split <- function(split, seed_base) {
    grid <- wavelength_classes[rep(seq_len(nrow(wavelength_classes)), each = landscapes_per_class), ]
    index <- seq_len(nrow(grid))
    data.frame(
      landscape_id = sprintf("%s_L%02d", split, index),
      split = split,
      wavelength_class = grid$wavelength_class,
      wavelength = grid$wavelength,
      landscape_seed = seed_base + index,
      stringsAsFactors = FALSE
    )
  }
  scenario <- function(name, depth_t0, depth_t1, interval, bottleneck = 0) {
    data.frame(scenario = name, depth_t0 = depth_t0, depth_t1 = depth_t1,
               interval = interval, bottleneck = bottleneck, stringsAsFactors = FALSE)
  }
  list(
    output_dir = output_dir,
    n_chr = 22L,
    founder_cn = 2L,
    n_waves = 10L,
    beta = 5e-5,
    # Karyotypes with a copy number outside [min_cn, max_cn] are not viable.
    abm = list(pop_size = 5e4, max_pop = 2e6, delta_t = 0.1, culling_survival = 0.01,
               min_cn = 1L, max_cn = 8L, record_days = 60, end_day = 2520),
    t0 = 1440,
    # Development landscapes are used freely. Held-out landscapes are generated
    # and scored only once, for the final confirmation.
    landscapes = rbind(make_split("dev", 71000L), make_split("heldout", 82000L)),
    patients_per_landscape = 8L,
    # Patients fitted by every method for the single-sample scoreboard. One per
    # landscape keeps a full development round within a few CPU-hours on a
    # workstation; landscapes, not patients, are the independent unit.
    scoreboard_patients = 1L,
    scenarios = rbind(
      scenario("very_sparse_50_50_180d", 50, 50, 180),
      scenario("sparse_200_200_180d", 200, 200, 180),
      scenario("sparse_200_200_360d", 200, 200, 360),
      scenario("typical_500_700_360d", 500, 700, 360),
      scenario("unbalanced_100_3000_360d", 100, 3000, 360),
      scenario("unbalanced_3000_100_360d", 3000, 100, 360),
      scenario("long_500_700_1080d", 500, 700, 1080),
      # After the first sample only 200 cells survive and regrow.
      scenario("bottleneck_500_700_360d", 500, 700, 360, bottleneck = 200)
    ),
    cohort_scenario = "typical_500_700_360d",
    # Scored nodes: all observed karyotypes, and up to `d1_cap` / `d2_cap`
    # uniformly sampled one- and two-step neighbours.
    scoring = list(min_cn = 1L, max_cn = 8L, d1_cap = 2000L, d2_cap = 2000L, reference_min_count = 5),
    # Roles in the acceptance rule (report.R). Any other method is a candidate.
    acceptance = list(zero = "zero", current = "alfak2_two_layer_fixed",
                      original = c("alfakR_minobs5", "alfakR_minobs20")),
    methods = list(
      zero = list(type = "zero"),
      alfakR_minobs5 = list(type = "alfakR", minobs = 5L, nboot = 45L, n0 = 5e4, nb = 2e6,
                            nn_prior = "empirical"),
      alfakR_minobs20 = list(type = "alfakR", minobs = 20L, nboot = 45L, n0 = 5e4, nb = 2e6,
                             nn_prior = "empirical"),
      # Graph hyperparameters fixed at the centre of the package's default grid.
      alfak2_two_layer_fixed = list(type = "alfak2_two_layer", local_shell_depth = 0L, max_nodes = 60000L,
                                    graph_hyperparameters = c(lambda_l = 1, lambda_e = 0.25, sigma_obs = 0.05)),
      # Package default: hyperparameters by leave-one-anchor-out selection.
      alfak2_two_layer_cv = list(type = "alfak2_two_layer", local_shell_depth = 0L, max_nodes = 60000L)
    )
  )
}
