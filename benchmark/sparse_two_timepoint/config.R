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
  # `shift` is the fraction of the landscape's waves redrawn at the first
  # sample. Core scenarios make the acceptance verdict; stress scenarios are
  # reported beside it.
  scenario <- function(name, depth_t0, depth_t1, interval, bottleneck = 0, shift = 0, role = "core") {
    data.frame(scenario = name, depth_t0 = depth_t0, depth_t1 = depth_t1,
               interval = interval, bottleneck = bottleneck, shift = shift, role = role,
               stringsAsFactors = FALSE)
  }
  landscapes <- rbind(make_split("dev", 71000L), make_split("heldout", 82000L))
  # What a cohort method needs to find the other samples of its cohort.
  cohort_method <- function(members, ...) {
    list(type = "alfak2_sparse_joint_cohort", birth_rate = 1, min_cn = 1L, max_cn = 8L, members = members,
         dataset_dir = file.path(output_dir, "datasets"), landscapes = landscapes[c("landscape_id", "split")],
         cohort_size = 8L, arguments = utils::modifyList(list(initial_state = "free"), list(...)))
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
    landscapes = landscapes,
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
      scenario("bottleneck_500_700_360d", 500, 700, 360, bottleneck = 200),
      # Few generations between the samples: 27 and 54 at the simulated
      # division rate, as for a slowly dividing tumour sampled months apart.
      scenario("short_500_700_30d", 500, 700, 30, role = "stress"),
      scenario("short_500_700_60d", 500, 700, 60, role = "stress"),
      # The landscape changes at the first sample, as when treatment starts
      # then: half of its waves are redrawn, which leaves fitness before and
      # after correlated. Truth is the landscape after.
      scenario("shift_500_700_360d", 500, 700, 360, shift = 0.5, role = "stress"),
      scenario("shift_3000_100_360d", 3000, 100, 360, shift = 0.5, role = "stress"),
      # Worst case for anything read from the first sample: the landscape is
      # turned about the dominant karyotype, so that which of its neighbours
      # is fit before says nothing about after.
      scenario("turned_500_700_360d", 500, 700, 360, shift = 1, role = "stress"),
      scenario("turned_3000_100_360d", 3000, 100, 360, shift = 1, role = "stress")
    ),
    cohort_scenario = "typical_500_700_360d",
    # Scored nodes: all observed karyotypes, and up to `d1_cap` / `d2_cap`
    # uniformly sampled one- and two-step neighbours.
    scoring = list(min_cn = 1L, max_cn = 8L, d1_cap = 2000L, d2_cap = 2000L, reference_min_count = 5),
    # Roles in the acceptance rule (report.R).
    acceptance = list(zero = "zero", current = "alfak2_two_layer_fixed",
                      original = c("alfakR_minobs5", "alfakR_minobs20"),
                      candidates = c("alfak2_sparse_joint", "alfak2_sparse_joint_balance"),
                      # Cohort methods exist for the cohort scenario only and are
                      # compared there with the single-sample model.
                      cohort = c("alfak2_sparse_joint_cohort", "alfak2_sparse_joint_cohort_unrelated",
                                 "alfak2_sparse_joint_cohort_shared_trend",
                                 "alfak2_sparse_joint_cohort_unrelated_shared_trend",
                                 "alfak2_sparse_joint_cohort_balance"),
                      cohort_baseline = "alfak2_sparse_joint"),
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
      alfak2_two_layer_cv = list(type = "alfak2_two_layer", local_shell_depth = 0L, max_nodes = 60000L),
      # Candidate. One division per day is assumed for every dataset; the
      # simulated populations divide about 0.9 times per day.
      # `initial_state = "free"` was the package default when these were
      # benchmarked; "balance" is the default now and has its own entry.
      alfak2_sparse_joint = list(type = "alfak2_sparse_joint", birth_rate = 1, min_cn = 1L, max_cn = 8L,
                                 arguments = list(initial_state = "free")),
      # Sensitivity to the assumed baseline birth rate: half and double.
      alfak2_sparse_joint_slow = list(type = "alfak2_sparse_joint", birth_rate = 0.5, min_cn = 1L, max_cn = 8L,
                                      arguments = list(initial_state = "free")),
      alfak2_sparse_joint_fast = list(type = "alfak2_sparse_joint", birth_rate = 2, min_cn = 1L, max_cn = 8L,
                                      arguments = list(initial_state = "free")),
      # Initial abundances of non-major karyotypes tied to mutation-selection
      # balance, which reads the first sample's levels as well.
      alfak2_sparse_joint_balance = list(type = "alfak2_sparse_joint", birth_rate = 1, min_cn = 1L, max_cn = 8L,
                                         arguments = list(initial_state = "balance")),
      # Baseline birth rate known only within a factor of 16. The simulated
      # rate, about 0.9, is near the centre, the top or the bottom of the range.
      alfak2_sparse_joint_range = list(type = "alfak2_sparse_joint_rates", birth_rate_range = c(0.25, 4),
                                       min_cn = 1L, max_cn = 8L, arguments = list(initial_state = "free")),
      alfak2_sparse_joint_range_low = list(type = "alfak2_sparse_joint_rates", birth_rate_range = c(0.0625, 1),
                                           min_cn = 1L, max_cn = 8L, arguments = list(initial_state = "free")),
      alfak2_sparse_joint_range_high = list(type = "alfak2_sparse_joint_rates", birth_rate_range = c(1, 16),
                                            min_cn = 1L, max_cn = 8L, arguments = list(initial_state = "free")),
      # Cohort fit of the eight patients of a landscape, scored on patient 1.
      # The trend is shared or per sample, whichever has the higher marginal
      # likelihood.
      alfak2_sparse_joint_cohort = cohort_method("landscape", trend = "select"),
      # Control: the same patient with patients of seven other landscapes,
      # which share the statistics of a landscape but not the landscape.
      alfak2_sparse_joint_cohort_unrelated = cohort_method("other_landscapes", trend = "select"),
      # The same two with one trend for the whole cohort, always.
      alfak2_sparse_joint_cohort_shared_trend = cohort_method("landscape", trend = "shared"),
      alfak2_sparse_joint_cohort_unrelated_shared_trend = cohort_method("other_landscapes", trend = "shared"),
      # The cohort fit with initial abundances at mutation-selection balance.
      alfak2_sparse_joint_cohort_balance = cohort_method("landscape", trend = "select", initial_state = "balance")
    )
  )
}
