# Few frequent karyotypes on rugged landscapes: the second weakness of
# fit_sparse_joint() and change C, a residual correlated along the karyotype
# lattice (README.md in this directory).
#
# The data come from the development simulator of
# benchmark/sparse_two_timepoint at the population sizes of the ALFA-K paper's
# ABM. The stage code and scores are those of the short-interval study
# (../short_interval/setup.R and stages.R), sourced unchanged; this study has
# its own landscapes, seeds and output directory.

lr_dir <- if (exists("lr_script_dir")) lr_script_dir else getwd()
si_script_dir <- normalizePath(file.path(lr_dir, "..", "short_interval"))
source(file.path(si_script_dir, "setup.R"))
source(file.path(si_script_dir, "stages.R"))

lr_default_output <- file.path(si_repo, "validation", "sparse_joint", "lattice_residual", "results")

# `set = "heldout"`: the same design on landscapes no one has looked at (seed
# base 96000, ids "lh_L.."), in results/heldout, for the final confirmation
# only.
lr_config <- function(output_dir = NULL, set = c("rugged", "heldout")) {
  set <- match.arg(set)
  if (is.null(output_dir)) output_dir <- file.path(lr_default_output, set)
  config <- si_config(output_dir)
  config$set <- set
  config$library_root <- lr_default_output
  # Wavelength 0.2 is the regime of the weakness; 0.4 and 0.8 are where a
  # correlated residual has neighbours that are alike. Eight landscapes per
  # class, and four more at wavelength 0.2 (numbered after the others): there
  # the number of frequent karyotypes varies from one to about ten between
  # landscapes, and half of the first eight had only one or two.
  classes <- data.frame(wavelength_class = c("w0.2", "w0.4", "w0.8"), wavelength = c(0.2, 0.4, 0.8),
                        stringsAsFactors = FALSE)
  grid <- classes[c(rep(1:3, each = 8L), rep(1L, 4L)), ]
  index <- seq_len(nrow(grid))
  # Seed bases 71000 and 82000 belong to the development benchmark, 93000 and
  # 94000 to the short-interval study.
  config$landscapes <- data.frame(
    landscape_id = sprintf(if (set == "rugged") "lr_L%02d" else "lh_L%02d", index), split = "lr",
    wavelength_class = grid$wavelength_class, wavelength = grid$wavelength,
    landscape_seed = (if (set == "rugged") 95000L else 96000L) + index, stringsAsFactors = FALSE
  )
  # The paper's ABM grows each passage from 2e5 to 2e7 cells; at the
  # development benchmark's 5e4 to 2e6, wavelength 0.2 gives a median of two
  # frequent karyotypes instead of the paper's three to four.
  config$abm$pop_size <- 2e5
  config$abm$max_pop <- 2e7
  config$abm$end_day <- 120
  config$abm$record_days <- 60
  # First sample at day 120, about the age of the paper's populations at their
  # windows; intervals of one, three and seven passages.
  intervals <- c(6, 18, 42)
  config$scenarios <- data.frame(
    scenario = sprintf("t0120_%02dd", intervals), depth_t0 = 1000, depth_t1 = 1000, interval = intervals,
    t0 = 120, bottleneck = 0, shift = 0, role = "diagnosis", stringsAsFactors = FALSE
  )
  sparse <- function(build, ...) {
    list(type = "sparse_joint", build = build, birth_rate = 1, min_cn = 1L, max_cn = 8L,
         arguments = list(...))
  }
  config$methods <- list(
    zero = list(type = "zero"),
    # The package as committed (af362af), at its defaults.
    base = sparse("af362af")
  )
  config
}
