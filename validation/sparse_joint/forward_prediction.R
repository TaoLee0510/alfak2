# Does a fit predict what the population looks like later?
#
#   Rscript validation/sparse_joint/forward_prediction.R [--split=dev] [--horizon=180] [--workers=2]
#                                                        [--unseal-heldout]
#
# For every scoreboard dataset whose simulated trajectory was recorded
# `horizon` days after the second sample, the sparse joint model is fitted to
# the two samples and its domain is run forward from the fitted composition
# at the second sample. The predicted composition is compared, as
# total-variation distance over the domain, with the true population at that
# later day. Four forecasts are compared:
#   no change        the fitted composition at the second sample, unchanged;
#   no selection     run forward with every karyotype at the baseline rate;
#   fitted           run forward with the fitted birth rates;
#   true fitness     run forward with the true rates: the best this forward
#                    model and this starting point can do.
# Set STT_R_LIB to the library holding the alfak2 build under test.

library_dir <- Sys.getenv("STT_R_LIB", unset = "")
if (nzchar(library_dir)) .libPaths(c(library_dir, .libPaths()))
for (part in c("landscape.R", "config.R", "generate.R", "profile.R", "score.R", "methods.R", "report.R")) {
  source(file.path("benchmark/sparse_two_timepoint", part))
}
args <- commandArgs(TRUE)
option <- function(name, default) {
  hit <- grep(paste0("^--", name, "(=|$)"), args, value = TRUE)
  if (!length(hit)) return(default)
  if (grepl("=", hit[[1L]], fixed = TRUE)) sub("^[^=]*=", "", hit[[1L]]) else TRUE
}
split <- option("split", "dev")
if (identical(split, "heldout") && !isTRUE(option("unseal-heldout", FALSE))) {
  stop("Held-out landscapes are reserved for the final confirmation; pass --unseal-heldout.", call. = FALSE)
}
horizon <- as.numeric(option("horizon", "180"))
workers <- as.integer(option("workers", "2"))
config <- stt_config()
method <- config$methods$alfak2_sparse_joint
index <- utils::read.csv(file.path(config$output_dir, sprintf("dataset_index_%s.csv", split)), stringsAsFactors = FALSE)
index <- index[index$scoreboard & index$role == "core", ]
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

total_variation <- function(p, q) 0.5 * sum(abs(p / sum(p) - q / sum(q)))
forecast <- function(dataset_id) {
  dataset <- readRDS(stt_dataset_path(config, dataset_id))
  trajectory <- readRDS(stt_truth_path(config, dataset$landscape_id, dataset$patient))
  later <- trajectory$snapshots[[as.character(dataset$t0 + dataset$interval + horizon)]]
  # The bottleneck scenario leaves the recorded trajectory at the first sample.
  if (is.null(later) || dataset$bottleneck > 0) return(NULL)
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  fit <- alfak2::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta,
                                  birth_rate = method$birth_rate, min_cn = method$min_cn, max_cn = method$max_cn)
  summary <- fit$summary
  karyotypes <- alfak2:::parse_karyotypes(summary$karyotype)
  kernel <- alfak2:::sparse_division_kernel(karyotypes, dataset$beta)
  self <- 2 * exp(kernel$log_self) - 1
  start <- log(summary$expected_t1)
  run <- function(log_birth) {
    exp(alfak2:::closed_form_forward_r(start, log_birth, kernel, horizon, order = 1L) - max(start))
  }
  true_rate <- stt_true_growth_rate(summary$karyotype, dataset$truth)
  truth_later <- later[summary$karyotype]
  truth_later[is.na(truth_later)] <- 0
  truth_now <- dataset$truth$frequency_t1[summary$karyotype]
  truth_now[is.na(truth_now)] <- 0
  forecasts <- list(
    no_change = exp(start),
    no_selection = run(rep(log(method$birth_rate), nrow(summary))),
    fitted = run(summary$log_birth_mean),
    # A karyotype that does not divide in the simulation gets a negligible rate.
    true_fitness = run(log(pmax(true_rate, 1e-6) / self))
  )
  data.frame(
    dataset_id = dataset_id, scenario = dataset$scenario, landscape_id = dataset$landscape_id,
    forecast = names(forecasts),
    distance = vapply(forecasts, total_variation, numeric(1), q = truth_later),
    # how far the true population itself moved inside the domain, and how much of it the domain holds
    true_change = total_variation(truth_now, truth_later),
    domain_mass_later = sum(truth_later) / sum(later),
    stringsAsFactors = FALSE, row.names = NULL
  )
}
rows <- parallel::mclapply(index$dataset_id, function(id) tryCatch(forecast(id), error = function(e) conditionMessage(e)),
                           mc.cores = workers, mc.preschedule = FALSE)
failed <- vapply(rows, is.character, logical(1))
if (any(failed)) cat("Failed:", sum(failed), "; first:", rows[[which(failed)[1L]]], "\n")
rows <- do.call(rbind, rows[!failed])
saveRDS(rows, file.path(config$output_dir, "report", sprintf("forward_prediction_%s.rds", split)))

# Landscapes weigh equally within a scenario.
per_landscape <- stats::aggregate(rows[c("distance", "true_change", "domain_mass_later")],
                                  rows[c("scenario", "landscape_id", "forecast")], mean)
table <- stats::aggregate(per_landscape[c("distance", "true_change", "domain_mass_later")],
                          per_landscape[c("scenario", "forecast")], mean)
wide <- stats::reshape(table[c("scenario", "forecast", "distance")], idvar = "scenario", timevar = "forecast",
                       direction = "wide")
names(wide) <- sub("distance.", "", names(wide), fixed = TRUE)
context <- table[table$forecast == "fitted", c("scenario", "true_change", "domain_mass_later")]
wide <- merge(wide, context)
options(width = 200)
cat(sprintf("Total-variation distance to the true composition %g days after the second sample (%s landscapes, %d datasets)\n",
            horizon, split, length(unique(rows$dataset_id))))
print(format(wide[c("scenario", "no_change", "no_selection", "fitted", "true_fitness", "true_change", "domain_mass_later")],
             digits = 3), row.names = FALSE)
pooled <- stats::aggregate(per_landscape["distance"], per_landscape[c("landscape_id", "forecast")], mean)
overall <- tapply(pooled$distance, pooled$forecast, mean)
cat("\npooled over scenarios:", paste(names(overall), signif(overall, 3), sep = " = ", collapse = "   "), "\n")
paired <- stats::reshape(pooled, idvar = "landscape_id", timevar = "forecast", direction = "wide")
difference <- paired$distance.fitted - paired$distance.no_change
set.seed(20261110L)
boot <- replicate(2000L, mean(sample(difference, replace = TRUE)))
cat(sprintf("fitted minus no change: %.3f (landscape bootstrap 95%% interval %.3f to %.3f)\n", mean(difference),
            stats::quantile(boot, 0.025), stats::quantile(boot, 0.975)))
