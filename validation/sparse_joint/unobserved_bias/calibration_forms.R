#!/usr/bin/env Rscript
# Forms of the reported uncertainty, fitted to the pooled default fits of
# calibration.R (results/calibration_nodes.rds) and compared by coverage.
#
#   shipped   floor 0.009, deficit 0.6 (the package's constants)
#   refit     the shipped form, its two constants refitted:
#               sd^2 = model^2 + (deficit * relative)^2 + (floor * (1 + shell))^2
#   split     a deficit for observed and one for unobserved karyotypes:
#               sd^2 = model^2 + (deficit_k * relative)^2 + (floor * (1 + shell))^2
#
# The baseline birth rate is 1 in every fit here, so the floor's factor
# birth_rate is left out. Constants are fitted by Gaussian likelihood with
# every dataset, and within it every kind of karyotype, weighted equally, on
# the sets named by --fit (default: dev core, diagnosis, rugged; the stress
# scenarios, where the landscape changes, are reported only, as before).
# --leave-out=<set> refits without one set and reports it.
#
#   Rscript validation/sparse_joint/unobserved_bias/calibration_forms.R [--leave-out=dev_core|diagnosis|rugged]

if (!exists("ub_dir")) ub_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})

ub_load_calibration_nodes <- function(dir = ub_dir) {
  x <- readRDS(file.path(dir, "results", "calibration_nodes.rds"))
  x$error <- x$pred - x$truth
  x$set <- ifelse(x$source == "dev", paste0("dev_", x$role), x$source)
  x$observed <- x$distance == 0L
  cell <- paste(x$dataset_id, x$kind)
  x$weight <- 1 / as.numeric(table(cell)[cell])
  x
}

ub_forms <- list(
  shipped = list(start = numeric(0), sd = function(p, x) {
    sqrt(x$model_sd^2 + (0.6 * x$pred)^2 + (0.009 * (1 + x$shell))^2)
  }, names = character(0)),
  refit = list(start = log(c(0.009, 0.6)), sd = function(p, x) {
    sqrt(x$model_sd^2 + (p[[2L]] * x$pred)^2 + (p[[1L]] * (1 + x$shell))^2)
  }, names = c("floor", "deficit")),
  split = list(start = log(c(0.009, 0.6, 0.6)), sd = function(p, x) {
    deficit <- ifelse(x$observed, p[[2L]], p[[3L]])
    sqrt(x$model_sd^2 + (deficit * x$pred)^2 + (p[[1L]] * (1 + x$shell))^2)
  }, names = c("floor", "deficit_observed", "deficit_unobserved"))
)

ub_fit_form <- function(form, x) {
  if (!length(form$start)) return(numeric(0))
  objective <- function(log_p) {
    v <- form$sd(exp(log_p), x)^2
    sum(x$weight * (log(v) + x$error^2 / v))
  }
  best <- stats::optim(form$start, objective, control = list(maxit = 4000, reltol = 1e-12))
  stats::setNames(exp(best$par), form$names)
}

# Dataset-weighted coverage and mean reported sd, by set, wavelength and kind.
ub_coverage <- function(form, p, x) {
  x$sd <- form$sd(p, x)
  x$covered <- abs(x$error) <= stats::qnorm(0.975) * x$sd
  per <- stats::aggregate(cbind(covered, sd) ~ set + wavelength_class + kind + dataset_id, x, mean)
  stats::aggregate(cbind(covered, sd) ~ set + wavelength_class + kind, per, mean)
}

ub_wide <- function(table, value) {
  w <- stats::reshape(table[c("set", "wavelength_class", "kind", value)], idvar = c("set", "wavelength_class"),
                      timevar = "kind", direction = "wide")
  names(w) <- sub(paste0(value, "."), "", names(w), fixed = TRUE)
  w[order(w$set, w$wavelength_class), ]
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  leave_out <- sub("^--leave-out=", "", c(grep("^--leave-out=", args, value = TRUE), "--leave-out=")[[1L]])
  x <- ub_load_calibration_nodes()
  fit_sets <- setdiff(c("dev_core", "diagnosis", "rugged"), leave_out)
  train <- x[x$set %in% fit_sets, ]
  options(width = 220)
  for (name in names(ub_forms)) {
    form <- ub_forms[[name]]
    p <- ub_fit_form(form, train)
    cat(sprintf("\n==== %s %s (fitted on %s)\n", name,
                paste(sprintf("%s=%.4f", names(p), p), collapse = " "), paste(fit_sets, collapse = ", ")))
    table <- ub_coverage(form, p, x)
    table$covered <- round(table$covered, 3)
    print(ub_wide(table, "covered"), row.names = FALSE)
  }
}
