#!/usr/bin/env Rscript
repo <- normalizePath(".")
directory <- file.path(repo, "validation/baseline_quadrature")
reference <- file.path(repo, "validation/posterior_reference_4to8/results")
source(file.path(repo, "validation/posterior_reference_4to8/reference.R"))
source(file.path(repo, "validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo, "validation/accuracy_round2/common.R"))
source(file.path(directory, "quadrature.R"))
runtime <- reference_load(repo, compile = FALSE)
cases <- readRDS(file.path(reference, "frozen_cases.rds"))
output <- file.path(directory, "results")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
files <- c(file.path(directory, c("PROTOCOL.md", "quadrature.R", "run.R")),
  list.files(file.path(repo, "R"), "[.]R$", full.names = TRUE), runtime$native,
  file.path(reference, "frozen_cases.rds"))
hash <- round2_sha256(files)
manifest <- data.frame(file = names(hash), sha256 = unname(hash))
freeze <- file.path(output, "source_sha256.csv")
if (file.exists(freeze) && !identical(read.csv(freeze), manifest)) stop("Frozen candidate source changed")
if (!file.exists(freeze)) write.csv(manifest, freeze, row.names = FALSE)
key <- round2_hash_object(manifest)
cache <- file.path("/private/tmp/alfak2_baseline_quad_20261104", substr(key, 1, 12))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)

# Independent prior-law check: mixture over the removed intercept must retain
# the complete original log-rate and physical birth prior moments.
for (n in c(4L, 8L)) {
  case <- cases[[which(vapply(cases, function(x) x$n == n, logical(1)))[1L]]]
  j <- which(colnames(case$prior$design) == "intercept"); sd0 <- case$prior$coef_sd[j]
  Krest <- case$covariance - matrix(sd0^2, n, n)
  rule <- oracle_normal_quadrature(17L)
  pieces <- lapply(sd0 * rule$nodes, function(a) {
    m <- case$prior$offset + a; b <- exp(m + .5 * diag(Krest))
    list(mean = b, covariance = tcrossprod(b) * expm1(Krest),
      latent_mean = m, latent_covariance = Krest)
  })
  mixture <- baseline_combine(pieces, rule$weights / sum(rule$weights))
  stopifnot(max(abs(mixture$mean - 1)) < 1e-10,
    max(abs(mixture$covariance - expm1(case$covariance))) < 1e-10,
    max(abs(mixture$latent_mean - case$prior$offset)) < 1e-12,
    max(abs(mixture$latent_covariance - case$covariance)) < 1e-12)
}
writeLines("PASS: same Gaussian prior, mean-one birth prior and complete latent/birth covariance at n=4/8",
  file.path(output, "prior_validation.txt"))
limit <- as.integer(Sys.getenv("ALFAK2_BASELINE_CASE_LIMIT", "16"))
statuses <- list()
for (i in seq_len(min(limit, length(cases)))) {
  case <- cases[[i]]
  ref_file <- file.path(reference, paste0(case$id, "_summary.rds"))
  if (!file.exists(ref_file)) stop("Reference has not returned for ", case$id)
  r <- readRDS(ref_file)
  if (!isTRUE(r$reference$convergence$passed)) stop("Uncertified MCMC reference for ", case$id)
  result_file <- file.path(output, paste0(case$id, "_integration.rds"))
  case_key <- list(source = key, case = round2_hash_object(case), reference = unname(round2_sha256(ref_file)))
  started <- proc.time()[["elapsed"]]
  if (file.exists(result_file)) {
    result <- readRDS(result_file)
    if (!identical(result$cache_key, case_key)) stop("Candidate result fingerprint mismatch")
  } else {
    result <- tryCatch(baseline_integrate(case, r$laplace, runtime,
      cache_directory = cache, cache_key = case_key), error = function(e)
        list(status = "integration_error", message = conditionMessage(e)))
    result$cache_key <- case_key
    result$elapsed_seconds <- proc.time()[["elapsed"]] - started
    saveRDS(result, result_file)
  }
  statuses[[case$id]] <- data.frame(case_id = case$id, status = result$status,
    elapsed_seconds = result$elapsed_seconds,
    final_order = if (!is.null(result$history) && nrow(result$history)) tail(result$history$order, 1L) else NA_integer_,
    conditional_nodes = if (!is.null(result$conditional_status)) nrow(result$conditional_status) else 0L,
    message = if (is.null(result$message)) "" else result$message)
  write.csv(do.call(rbind, statuses), file.path(output, "integration_status.csv"), row.names = FALSE)
  if (!is.null(result$history)) write.csv(result$history,
    file.path(output, paste0(case$id, "_refinement.csv")), row.names = FALSE)
  if (!is.null(result$conditional_status)) write.csv(result$conditional_status,
    file.path(output, paste0(case$id, "_conditional_status.csv")), row.names = FALSE)
  cat(case$id, result$status, "order", statuses[[case$id]]$final_order,
    "time", result$elapsed_seconds, "s\n")
}
