#!/usr/bin/env Rscript
# Fixed-data approximation diagnostic; no reserved accuracy data are consumed.
repo <- normalizePath(".")
directory <- file.path(repo, "validation/posterior_reference_4to8")
source(file.path(directory, "reference.R"))
source(file.path(repo, "validation/accuracy_round2/common.R"))
runtime <- reference_load(repo)
output <- Sys.getenv("ALFAK2_REFERENCE_OUTPUT", file.path(directory, "results"))
raw <- Sys.getenv("ALFAK2_REFERENCE_RAW", "/private/tmp/alfak2_reference_4to8_20261104")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(raw, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(file.path(output, "validation/validation_tests.csv"))) stop("Run validate.R first")
validation <- read.csv(file.path(output, "validation/validation_tests.csv"))
if (any(validation$failed > 0 | validation$error | validation$warning > 0)) stop("Reference validation has not passed")
cases <- list()
for (n in c(4L, 8L)) for (profile in c("balanced", "sparse")) for (beta in c(0, .005)) for (growth in c(FALSE, TRUE)) {
  case <- reference_case(n, profile, beta, growth, runtime)
  cases[[case$id]] <- case
}
case_file <- file.path(output, "frozen_cases.rds")
if (file.exists(case_file) && !identical(readRDS(case_file), cases)) stop("Frozen reference cases changed")
if (!file.exists(case_file)) saveRDS(cases, case_file, version = 2)
code <- c(file.path(directory, c("PROTOCOL.md", "reference.R", "reference.cpp", "run.R", "validate.R")),
  list.files(file.path(repo, "R"), "[.]R$", full.names = TRUE), runtime$native, case_file)
hash <- round2_sha256(code)
freeze <- data.frame(file = names(hash), sha256 = unname(hash), stringsAsFactors = FALSE)
freeze_file <- file.path(output, "source_sha256.csv")
if (file.exists(freeze_file) && !identical(read.csv(freeze_file), freeze)) stop("Source/case fingerprint changed; preserve the earlier run")
if (!file.exists(freeze_file)) write.csv(freeze, freeze_file, row.names = FALSE)
computation_key <- round2_hash_object(freeze)
saveRDS(sessionInfo(), file.path(output, "session_info.rds"))
status <- list()
for (i in seq_along(cases)) {
  case <- cases[[i]]
  target_file <- file.path(raw, paste0(case$id, "_complete.rds"))
  cache_key <- list(computation_sha256 = computation_key, case_sha256 = round2_hash_object(case))
  if (file.exists(target_file)) {
    result <- readRDS(target_file)
    if (!identical(result$cache_key, cache_key)) stop("Raw-chain cache fingerprint mismatch; use a new ALFAK2_REFERENCE_RAW directory")
  } else {
    started <- proc.time()[["elapsed"]]
    fit <- reference_fit(case, runtime)
    if (!identical(fit$diagnostics$covariance_status, "TMB_sdreport")) stop("Untrusted Laplace comparator: ", case$id)
    fit$optimizer <- NULL
    chains <- NULL; history <- list(); retained <- 0L
    for (length in c(20000L, 40000L, 80000L, 160000L)) {
      chains <- reference_chains(case, runtime, retained = length - retained,
        seed = 720000L + 1000L * i + length, previous = chains)
      retained <- length
      d <- reference_diagnostics(chains, case)
      history[[as.character(length)]] <- d$convergence
      cat(case$id, "retained", length, "Rhat", d$convergence$max_rhat,
        "min bulk/tail ESS", d$convergence$min_bulk_ess, d$convergence$min_tail_ess,
        "passed", d$convergence$passed, "\n")
      if (isTRUE(d$convergence$passed)) break
    }
    result <- list(cache_key = cache_key, case_id = case$id, case_index = i, laplace = fit,
      reference = d, chains = chains, convergence_history = do.call(rbind, history),
      elapsed_seconds = proc.time()[["elapsed"]] - started)
    saveRDS(result, target_file)
  }
  compact <- result; compact$chains <- NULL
  compact$raw_chain_file <- target_file
  compact$raw_chain_sha256 <- unname(round2_sha256(target_file))
  saveRDS(compact, file.path(output, paste0(case$id, "_summary.rds")))
  write.csv(result$reference$summary, file.path(output, paste0(case$id, "_diagnostics.csv")), row.names = FALSE)
  write.csv(result$reference$product_summary, file.path(output, paste0(case$id, "_covariance_diagnostics.csv")), row.names = FALSE)
  write.csv(result$convergence_history, file.path(output, paste0(case$id, "_convergence_history.csv")), row.names = FALSE)
  status[[case$id]] <- cbind(result$reference$convergence,
    data.frame(n = case$n, profile = case$profile, beta = case$beta, growth = case$growth,
      laplace_covariance = result$laplace$diagnostics$covariance_status,
      elapsed_seconds = result$elapsed_seconds, raw_chain_file = target_file,
      raw_chain_sha256 = compact$raw_chain_sha256))
  write.csv(do.call(rbind, status), file.path(output, "case_status.csv"), row.names = FALSE)
  cat(sprintf("Completed case %d/%d: %s (%.2fs)\n", i, length(cases), case$id, result$elapsed_seconds))
}
cat("All fixed-data cases returned; inspect convergence before certifying any reference.\n")
