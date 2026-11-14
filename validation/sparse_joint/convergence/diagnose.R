#!/usr/bin/env Rscript
# Stage C, step 1 (handoff HANDOFF_CLAUDE_2026-11-13_multi_timepoint.md): why
# the outer nlminb of fit_sparse_joint() stops without convergence.
#
#   Rscript validation/sparse_joint/convergence/diagnose.R [--build=ff13cf4] [--workers=6] [--ids=REGEX]
#
# Refits the 120 development-benchmark items of the stage-0 reference set
# (patient 1 of the 12 development landscapes, 8 core and 2 short-interval
# scenarios) with the package build `results/lib/alfak2_<build>` of
# ../multi_timepoint, at the defaults, and probes every candidate history
# where its scale optimisation ended:
#
# - nlminb's own result (code, message, iterations, evaluations);
# - "warm": the objective at the end point with TMB's inner start
#   (last.par.best, as in the fit), the smallest eigenvalues of the inner
#   Hessian there and the latent variables along the flattest direction;
# - "fresh": the objective at the same point with the inner optimisation
#   started from the initial latent values, the smallest eigenvalue, the TMB
#   gradient and a central finite-difference gradient (step 1e-3), and the
#   largest decrease found by coordinate steps of +-0.05;
# - "start": the smallest eigenvalue at the prior-centred start.
#
# The probe restores TMB's state afterwards, so the fit goes on as it would;
# its objective is checked against the stage-0 record. One file per item in
# results/<build>/diagnose/. Nothing here changes the package.

cv_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
cv_repo <- normalizePath(file.path(cv_script_dir, "..", "..", ".."))
cv_mt <- file.path(cv_repo, "validation", "sparse_joint", "multi_timepoint", "results")
cv_benchmark <- file.path(cv_repo, "benchmark", "results", "sparse_two_timepoint")

cv_args <- function() {
  args <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    key <- gsub("-", "_", sub("=.*$", "", body))
    args[[key]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  args
}

# Smallest eigenvalues of the inner Hessian at a mode, and optionally the
# variables that load most on the flattest direction.
cv_curvature <- function(env, mode, vectors = FALSE) {
  H <- as.matrix(env$spHess(mode, random = TRUE))
  e <- eigen(H, symmetric = TRUE, only.values = !vectors)
  k <- length(e$values)
  # Joint negative log density at the mode (without the log determinant) and
  # how well the inner problem converged there.
  joint <- tryCatch(as.numeric(env$f(mode, order = 0)), error = function(e) NA_real_)
  inner_grad <- tryCatch(max(abs(env$f(mode, order = 1)[env$random])), error = function(e) NA_real_)
  out <- list(min_eig = e$values[[k]], second_eig = if (k > 1L) e$values[[k - 1L]] else NA_real_,
              max_eig = e$values[[1L]], joint = joint, inner_grad = inner_grad,
              half_logdet = sum(log(pmax(e$values, .Machine$double.xmin))) / 2)
  if (vectors) {
    rnd <- env$random
    nm <- names(env$par)[rnd]
    within <- stats::ave(seq_along(nm), nm, FUN = seq_along)
    v <- e$vectors[, k]
    top <- order(-abs(v))[1:3]
    n_major <- env$data$n_major
    out$direction <- data.frame(variable = nm[top], index = within[top], loading = v[top],
                                major = nm[top] %in% c("zeta", "g") & within[top] <= n_major,
                                value = unname(mode[rnd][top]),
                                y0 = ifelse(nm[top] %in% c("zeta", "g"), env$data$y0[pmin(within[top], length(env$data$y0))], NA),
                                y1 = ifelse(nm[top] %in% c("zeta", "g"), env$data$y1[pmin(within[top], length(env$data$y1))], NA),
                                stringsAsFactors = FALSE)
  }
  out
}

cv_probe <- function(objective, optimum, lower, upper) {
  env <- objective$env
  saved <- mget(c("last.par.best", "value.best", "last.par", "random.start"), envir = env)
  on.exit(for (k in names(saved)) assign(k, saved[[k]], envir = env), add = TRUE)
  started <- proc.time()[["elapsed"]]
  par <- optimum$par
  rnd <- env$random
  record <- list(history = env$data$history, n_random = length(rnd), n_states = length(env$data$y0),
                 n_major = env$data$n_major, par = par, lower = lower, upper = upper,
                 nlminb = optimum[intersect(names(optimum), c("objective", "convergence", "message", "iterations",
                                                              "evaluations"))])
  record$at_lower <- par <= lower + 1e-6
  record$at_upper <- par >= upper - 1e-6
  # Warm: as the fit sees it.
  f_warm <- tryCatch(as.numeric(objective$fn(par)), error = function(e) NA_real_)
  record$warm <- c(list(f = f_warm), tryCatch(cv_curvature(env, env$last.par, vectors = TRUE),
                                              error = function(e) list(error = conditionMessage(e))))
  record$warm$tmb_grad <- tryCatch(as.numeric(objective$gr(par)), error = function(e) rep(NA_real_, length(par)))
  # Fresh: the inner optimisation always starts from the initial values.
  assign(".cv_fixed_start", env$par[rnd], envir = env)
  assign("random.start", expression(.cv_fixed_start), envir = env)
  fresh_fn <- function(p) tryCatch(as.numeric(objective$fn(p)), error = function(e) NA_real_)
  f_fresh <- fresh_fn(par)
  record$fresh <- c(list(f = f_fresh), tryCatch(cv_curvature(env, env$last.par),
                                                error = function(e) list(error = conditionMessage(e))))
  record$fresh$tmb_grad <- tryCatch(as.numeric(objective$gr(par)), error = function(e) rep(NA_real_, length(par)))
  h <- 1e-3
  fd <- vapply(seq_along(par), function(i) {
    e <- replace(numeric(length(par)), i, h)
    (fresh_fn(par + e) - fresh_fn(par - e)) / (2 * h)
  }, numeric(1))
  projected <- fd
  projected[record$at_lower & fd > 0] <- 0
  projected[record$at_upper & fd < 0] <- 0
  record$fresh$fd_grad <- fd
  record$fresh$projected_grad <- projected
  steps <- unlist(lapply(seq_along(par), function(i) lapply(c(-0.05, 0.05), function(s) {
    p <- par
    p[[i]] <- min(max(p[[i]] + s, lower), upper)
    p
  })), recursive = FALSE)
  record$fresh$probe_decrease <- max(f_fresh - vapply(steps, fresh_fn, numeric(1)), na.rm = TRUE)
  # Start.
  f_start <- fresh_fn(objective$par)
  record$start <- c(list(f = f_start), tryCatch(cv_curvature(env, env$last.par),
                                                error = function(e) list(error = conditionMessage(e))))
  record$probe_seconds <- proc.time()[["elapsed"]] - started
  .cv_runs[[length(.cv_runs) + 1L]] <<- record
  invisible(NULL)
}

cv_item <- function(id) {
  dataset <- readRDS(file.path(cv_benchmark, "datasets", paste0(id, ".rds")))
  counts <- dataset$counts
  storage.mode(counts) <- "integer"
  .cv_runs <<- list()
  started <- proc.time()[["elapsed"]]
  fit <- alfak2::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta, birth_rate = 1,
                                  min_cn = 1L, max_cn = 8L)
  d <- fit$diagnostics
  list(id = id, seconds = proc.time()[["elapsed"]] - started, candidates = .cv_runs,
       selected = list(history = d$history, objective = d$objective, convergence = d$convergence,
                       history_profile = d$history_profile))
}

cv_main <- function() {
  args <- cv_args()
  build <- if (is.null(args$build)) "ff13cf4" else args$build
  workers <- if (is.null(args$workers)) 6L else as.integer(args$workers)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  .libPaths(c(file.path(cv_mt, "lib", paste0("alfak2_", build)), .libPaths()))
  loadNamespace("alfak2")
  cat("alfak2 from", dirname(find.package("alfak2")), "\n")
  trace("sparse_joint_optimise_scales", where = asNamespace("alfak2"), print = FALSE,
        exit = quote(cv_probe(objective, optimum, lower, upper)))
  reference <- file.path(cv_mt, "reference", "ff13cf4", "dev")
  ids <- sub("[.]rds$", "", list.files(reference, pattern = "[.]rds$"))
  if (!is.null(args$ids)) ids <- grep(args$ids, ids, value = TRUE)
  # Largest first, so that the workers finish together.
  size <- vapply(ids, function(id) readRDS(file.path(reference, paste0(id, ".rds")))$diagnostics$n_states, numeric(1))
  ids <- ids[order(-size)]
  output <- if (is.null(args$output_dir)) file.path(cv_script_dir, "results", build, "diagnose") else args$output_dir
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  ids <- ids[!file.exists(file.path(output, paste0(ids, ".rds")))]
  cat("Items:", length(ids), "\n")
  status <- parallel::mclapply(ids, function(id) {
    result <- tryCatch(cv_item(id), error = function(e) list(id = id, error = conditionMessage(e)))
    ref <- readRDS(file.path(reference, paste0(id, ".rds")))
    result$reference_objective <- ref$diagnostics$objective
    path <- file.path(output, paste0(id, ".rds"))
    temporary <- paste0(path, ".tmp", Sys.getpid())
    saveRDS(result, temporary)
    file.rename(temporary, path)
    if (is.null(result$error)) "ok" else paste(id, result$error)
  }, mc.cores = workers, mc.preschedule = FALSE)
  failed <- Filter(function(x) !identical(x, "ok"), status)
  cat("Done.", length(failed), "item(s) failed.\n")
  if (length(failed)) cat(unlist(failed), sep = "\n")
}

if (sys.nframe() == 0L) cv_main()
