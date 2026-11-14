#!/usr/bin/env Rscript
# Stage M, step 0: where the Laplace approximation fails, and whether
# integrating the initial abundance of non-major states out would remove it.
#
#   Rscript validation/sparse_joint/model_m/m0.R --stage=fits --set=degenerate|healthy|cliff [--workers=5]
#   Rscript validation/sparse_joint/model_m/m0.R --stage=decompose|cliffs|quadrature [--workers=5]
#   Rscript validation/sparse_joint/model_m/m0.R --stage=report
#
# Sets (development benchmark only):
#   degenerate  the 91 candidate histories C1 classed degenerate (ff13cf4);
#   healthy     40 that C1 classed optimum, 4 per scenario (ff13cf4), for
#               contrast;
#   cliff       the candidates cv1 did not accept after its restarts, in the
#               dev items it did not accept (cv1).
# Results in results/m0/ (git-ignored).

mm_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
mm_repo <- normalizePath(file.path(mm_script_dir, "..", "..", ".."))
mm_out <- file.path(mm_script_dir, "results", "m0")
mm_lib <- file.path(mm_repo, "validation", "sparse_joint", "multi_timepoint", "results", "lib")
mm_datasets <- file.path(mm_repo, "benchmark", "results", "sparse_two_timepoint", "datasets")
mm_c1 <- file.path(mm_repo, "validation", "sparse_joint", "convergence", "results", "ff13cf4", "diagnose",
                   "candidates.csv")
mm_c3 <- file.path(mm_repo, "validation", "sparse_joint", "multi_timepoint", "results", "reference_hpc", "cv1", "dev")
source(file.path(mm_script_dir, "reference.R"))

mm_args <- function() {
  out <- list()
  for (arg in commandArgs(trailingOnly = TRUE)) {
    body <- sub("^--", "", arg)
    out[[gsub("-", "_", sub("=.*$", "", body))]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
  }
  out
}

mm_load <- function(build) {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  .libPaths(c(file.path(mm_lib, paste0("alfak2_", build)), .libPaths()))
  loadNamespace("alfak2")
  cat("alfak2 from", dirname(find.package("alfak2")), "\n")
}

mm_parallel <- function(x, f, workers) {
  parallel::mclapply(x, function(item) tryCatch(f(item), error = function(e) {
    list(error = conditionMessage(e), item = item)
  }), mc.cores = workers, mc.preschedule = FALSE)
}

mm_save <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(object, temporary)
  file.rename(temporary, path)
}

# --- fits -----------------------------------------------------------------

mm_set <- function(set) {
  if (set %in% c("degenerate", "healthy")) {
    c1 <- utils::read.csv(mm_c1, stringsAsFactors = FALSE)
    if (set == "degenerate") return(c1[c1$class == "degenerate", c("id", "history", "f_nlminb", "eig_warm")])
    pool <- c1[c1$class == "optimum" & c1$converged, ]
    set.seed(20261114)
    keep <- unlist(lapply(split(seq_len(nrow(pool)), pool$scenario), function(i) i[sample.int(length(i), min(4L, length(i)))]))
    return(pool[sort(keep), c("id", "history", "f_nlminb", "eig_warm")])
  }
  if (set == "cliff") {
    rows <- list()
    for (file in list.files(mm_c3, pattern = "[.]rds$")) {
      record <- readRDS(file.path(mm_c3, file))
      if (!identical(as.numeric(record$diagnostics$convergence), 1)) next
      profile <- record$diagnostics$history_profile
      histories <- if (is.null(profile)) record$diagnostics$history else profile$history[profile$convergence == 1]
      rows[[length(rows) + 1L]] <- data.frame(id = sub("[.]rds$", "", file), history = histories, f_nlminb = NA_real_,
                                              eig_warm = NA_real_, stringsAsFactors = FALSE)
    }
    return(do.call(rbind, rows))
  }
  stop("unknown set ", set)
}

mm_stage_fits <- function(set, workers) {
  build <- if (set == "cliff") "cv1" else "ff13cf4"
  mm_load(build)
  todo <- mm_set(set)
  out <- file.path(mm_out, "fits", set)
  todo$path <- file.path(out, sprintf("%s__%g.rds", todo$id, todo$history))
  todo <- todo[!file.exists(todo$path), ]
  cat("Set", set, "build", build, ":", nrow(todo), "candidate(s) to fit\n")
  status <- mm_parallel(split(todo, seq_len(nrow(todo))), function(row) {
    dataset <- readRDS(file.path(mm_datasets, paste0(row$id, ".rds")))
    counts <- dataset$counts
    storage.mode(counts) <- "integer"
    history <- row$history
    candidates <- mm_capture(quote(alfak2::fit_sparse_joint(counts, dt = dataset$interval, beta = dataset$beta,
                                                             birth_rate = 1, min_cn = 1L, max_cn = 8L,
                                                             history = history)))
    candidate <- candidates[[length(candidates)]]
    candidate$id <- row$id
    candidate$set <- set
    candidate$build <- build
    candidate$c1_objective <- row$f_nlminb
    mm_save(candidate, row$path)
    "ok"
  }, workers)
  failed <- Filter(function(x) !identical(x, "ok"), status)
  cat("Done.", length(failed), "failed.\n")
  for (f in failed) cat(" ", f$item$id, f$item$history, f$error, "\n")
}

# --- decompose ------------------------------------------------------------

mm_blocks <- function(names, n_major) {
  index <- ave(seq_along(names), names, FUN = seq_along)
  ifelse(names %in% c("zeta", "g"), paste0(names, ifelse(index <= n_major, "_major", "_minor")), names)
}

# The Hessian of the joint density at the mode, split by term (see the
# handoff entry of 2026-11-14 02:06): prior with the balance prior at zero
# residual, the balance prior's second-order part, the likelihood with
# counts at their expectation, and the residual part of the likelihood.
mm_hessian_parts <- function(candidate, objective, p) {
  d <- candidate$tmb_data
  n <- length(d$y0)
  mode <- candidate$mode
  H <- function(data) as.matrix(mm_make_objective(candidate, data)$env$spHess(mode, random = TRUE))
  report <- objective$report(mode)
  prior_only <- d
  prior_only$weight <- rep(0, n)
  prior_gn <- prior_only
  prior_gn$balance_offset <- report$balance_residual
  expected <- d
  expected$balance_offset <- report$balance_residual
  expected$y0 <- exp(p$zeta)
  expected$y1 <- exp(report$log_intensity1)
  H1 <- as.matrix(objective$env$spHess(mode, random = TRUE))
  H2 <- H(prior_only)
  H3 <- H(prior_gn)
  H4 <- H(expected)
  list(total = H1, prior_gn = H3, balance_second = H2 - H3, likelihood_expected = H4 - H3,
       likelihood_residual = (H1 - H2) - (H4 - H3), report = report)
}

mm_decompose_one <- function(path) {
  candidate <- readRDS(path)
  d <- candidate$tmb_data
  objective <- mm_make_objective(candidate)
  p <- mm_unpack(objective, candidate$mode)
  joint_tmb <- as.numeric(objective$env$f(candidate$mode, order = 0))
  joint_r <- mm_joint_nll(p, d)
  parts <- mm_hessian_parts(candidate, objective, p)
  e <- eigen(parts$total, symmetric = TRUE)
  k <- length(e$values)
  v <- e$vectors[, k]
  along <- vapply(parts[c("total", "prior_gn", "balance_second", "likelihood_expected", "likelihood_residual")],
                  function(M) as.numeric(t(v) %*% M %*% v), numeric(1))
  names_random <- names(candidate$mode)[objective$env$random]
  block <- mm_blocks(names_random, d$n_major)
  loading <- tapply(v^2, block, sum)
  # Leading non-major states of the flat direction.
  zi <- which(names_random == "zeta")
  gi <- which(names_random == "g")
  n <- length(d$y0)
  minor <- seq_len(n) > d$n_major
  top <- order(-abs(v[zi]) * minor)[1:3]
  report <- parts$report
  net <- report$net
  from <- d$edge_from + 1L
  to <- d$edge_to + 1L
  log_phi <- mm_log_phi1(net[to], net[from], d$dt)
  inflow <- log(2) + d$edge_log_prob + p$g[from] + p$zeta[from] + log_phi
  own <- p$zeta + net * d$dt
  states <- do.call(rbind, lapply(top, function(j) {
    src <- which(to == j)
    total <- mm_logsumexp(c(own[j], inflow[src]))
    share_own <- exp(own[j] - total)
    share_major <- sum(exp(inflow[src][from[src] <= d$n_major] - total))
    lambda1 <- exp(report$log_intensity1[j])
    data.frame(state = j, loading = v[zi][j], y0 = d$y0[j], y1 = d$y1[j], expected0 = exp(p$zeta[j]),
               expected1 = lambda1, share_own = share_own, share_major = share_major,
               share_minor = 1 - share_own - share_major,
               diag_total = parts$total[zi[j], zi[j]], diag_prior = parts$prior_gn[zi[j], zi[j]],
               diag_residual = parts$likelihood_residual[zi[j], zi[j]],
               residual_formula = -(d$y1[j] - lambda1) * share_own * (1 - share_own),
               balance_residual = report$balance_residual[j])
  }))
  list(id = candidate$id, history = candidate$history, set = candidate$set, n_states = n, n_major = d$n_major,
       n_random = length(objective$env$random), objective = candidate$objective, c1_objective = candidate$c1_objective,
       joint_tmb = joint_tmb, joint_r = joint_r, eigen = rev(e$values)[1:3], along = along, loading = loading,
       states = states, vector_g = sum(v[gi]^2), vector_top = data.frame(variable = block[order(-abs(v))[1:8]],
                                                                         loading = v[order(-abs(v))[1:8]]))
}

mm_stage_decompose <- function(workers) {
  mm_load("ff13cf4")
  paths <- unlist(lapply(c("degenerate", "healthy"), function(set) {
    list.files(file.path(mm_out, "fits", set), pattern = "[.]rds$", full.names = TRUE)
  }))
  out <- file.path(mm_out, "decompose")
  todo <- paths[!file.exists(file.path(out, basename(dirname(paths)), basename(paths)))]
  cat("Decompose:", length(todo), "of", length(paths), "\n")
  status <- mm_parallel(todo, function(path) {
    result <- mm_decompose_one(path)
    mm_save(result, file.path(out, basename(dirname(path)), basename(path)))
    "ok"
  }, workers)
  failed <- Filter(function(x) !identical(x, "ok"), status)
  cat("Done.", length(failed), "failed.\n")
  for (f in failed) cat(" ", f$item, f$error, "\n")
}

# --- cliffs ---------------------------------------------------------------

mm_cliff_one <- function(path) {
  candidate <- readRDS(path)
  d <- candidate$tmb_data
  objective <- mm_make_objective(candidate)
  env <- objective$env
  random <- env$random
  names_random <- names(candidate$mode)[random]
  block <- mm_blocks(names_random, d$n_major)
  zi <- which(names_random == "zeta")
  gi <- which(names_random == "g")
  assign(".mm_start", candidate$mode[random], envir = env)
  assign("random.start", expression(.mm_start), envir = env)
  at <- function(par) {
    f <- tryCatch(as.numeric(objective$fn(par)), error = function(e) NA_real_)
    mode <- env$last.par
    list(f = f, mode = mode, joint = if (is.finite(f)) as.numeric(env$f(mode, order = 0)) else NA_real_)
  }
  base <- at(candidate$par)
  rule <- list(type = "grid", step = 0.2, width = 8)
  marginal_at <- function(mode) mm_marginal_nll(mm_unpack(objective, mode), d, rule)
  base_marginal <- marginal_at(base$mode)
  h <- 1e-3
  rows <- list()
  for (k in seq_along(candidate$par)) {
    for (sign in c(-1, 1)) {
      side <- at(candidate$par + sign * h * (seq_along(candidate$par) == k))
      delta <- (side$mode - base$mode)[random]
      moved_zeta <- order(-abs(delta[zi]))[1:3]
      rows[[length(rows) + 1L]] <- data.frame(
        id = candidate$id, history = candidate$history, scale = k, sign = sign,
        f_change = side$f - base$f, joint_change = side$joint - base$joint,
        marginal_change = if (is.finite(side$f)) marginal_at(side$mode) - base_marginal else NA_real_,
        max_zeta_minor = max(abs(delta[block == "zeta_minor"])), max_zeta_major = max(abs(delta[block == "zeta_major"])),
        max_g = max(abs(delta[gi])), zeta_moved_05 = sum(abs(delta[zi]) > 0.5),
        top_major = moved_zeta[[1L]] <= d$n_major, top_y0 = d$y0[moved_zeta[[1L]]], top_y1 = d$y1[moved_zeta[[1L]]],
        top_change = delta[zi][moved_zeta[[1L]]], stringsAsFactors = FALSE)
    }
  }
  list(id = candidate$id, history = candidate$history, check = candidate$check, f = base$f,
       steps = do.call(rbind, rows))
}

mm_stage_cliffs <- function(workers) {
  mm_load("cv1")
  paths <- list.files(file.path(mm_out, "fits", "cliff"), pattern = "[.]rds$", full.names = TRUE)
  out <- file.path(mm_out, "cliffs")
  todo <- paths[!file.exists(file.path(out, basename(paths)))]
  cat("Cliffs:", length(todo), "of", length(paths), "\n")
  status <- mm_parallel(todo, function(path) {
    mm_save(mm_cliff_one(path), file.path(out, basename(path)))
    "ok"
  }, workers)
  failed <- Filter(function(x) !identical(x, "ok"), status)
  cat("Done.", length(failed), "failed.\n")
  for (f in failed) cat(" ", f$item, f$error, "\n")
}

# --- quadrature -----------------------------------------------------------

mm_rules <- list(
  gh20 = list(type = "gh", nodes = 20L), gh40 = list(type = "gh", nodes = 40L), gh80 = list(type = "gh", nodes = 80L),
  grid_0.3_7 = list(type = "grid", step = 0.3, width = 7), grid_0.2_8 = list(type = "grid", step = 0.2, width = 8),
  grid_0.1_10 = list(type = "grid", step = 0.1, width = 10), mixture20 = list(type = "mixture", nodes = 20L)
)

mm_quadrature_errors <- function(p, d, label) {
  exact <- mm_marginal_nll(p, d, list(type = "exact"), parts = TRUE)
  minor <- exact$inputs$minor
  out <- data.frame(label = label, state = minor, y0 = d$y0[minor], y1 = d$y1[minor], weight = d$weight[minor],
                    m = exact$inputs$m[minor], s = exact$inputs$s[minor], data_centre = exact$inputs$centre[minor],
                    exact = exact$terms, stringsAsFactors = FALSE)
  for (name in names(mm_rules)) out[[name]] <- mm_marginal_nll(p, d, mm_rules[[name]], parts = TRUE)$terms - exact$terms
  out
}

mm_quadrature_one <- function(path) {
  candidate <- readRDS(path)
  d <- candidate$tmb_data
  objective <- mm_make_objective(candidate)
  p <- mm_unpack(objective, candidate$mode)
  set.seed(7)
  # At the mode, and moved as an optimiser might move it.
  variants <- list(mode = p)
  for (r in 1:2) {
    q <- p
    q$theta1 <- q$theta1 + c(-1, 1)[r]
    # Net growth over the interval and the history moves by about 0.3.
    q$g <- q$g + stats::rnorm(length(q$g), 0, 0.3 / max(d$dt, d$history))
    major <- seq_len(length(q$zeta)) <= d$n_major
    q$zeta[major] <- q$zeta[major] + stats::rnorm(sum(major), 0, 0.5)
    variants[[paste0("moved", r)]] <- q
  }
  rows <- do.call(rbind, lapply(names(variants), function(v) mm_quadrature_errors(variants[[v]], d, v)))
  cbind(id = candidate$id, history = candidate$history, set = candidate$set, rows, stringsAsFactors = FALSE)
}

# Extreme states: many cells where the balance value is far below, and the
# reverse, under both priors (balance sd 1, free sd 3).
mm_quadrature_stress <- function() {
  grid <- expand.grid(y0 = c(0, 1, 5, 30, 149), y1 = c(0, 2, 10, 100, 300), m = c(-10, -4, -1, 3), s = c(1, 3),
                      log_inflow = c(-Inf, -3, 1), a = c(-0.3, 0.3))
  K <- nrow(grid)
  inputs <- list(minor = seq_len(K), a = grid$a, log_inflow = grid$log_inflow, c1 = 0)
  d <- list(y0 = grid$y0, y1 = grid$y1, weight = rep(1, K))
  l <- mm_state_loglik(inputs, d)
  prec0 <- d$y0 + 0.5
  own <- pmax(d$y1 - exp(grid$log_inflow), 0)
  prec1 <- own + 0.5
  centre <- (prec0 * log(prec0) + prec1 * (log(prec1) - grid$a)) / (prec0 + prec1)
  spread <- 1 / sqrt(prec0 + prec1)
  exact <- mm_quadrature(l, grid$m, grid$s, list(type = "exact"))
  for (name in names(mm_rules)) {
    grid[[name]] <- mm_quadrature(l, grid$m, grid$s, mm_rules[[name]], centre, spread) - exact
  }
  # The truth itself, against stats::integrate, on a sample.
  set.seed(3)
  check <- sample.int(K, 60)
  grid$integrate_check <- NA_real_
  grid$integrate_check[check] <- vapply(check, function(k) {
    f <- function(z) exp(stats::dnorm(z, grid$m[k], grid$s[k], log = TRUE) + l(matrix(z, 1), k) - exact[k])
    log(stats::integrate(f, -Inf, Inf, rel.tol = 1e-10, subdivisions = 2000L)$value)
  }, numeric(1))
  grid$exact <- exact
  grid
}

mm_stage_quadrature <- function(workers) {
  mm_load("ff13cf4")
  paths <- unlist(lapply(c("degenerate", "healthy", "cliff"), function(set) {
    list.files(file.path(mm_out, "fits", set), pattern = "[.]rds$", full.names = TRUE)
  }))
  out <- file.path(mm_out, "quadrature")
  todo <- paths[!file.exists(file.path(out, basename(dirname(paths)), basename(paths)))]
  cat("Quadrature:", length(todo), "of", length(paths), "\n")
  status <- mm_parallel(todo, function(path) {
    mm_save(mm_quadrature_one(path), file.path(out, basename(dirname(path)), basename(path)))
    "ok"
  }, workers)
  failed <- Filter(function(x) !identical(x, "ok"), status)
  cat("Done.", length(failed), "failed.\n")
  for (f in failed) cat(" ", f$item, f$error, "\n")
  mm_save(mm_quadrature_stress(), file.path(out, "stress.rds"))
}

mm_main <- function() {
  args <- mm_args()
  workers <- if (is.null(args$workers)) 5L else as.integer(args$workers)
  stage <- args$stage
  if (identical(stage, "fits")) return(mm_stage_fits(args$set, workers))
  if (identical(stage, "decompose")) return(mm_stage_decompose(workers))
  if (identical(stage, "cliffs")) return(mm_stage_cliffs(workers))
  if (identical(stage, "quadrature")) return(mm_stage_quadrature(workers))
  if (identical(stage, "report")) {
    source(file.path(mm_script_dir, "m0_report.R"))
    return(mm_report())
  }
  stop("--stage must be fits, decompose, cliffs, quadrature or report")
}

if (sys.nframe() == 0L) mm_main()
