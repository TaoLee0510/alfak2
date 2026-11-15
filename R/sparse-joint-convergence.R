# Convergence of the scale optimisation of the sparse joint models.
#
# The scales (dosage, curvature, interaction, residual) are fitted by Laplace
# marginal likelihood: an outer nlminb over the log scales, an inner Newton
# optimisation of the latent variables. The inner problem is not convex. A
# karyotype with no cell in the first sample and one in the second can owe it
# to having been there already or to mutation from a neighbour, and between
# the two explanations the joint density is nearly flat. Where the inner
# Hessian is nearly singular the log determinant of the Laplace approximation
# runs off towards minus infinity: the objective drops for no reason in the
# data, and TMB's gradient there is out by many orders of magnitude. The inner
# optimisation also starts from the best mode so far (TMB's `last.par.best`),
# so the objective nlminb sees depends on the path it took. nlminb then stops
# with "false convergence" at a point that is no optimum, or at one that is.
#
# nlminb's verdict is therefore not taken as it is. A point is accepted when
#   - the objective there, from the fit's own inner mode, reproduces the one
#     nlminb reported (to `objective`);
#   - the inner Hessian at that mode is healthy (smallest eigenvalue at least
#     `eigen`);
#   - an inner optimisation started from the initial latent values finds no
#     mode with a better joint density (by more than `objective`);
#   - the central finite-difference gradient (step `step`) of the objective
#     that follows the fit's mode (inner optimisations started from it),
#     projected on the bounds, is at most `gradient`. TMB's own gradient is
#     not used: the inner Hessian is badly conditioned (condition numbers to
#     1e8), the inner optimisation stops at gradients near 1e-7, and TMB's
#     gradient, which takes the inner mode as exact, is then out by tens even
#     where the objective is flat.
# The fit is exactly what it was without the check (alfak2 ff13cf4): nlminb's
# point, or the prior-centred scales if nlminb fails; the check only labels
# it, and a point that is not accepted is reported as not converged, with the
# reasons and the flat direction of the inner Hessian. Searching again is
# internal and off by default (option `alfak2.convergence_restarts`): the
# scales are searched again with the inner start fixed, so that the objective
# is a function of the scales alone, and points with a near-singular inner
# Hessian rejected, by a bounded compass search from the end of the first run
# and, if that is not accepted, from the prior-centred start; the best
# accepted point is kept, or, when none is, the best point with a healthy
# inner Hessian.
#
# Tolerances from the 420 candidate fits of the development benchmark of
# benchmark/sparse_two_timepoint at alfak2 ff13cf4
# (validation/sparse_joint/convergence/diagnose.R): the smallest inner
# eigenvalue was at least 4.4e-3 where nlminb converged and at most 6.9e-4
# where the Laplace approximation had degenerated; the projected
# finite-difference gradient was at most 0.029 at optima where nlminb
# converged and 0.044 where it stopped at an optimum without saying so.
sparse_joint_convergence_tolerance <- c(eigen = 2e-3, objective = 1e-3, gradient = 0.05, step = 1e-3)

# TMB's state between evaluations (TMB 1.9: MakeADFun's environment).
sparse_joint_tmb_state_names <- c("last.par", "last.par1", "last.par2", "last.par.ok", "last.par.best",
                                  "value.best", "random.start")

sparse_joint_tmb_state <- function(env) {
  mget(intersect(sparse_joint_tmb_state_names, ls(env, all.names = TRUE)), envir = env)
}

sparse_joint_restore_tmb_state <- function(env, state) {
  for (name in names(state)) assign(name, state[[name]], envir = env)
  invisible(NULL)
}

# Every inner optimisation starts from `start`, by default the initial latent
# values.
sparse_joint_fix_inner_start <- function(env, start = env$par[env$random]) {
  assign(".alfak2_inner_start", start, envir = env)
  assign("random.start", expression(.alfak2_inner_start), envir = env)
  invisible(NULL)
}

# Smallest eigenvalue of the inner Hessian at `mode`, by inverse iteration on
# its sparse Cholesky factor: the Rayleigh quotient, which approaches it from
# above, fast when it is well separated from the next (as a near-zero one is).
# Zero if there is no Cholesky factor (not positive definite), NA if the
# Hessian cannot be computed. The start vector is fixed, so the result does
# not depend on, or change, the random number stream.
sparse_joint_inner_min_eigen <- function(env, mode = env$last.par, iterations = 200L, tol = 1e-6) {
  sparse_joint_inner_flattest(env, mode, iterations, tol)$value
}

# The same, with the direction: `value` as above and `vector`, the unit
# eigenvector of the inner Hessian's smallest eigenvalue (NULL if there is
# none). Without a Cholesky factor `value` is zero and `vector` points along
# the most negative eigenvalue, from inverse iteration on the Hessian shifted
# until it has one.
sparse_joint_inner_flattest <- function(env, mode = env$last.par, iterations = 200L, tol = 1e-6) {
  H <- tryCatch(env$spHess(mode, random = TRUE), error = function(e) NULL)
  # Not computable (for instance a depth offset run off to minus infinity):
  # NA, which the check reads as near-singular.
  if (is.null(H)) return(list(value = NA_real_, vector = NULL))
  entries <- if (methods::is(H, "CsparseMatrix")) H@x else as.numeric(as.matrix(H))
  if (any(!is.finite(entries))) return(list(value = NA_real_, vector = NULL))
  # CHOLMOD only warns when the matrix is not positive definite.
  cholesky <- function(M) tryCatch(Matrix::Cholesky(M, LDL = FALSE), warning = function(w) NULL,
                                   error = function(e) NULL)
  iterate <- function(factor) {
    x <- sin(seq_len(nrow(H)))
    x <- x / sqrt(sum(x^2))
    lambda <- Inf
    for (k in seq_len(iterations)) {
      y <- as.numeric(Matrix::solve(factor, x))
      if (any(!is.finite(y))) return(list(value = 0, vector = NULL))
      x <- y / sqrt(sum(y^2))
      rayleigh <- sum(x * as.numeric(H %*% x))
      if (!is.finite(rayleigh)) return(list(value = NA_real_, vector = NULL))
      if (abs(rayleigh - lambda) <= tol * abs(rayleigh)) return(list(value = rayleigh, vector = x))
      lambda <- rayleigh
    }
    list(value = lambda, vector = x)
  }
  factor <- cholesky(H)
  if (!is.null(factor)) return(iterate(factor))
  scale <- max(abs(Matrix::diag(H)), 1)
  for (shift in scale * 10^seq(-8, 2)) {
    factor <- cholesky(H + Matrix::Diagonal(nrow(H), shift))
    if (!is.null(factor)) return(list(value = 0, vector = iterate(factor)$vector))
  }
  list(value = 0, vector = NULL)
}

# What each latent variable of the inner problem is, to name the flat
# direction of a fit that is not accepted. `names` are their parameter names,
# in the order of the inner problem. `states` describes the states that
# `zeta` and `g` follow, one row each in their order (`karyotype`, `cells`:
# its counts in each sample, as text, `major`, and in a cohort `sample`);
# `nodes` likewise the karyotypes of a cohort's shared birth rates `m`.
sparse_joint_latent_labels <- function(names, states, nodes = NULL, state_index = NULL) {
  # `state_index`: the state of each latent variable, where it is not its
  # position among those of its name (the pairwise and state-space methods).
  index <- if (is.null(state_index)) as.integer(stats::ave(seq_along(names), names, FUN = seq_along)) else state_index
  variable <- c(zeta = "initial abundance", g = "birth rate", m = "birth rate", coef = "trend coefficient",
                theta1 = "depth", theta = "depth", zeta_later = "later abundance")[names]
  out <- data.frame(variable = ifelse(is.na(variable), names, variable), karyotype = NA_character_,
                    sample = NA_character_, cells = NA_character_, major = NA, stringsAsFactors = FALSE)
  fill <- function(out, rows, table) {
    if (is.null(table) || !any(rows)) return(out)
    for (column in intersect(c("karyotype", "sample", "cells", "major"), names(table))) {
      out[[column]][rows] <- table[[column]][index[rows]]
    }
    out
  }
  out <- fill(out, names %in% c("zeta", "g", "zeta_later"), states)
  fill(out, names == "m", nodes)
}

# The latent variables along a flattest direction `vector` that make up at
# least `min_share` of it (the largest at least, `max_rows` at most), with
# their `share`.
sparse_joint_flat_direction <- function(vector, labels, min_share = 0.05, max_rows = 5L) {
  if (is.null(vector) || is.null(labels) || length(vector) != nrow(labels)) return(NULL)
  share <- vector^2 / sum(vector^2)
  keep <- order(share, decreasing = TRUE)[seq_len(min(max_rows, max(1L, sum(share >= min_share))))]
  out <- labels[keep, , drop = FALSE]
  out$share <- share[keep]
  rownames(out) <- NULL
  out
}

# Shares of a flattest direction `vector` by kind of latent variable
# (`variable`, and `major` for those of a karyotype), largest first, those
# with at least `min_share` (the largest at least).
sparse_joint_flat_composition <- function(vector, labels, min_share = 0.01) {
  if (is.null(vector) || is.null(labels) || length(vector) != nrow(labels)) return(NULL)
  share <- vector^2 / sum(vector^2)
  key <- paste(labels$variable, labels$major)
  groups <- unique(key)
  first <- match(groups, key)
  out <- data.frame(variable = labels$variable[first], major = labels$major[first],
                    share = as.numeric(tapply(share, factor(key, levels = groups), sum)), stringsAsFactors = FALSE)
  out <- out[order(out$share, decreasing = TRUE), , drop = FALSE]
  out <- out[seq_len(max(1L, sum(out$share >= min_share))), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Words on a flat direction, from its leading latent variables `flat` (those
# with at least `min_share` of it, the largest at least, `max_terms` at most)
# and its `composition`: "in a direction made up of ..., led by ..." or
# "along ...".
sparse_joint_describe_flat <- function(flat, composition = NULL, min_share = 0.25, max_terms = 3L) {
  if (is.null(flat) || !nrow(flat)) return(NULL)
  leading <- function(table, max_rows) {
    table[seq_len(min(max_rows, max(1L, sum(table$share >= min_share)))), , drop = FALSE]
  }
  rows <- leading(flat, max_terms)
  lead <- vapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    what <- if (is.na(r$karyotype)) paste("a", r$variable) else {
      paste0("the ", r$variable, " of karyotype ", r$karyotype, if (!is.na(r$sample)) paste0(" in sample ", r$sample))
    }
    detail <- c(if (!is.na(r$major)) if (isTRUE(r$major)) "major" else "not major",
                if (!is.na(r$cells)) paste("cells", r$cells), sprintf("%.0f%%", 100 * r$share))
    paste0(what, " (", paste(detail, collapse = ", "), ")")
  }, character(1))
  lead <- paste(lead, collapse = " and ")
  if (is.null(composition) || !nrow(composition)) return(paste("along", lead))
  plural <- c("initial abundance" = "initial abundances", "birth rate" = "birth rates",
              "trend coefficient" = "trend coefficients", "depth" = "depths")
  groups <- leading(composition, 2L)
  groups <- vapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, ]
    name <- if (is.na(plural[g$variable])) g$variable else plural[[g$variable]]
    if (!is.na(g$major)) name <- paste(name, if (isTRUE(g$major)) "of major karyotypes" else "of karyotypes that are not major")
    sprintf("%s (%.0f%%)", name, 100 * g$share)
  }, character(1))
  paste0("in a direction made up of ", paste(groups, collapse = " and "), ", led by ", lead)
}

# Projected gradient: components that push against an active bound are zero.
sparse_joint_projected_gradient <- function(gradient, par, lower, upper) {
  lower <- rep_len(lower, length(par))
  upper <- rep_len(upper, length(par))
  gradient[par <= lower + 1e-8 & gradient > 0] <- 0
  gradient[par >= upper - 1e-8 & gradient < 0] <- 0
  gradient
}

# Central finite-difference gradient; one-sided at a bound.
sparse_joint_fd_gradient <- function(fn, par, lower, upper, step) {
  lower <- rep_len(lower, length(par))
  upper <- rep_len(upper, length(par))
  vapply(seq_along(par), function(i) {
    up <- par
    down <- par
    up[[i]] <- min(par[[i]] + step, upper[[i]])
    down[[i]] <- max(par[[i]] - step, lower[[i]])
    (fn(up) - fn(down)) / (up[[i]] - down[[i]])
  }, numeric(1))
}

# Whether `par`, with the objective `value` the optimiser reported there, is
# an optimum to accept. The fit's own inner mode is the one TMB reaches from
# its current state, as the posterior will. TMB's state is restored
# afterwards. When the point is not accepted, `flat_direction` lists the
# latent variables (described by `labels`, from sparse_joint_latent_labels())
# along the flattest direction of the inner Hessian at the fit's mode, and
# `flat_composition` its make-up by kind of latent variable.
sparse_joint_check_scales <- function(objective, par, value, lower, upper, tolerance, labels = NULL) {
  env <- objective$env
  state <- sparse_joint_tmb_state(env)
  on.exit(sparse_joint_restore_tmb_state(env, state), add = TRUE)
  evaluate <- function(p) tryCatch(as.numeric(objective$fn(p)), error = function(e) NA_real_)
  joint <- function() tryCatch(as.numeric(env$f(env$last.par, order = 0)), error = function(e) NA_real_)
  refit <- evaluate(par)
  mode <- env$last.par
  flattest <- sparse_joint_inner_flattest(env)
  eigen_fit <- flattest$value
  joint_fit <- joint()
  # The inner problem can have more than one mode: from the initial latent
  # values, is there a better one?
  sparse_joint_fix_inner_start(env)
  fixed <- evaluate(par)
  joint_fixed <- joint()
  # The gradient of the objective that follows the fit's mode: every step
  # starts the inner optimisation from it.
  sparse_joint_fix_inner_start(env, mode[env$random])
  gradient <- sparse_joint_fd_gradient(evaluate, par, lower, upper, tolerance[["step"]])
  projected <- max(abs(sparse_joint_projected_gradient(gradient, par, lower, upper)))
  reasons <- c(
    if (!is.finite(value) || !is.finite(refit)) "objective not finite",
    if (is.finite(value) && is.finite(refit) && abs(refit - value) > tolerance[["objective"]]) {
      "objective not reproduced at the fit's mode"
    },
    if (!isTRUE(eigen_fit >= tolerance[["eigen"]])) "near-singular inner Hessian",
    if (isTRUE(joint_fixed < joint_fit - tolerance[["objective"]])) "a better inner mode exists",
    if (any(!is.finite(gradient))) "objective not finite next to the point",
    if (all(is.finite(gradient)) && !isTRUE(projected <= tolerance[["gradient"]])) "gradient not small"
  )
  list(pass = !length(reasons), reasons = reasons, objective_difference = fixed - value,
       joint_difference = joint_fixed - joint_fit, min_eigen = eigen_fit, projected_gradient = projected,
       flat_direction = if (length(reasons)) sparse_joint_flat_direction(flattest$vector, labels),
       flat_composition = if (length(reasons)) sparse_joint_flat_composition(flattest$vector, labels))
}

# The objective with near-singular points rejected (NaN). Counts its
# evaluations and keeps the inner mode of the last one.
sparse_joint_guarded_objective <- function(objective, tolerance) {
  env <- objective$env
  memo <- new.env(parent = emptyenv())
  memo$evaluations <- 0L
  memo$mode <- NULL
  fn <- function(p) {
    memo$evaluations <- memo$evaluations + 1L
    f <- tryCatch(as.numeric(objective$fn(p)), error = function(e) NaN)
    memo$mode <- env$last.par
    if (!is.finite(f) || !isTRUE(sparse_joint_inner_min_eigen(env) >= tolerance[["eigen"]])) NaN else f
  }
  list(fn = fn, memo = memo)
}

# Bounded compass search: coordinate steps that halve when none improves.
# `improved_hook()` is called after each evaluation that improves the best
# point.
sparse_joint_compass_search <- function(fn, par, lower, upper, step = 0.25, min_step = 1e-3, max_eval = 150L,
                                        improved_hook = function() NULL) {
  lower <- rep_len(lower, length(par))
  upper <- rep_len(upper, length(par))
  value <- fn(par)
  if (is.finite(value)) improved_hook()
  evaluations <- 1L
  while (step >= min_step && evaluations < max_eval) {
    improved <- FALSE
    for (i in seq_along(par)) for (direction in c(-1, 1)) {
      candidate <- par
      candidate[[i]] <- min(max(par[[i]] + direction * step, lower[[i]]), upper[[i]])
      if (candidate[[i]] == par[[i]]) next
      f <- fn(candidate)
      evaluations <- evaluations + 1L
      if (is.finite(f) && (!is.finite(value) || f < value - 1e-10)) {
        par <- candidate
        value <- f
        improved <- TRUE
        improved_hook()
      }
    }
    if (!improved) step <- step / 2
  }
  list(par = par, objective = value, evaluations = evaluations,
       convergence = if (step < min_step) 0L else 1L,
       message = if (step < min_step) "step below its minimum" else "evaluation limit reached")
}

# Scale parameters by Laplace marginal likelihood, kept inside a finite box:
# outside it the inner problem can become degenerate. Returns nlminb's list
# with `convergence` 0 when the point was accepted and 1 when not, and
# `check`: how it was reached and checked. `labels` describe the latent
# variables (sparse_joint_latent_labels()), to name the flat direction of a
# point that is not accepted.
sparse_joint_optimise_scales <- function(objective, control, lower = log(1e-4), upper = log(0.5),
                                         tolerance = sparse_joint_convergence_tolerance, labels = NULL) {
  env <- objective$env
  first <- tryCatch(
    suppressWarnings(stats::nlminb(objective$par, objective$fn, objective$gr,
                                   lower = lower, upper = upper, control = control)),
    error = function(e) NULL
  )
  first_evaluations <- if (is.null(first)) NA_integer_ else as.integer(first$evaluations[[1L]])
  record <- function(optimum, stage, check, restarts, evaluations) {
    optimum$check <- list(
      status = if (check$pass) "converged" else "not converged", stage = stage,
      reasons = check$reasons, restarts = restarts, evaluations = evaluations,
      nlminb_convergence = if (is.null(first)) NA_integer_ else as.integer(first$convergence),
      nlminb_message = if (is.null(first)) "nlminb failed" else first$message,
      objective_difference = check$objective_difference, joint_difference = check$joint_difference,
      min_eigen = check$min_eigen, projected_gradient = check$projected_gradient,
      flat_direction = check$flat_direction, flat_composition = check$flat_composition
    )
    optimum$convergence <- if (check$pass) 0L else 1L
    optimum$message <- if (check$pass) {
      if (stage == "first") "converged" else paste("converged after restart:", stage)
    } else paste("not converged:", paste(check$reasons, collapse = "; "))
    optimum
  }
  check <- list(pass = FALSE, reasons = "nlminb failed", objective_difference = NA_real_,
                joint_difference = NA_real_, min_eigen = NA_real_, projected_gradient = NA_real_,
                flat_direction = NULL, flat_composition = NULL)
  if (!is.null(first)) {
    check <- sparse_joint_check_scales(objective, first$par, first$objective, lower, upper, tolerance, labels)
    if (check$pass) return(record(first, "first", check, 0L, first_evaluations))
    # Searching again is internal and off by default: on the development
    # benchmark it accepted few more points (77% of candidates against 73%)
    # at twice the time, and did not change accuracy.
    if (!isTRUE(getOption("alfak2.convergence_restarts", FALSE))) {
      return(record(first, "first", check, 0L, first_evaluations))
    }
  } else if (!isTRUE(getOption("alfak2.convergence_restarts", FALSE))) {
    # nlminb failed: the prior-centred scales, as before this check (alfak2
    # ff13cf4), checked for the record.
    value <- as.numeric(objective$fn(objective$par))
    fallback <- list(par = objective$par, objective = value, convergence = 99L,
                     message = "scale optimisation failed")
    check <- sparse_joint_check_scales(objective, objective$par, value, lower, upper, tolerance, labels)
    check$pass <- FALSE
    check$reasons <- c("nlminb failed", check$reasons)
    return(record(fallback, "none", check, 0L, first_evaluations))
  }

  # Search again: from the end of the first run, then from the prior-centred
  # start. TMB's best mode is restored afterwards.
  saved <- sparse_joint_tmb_state(env)
  sparse_joint_fix_inner_start(env)
  guarded <- sparse_joint_guarded_objective(objective, tolerance)
  starts <- Filter(Negate(is.null), list(end_of_first_run = if (!is.null(first)) first$par,
                                         prior_centred = objective$par))
  results <- list()
  accepted <- list()
  # Each search starts its inner optimisations from the initial latent values;
  # once it holds a healthy point, from that point's mode, which is far
  # cheaper. The check then asks whether the initial values lead elsewhere.
  for (stage in names(starts)) {
    sparse_joint_fix_inner_start(env)
    run <- sparse_joint_compass_search(guarded$fn, starts[[stage]], lower, upper, improved_hook = function() {
      sparse_joint_fix_inner_start(env, guarded$memo$mode[env$random])
    })
    if (!is.finite(run$objective)) next
    # The best point's own mode, from which the check and the fit start.
    run$mode <- env$.alfak2_inner_start
    run$check <- sparse_joint_check_scales(objective, run$par, run$objective, lower, upper, tolerance, labels)
    results[[stage]] <- run
    if (run$check$pass) {
      accepted[[stage]] <- run
      break
    }
  }
  evaluations <- sum(first_evaluations, guarded$memo$evaluations, na.rm = TRUE)
  pool <- if (length(accepted)) accepted else results
  if (!length(pool)) {
    # Nothing with a healthy inner Hessian: nlminb's own result, or the
    # prior-centred scales if nlminb failed, as before this check.
    sparse_joint_restore_tmb_state(env, saved)
    if (is.null(first)) {
      first <- list(par = objective$par, objective = as.numeric(objective$fn(objective$par)), convergence = 99L,
                    message = "scale optimisation failed")
    }
    check$reasons <- c("no point with a healthy inner Hessian", check$reasons)
    return(record(first, "none", check, length(results), evaluations))
  }
  stage <- names(pool)[[which.min(vapply(pool, function(run) run$objective, numeric(1)))]]
  best <- pool[[stage]]
  # Leave TMB at the chosen point's mode, with its usual inner start: the
  # posterior is evaluated from there and stays at that mode.
  sparse_joint_fix_inner_start(env, best$mode)
  value <- as.numeric(objective$fn(best$par))
  mode <- env$last.par
  sparse_joint_restore_tmb_state(env, saved)
  assign("last.par.best", mode, envir = env)
  assign("value.best", value, envir = env)
  best$objective <- value
  record(best, stage, best$check, length(results), evaluations)
}

# The candidate (history or trend) to keep: the lowest finite objective. With
# the internal option `alfak2.choose_accepted`, the lowest among those whose
# optimisation was accepted (all, if none was): one that was not accepted may
# owe its objective to a degenerate Laplace approximation. Off by default: it
# can also pass over a much better model whose optimisation was not accepted.
sparse_joint_choose_candidate <- function(objective, convergence) {
  usable <- is.finite(objective)
  pool <- which(usable)
  if (isTRUE(getOption("alfak2.choose_accepted", FALSE)) && any(usable & convergence == 0)) {
    pool <- which(usable & convergence == 0)
  }
  if (!length(pool)) return(NA_integer_)
  pool[[which.min(objective[pool])]]
}

# Warning for a fit whose scale optimisation was not accepted; `check` is its
# `diagnostics$optimisation`.
sparse_joint_convergence_warning <- function(check) {
  if (is.null(check) || identical(check$status, "converged")) return(invisible(NULL))
  flat <- if ("near-singular inner Hessian" %in% check$reasons) {
    sparse_joint_describe_flat(check$flat_direction, check$flat_composition)
  }
  warning(structure(class = c("alfak2_convergence_warning", "warning", "condition"), list(
    message = paste0("The scale optimisation did not converge (", paste(check$reasons, collapse = "; "),
                     "); the best point found is used.",
                     if (!is.null(flat)) paste0(" The inner Hessian is nearly flat ", flat,
                                                "; the data hardly constrain these."),
                     " See `diagnostics$optimisation`."),
    call = NULL
  )))
}
