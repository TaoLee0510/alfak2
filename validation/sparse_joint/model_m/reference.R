# Stage M (handoff HANDOFF_CLAUDE_2026-11-13_multi_timepoint.md): the sparse
# joint model with the initial abundance of non-major states integrated out.
#
# Functions only; sourced by m0.R. Two densities of the pair model
# (src/sparse_joint_model.h), as negative logs:
#
#   mm_joint_nll()     the current joint density of all latent variables, a
#                      line-by-line translation of the template, checked
#                      against TMB;
#   mm_marginal_nll()  the stage-M density: for each non-major state i the
#                      initial log abundance z_i is integrated out,
#
#     - log int N(z; m_i, s_i) Pois(y0_i; e^z)^w_i Pois(y1_i; lambda_i(z))^w_i dz,
#
#     lambda_i(z) = exp(c + log(e^(z + a_i) + e^(B_i))), a_i = net_i dt,
#     c = theta1 - net_ref dt, B_i the log inflow from i's neighbours. Under
#     the balance prior m_i is the balance value and s_i = balance_sd; a
#     non-major state that feeds a neighbour does so at its balance value
#     (the inflow is first order in the mutation rate, and a non-major
#     source carries at most a few per cent of any count).
#
# The integral is one-dimensional per state; mm_quadrature() holds the rules
# compared in M0.

mm_logsumexp <- function(x) {
  top <- max(x)
  if (!is.finite(top)) return(top)
  top + log(sum(exp(x - top)))
}

# Row-wise log(exp(a) + exp(b)) for vectors, -Inf allowed.
mm_logadd <- function(a, b) {
  top <- pmax(a, b)
  out <- top + log1p(exp(-abs(a - b)))
  out[is.infinite(a) & a < 0] <- b[is.infinite(a) & a < 0]
  out[is.infinite(b) & b < 0] <- a[is.infinite(b) & b < 0]
  out
}

# Grouped log-sum-exp: one value per level of `group` in 1..n, -Inf if empty.
mm_group_logsumexp <- function(x, group, n) {
  out <- rep(-Inf, n)
  if (!length(x)) return(out)
  top <- tapply(x, group, max)
  idx <- as.integer(names(top))
  sums <- tapply(exp(x - top[as.character(group)]), group, sum)
  out[idx] <- as.numeric(top) + log(as.numeric(sums[names(top)]))
  out
}

mm_log_phi1 <- function(a, b, dt) {
  gap <- (a - b) * dt
  x <- sqrt(gap * gap + 1e-12)
  0.5 * (a + b) * dt + 0.5 * x + log(dt) + log((1 - exp(-x)) / x)
}

# TMB's dpois, for non-integer counts as well.
mm_dpois <- function(y, log_lambda) y * log_lambda - exp(log_lambda) - lgamma(y + 1)

# Parameters as a list, mapped ones included, from a TMB objective and a
# full parameter vector (fixed and random, as env$last.par).
mm_unpack <- function(objective, par) lapply(objective$env$parList(par = par), as.numeric)

# Everything but the counts of non-major states, shared by both densities.
mm_parts <- function(p, d) {
  n <- length(d$y0)
  major <- seq_len(n) <= d$n_major
  scale <- exp(p$log_scale)
  prior <- -sum(stats::dnorm(p$log_scale, d$log_scale_prior_mean, d$log_scale_prior_sd, log = TRUE))
  group <- d$coef_group
  coef_sd <- ifelse(group == 0, d$intercept_sd, ifelse(group == 4, d$shell_sd, scale[pmax(group, 1L)]))
  coef_mean <- ifelse(group == 0, d$intercept_mean, 0)
  prior <- prior - sum(stats::dnorm(p$coef, coef_mean, coef_sd, log = TRUE))
  trend <- drop(d$design %*% p$coef)
  prior <- prior - sum(stats::dnorm(p$g, trend, scale[[4L]], log = TRUE))
  net <- (2 * exp(d$log_self) - 1) * exp(p$g)
  from <- d$edge_from + 1L
  to <- d$edge_to + 1L
  # Balance value of each non-major state fed by a major.
  log_start <- rep(-Inf, n)
  fed <- rep(FALSE, n)
  if (d$balance == 1L) {
    e <- which(major[from] & !major[to])
    term <- log(2) + d$edge_log_prob[e] + p$g[from[e]] + p$zeta[from[e]] - net[from[e]] * d$history +
      mm_log_phi1(net[to[e]], net[from[e]], d$history)
    log_start <- mm_group_logsumexp(term, to[e], n)
    fed[unique(to[e])] <- TRUE
  }
  list(n = n, major = major, prior = prior, net = net, from = from, to = to, log_start = log_start, fed = fed,
       centre = net[[d$reference + 1L]] * d$dt)
}

mm_joint_nll <- function(p, d) {
  q <- mm_parts(p, d)
  nll <- q$prior
  if (d$balance == 1L) {
    nll <- nll - sum(stats::dnorm(p$zeta[q$fed] - q$log_start[q$fed] - d$balance_offset[q$fed], 0, d$balance_sd,
                                  log = TRUE))
    nll <- nll - sum(stats::dnorm(p$zeta[!q$fed], d$zeta_prior_mean[!q$fed], d$zeta_prior_sd[!q$fed], log = TRUE))
  } else {
    nll <- nll - sum(stats::dnorm(p$zeta, d$zeta_prior_mean, d$zeta_prior_sd, log = TRUE))
  }
  inflow <- log(2) + d$edge_log_prob + p$g[q$from] + p$zeta[q$from] +
    mm_log_phi1(q$net[q$to], q$net[q$from], d$dt)
  log_n1 <- mm_logadd(p$zeta + q$net * d$dt, mm_group_logsumexp(inflow, q$to, q$n))
  log_intensity1 <- log_n1 - q$centre + p$theta1
  w <- d$weight
  keep <- w != 0
  nll - sum(w[keep] * (mm_dpois(d$y0[keep], p$zeta[keep]) + mm_dpois(d$y1[keep], log_intensity1[keep])))
}

# Quadrature for log E[exp(l(m + s u))], u ~ N(0, 1), for K states at once.
# `l(z)` takes a K-row matrix of points and returns the log integrand there.
#   gh:      Gauss-Hermite on the prior, `nodes` points;
#   grid:    trapezoid on the prior, u in [-width, width] with step `step`;
#   mixture: Gauss-Hermite on the prior and on a normal centred where the
#            counts alone would put z (`centre`, `spread`), combined through
#            their mixture density;
#   exact:   a trapezoid fine enough to serve as the truth (checked against
#            stats::integrate in M0).
mm_gh <- function(k) {
  # Golub-Welsch for the probabilists' Hermite polynomials.
  if (k == 1L) return(list(x = 0, w = 1))
  J <- matrix(0, k, k)
  off <- sqrt(seq_len(k - 1L))
  J[cbind(1:(k - 1L), 2:k)] <- off
  J[cbind(2:k, 1:(k - 1L))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1L, ]^2)
}

mm_row_logsumexp <- function(x) {
  top <- apply(x, 1L, max)
  top + log(rowSums(exp(x - top)))
}

mm_quadrature <- function(l, m, s, rule = list(type = "gh", nodes = 20L), centre = NULL, spread = NULL) {
  K <- length(m)
  if (!K) return(numeric())
  if (rule$type %in% c("gh", "grid", "exact")) {
    if (rule$type == "gh") {
      q <- mm_gh(rule$nodes)
      u <- q$x
      lw <- log(q$w)
    } else {
      step <- if (rule$type == "exact") 2e-3 else rule$step
      width <- if (rule$type == "exact") 40 else rule$width
      u <- seq(-width, width, by = step)
      lw <- stats::dnorm(u, log = TRUE) + log(step)
    }
    if (rule$type == "exact" && K > 50L) {
      # In blocks, to bound memory.
      blocks <- split(seq_len(K), ceiling(seq_len(K) / 50L))
      return(unlist(lapply(blocks, function(b) {
        mm_row_logsumexp(sweep(l(outer(m[b], u * 0, "+") + outer(s[b], u), b), 2L, lw, "+"))
      }), use.names = FALSE))
    }
    z <- outer(m, u * 0, "+") + outer(s, u)
    return(mm_row_logsumexp(sweep(l(z, seq_len(K)), 2L, lw, "+")))
  }
  if (rule$type == "mixture") {
    q <- mm_gh(rule$nodes)
    z <- cbind(outer(m, q$x * 0, "+") + outer(s, q$x), outer(centre, q$x * 0, "+") + outer(spread, q$x))
    M <- matrix(m, K, ncol(z))
    S <- matrix(s, K, ncol(z))
    C <- matrix(centre, K, ncol(z))
    P <- matrix(spread, K, ncol(z))
    lprior <- stats::dnorm(z, M, S, log = TRUE)
    lq <- mm_logadd(log(0.5) + lprior, log(0.5) + stats::dnorm(z, C, P, log = TRUE))
    x <- sweep(lprior + l(z, seq_len(K)) - lq, 2L, log(0.5) + log(c(q$w, q$w)), "+")
    return(mm_row_logsumexp(x))
  }
  stop("unknown rule")
}

# Per-state inputs of the integral at parameters `p` (zeta of non-major
# states, if present, are ignored).
mm_state_inputs <- function(p, d) {
  q <- mm_parts(p, d)
  minor <- which(!q$major)
  zeta_source <- p$zeta
  if (d$balance == 1L) {
    zeta_source[minor] <- ifelse(q$fed[minor], q$log_start[minor], d$zeta_prior_mean[minor])
  } else {
    zeta_source[minor] <- d$zeta_prior_mean[minor]
  }
  inflow <- log(2) + d$edge_log_prob + p$g[q$from] + zeta_source[q$from] +
    mm_log_phi1(q$net[q$to], q$net[q$from], d$dt)
  log_inflow <- mm_group_logsumexp(inflow, q$to, q$n)
  fed <- d$balance == 1L & q$fed
  m <- ifelse(fed, q$log_start + d$balance_offset, d$zeta_prior_mean)
  s <- ifelse(fed, d$balance_sd, d$zeta_prior_sd)
  c1 <- p$theta1 - q$centre
  # Where the counts alone put z: the count at the first time point, and the
  # count at the second read as the state's own growth; the prior where they say nothing.
  own <- pmax(d$y1 - exp(c1 + log_inflow), 0)
  info <- d$y0 + own
  xlogx <- function(x) ifelse(x > 0, x * log(x), 0)
  centre <- ifelse(info > 0, (xlogx(d$y0) + xlogx(own) - own * (c1 + q$net * d$dt)) / pmax(info, 1e-300), m)
  spread <- ifelse(info > 0, 1 / sqrt(pmax(info, 1e-300)), s)
  list(q = q, zeta_source = zeta_source, log_inflow = log_inflow, a = q$net * d$dt, c1 = c1, m = m, s = s,
       centre = centre, spread = spread, minor = minor)
}

mm_state_loglik <- function(inputs, d) {
  force(inputs)
  # z: one row per non-major state in `rows` (indices into inputs$minor).
  function(z, rows) {
    i <- inputs$minor[rows]
    w <- d$weight[i]
    log_l1 <- inputs$c1 + mm_logadd(z + inputs$a[i], matrix(inputs$log_inflow[i], nrow(z), ncol(z)))
    w * (d$y0[i] * z - exp(z) - lgamma(d$y0[i] + 1) + d$y1[i] * log_l1 - exp(log_l1) - lgamma(d$y1[i] + 1))
  }
}

# Stage-M negative log density; `parts = TRUE` returns the terms per state.
mm_marginal_nll <- function(p, d, rule = list(type = "gh", nodes = 20L), parts = FALSE) {
  inputs <- mm_state_inputs(p, d)
  q <- inputs$q
  major <- which(q$major)
  nll <- q$prior
  if (d$balance == 1L) {
    nll <- nll - sum(stats::dnorm(p$zeta[major], d$zeta_prior_mean[major], d$zeta_prior_sd[major], log = TRUE))
  } else {
    nll <- nll - sum(stats::dnorm(p$zeta[major], d$zeta_prior_mean[major], d$zeta_prior_sd[major], log = TRUE))
  }
  log_n1 <- mm_logadd(p$zeta[major] + q$net[major] * d$dt, inputs$log_inflow[major])
  w <- d$weight[major]
  nll <- nll - sum(w * (mm_dpois(d$y0[major], p$zeta[major]) + mm_dpois(d$y1[major], inputs$c1 + log_n1)))
  minor <- inputs$minor
  terms <- -mm_quadrature(mm_state_loglik(inputs, d), inputs$m[minor], inputs$s[minor], rule,
                          inputs$centre[minor], inputs$spread[minor])
  if (parts) return(list(base = nll, terms = terms, inputs = inputs))
  nll + sum(terms)
}

# --- Capturing fits -------------------------------------------------------

# Fit `counts` with the loaded build, one candidate history per call to the
# scale optimiser, and keep each candidate's problem, scales and mode.
mm_capture <- function(fit_call, envir = parent.frame()) {
  grab <- new.env()
  grab$runs <- list()
  ns <- asNamespace("alfak2")
  suppressMessages(trace("sparse_joint_problem", where = ns, print = FALSE,
                         exit = quote(assign("problem", returnValue(), envir = .mm_grab))))
  suppressMessages(trace("sparse_joint_optimise_scales", where = ns, print = FALSE,
                         exit = quote(.mm_grab$runs[[length(.mm_grab$runs) + 1L]] <-
                                        list(problem = .mm_grab$problem, objective = objective,
                                             result = returnValue(),
                                             mode = objective$env$last.par.best + 0,
                                             value = objective$env$value.best + 0))))
  on.exit({
    suppressMessages(untrace("sparse_joint_problem", where = ns))
    suppressMessages(untrace("sparse_joint_optimise_scales", where = ns))
  }, add = TRUE)
  assign(".mm_grab", grab, envir = globalenv())
  started <- proc.time()[["elapsed"]]
  fit <- withCallingHandlers(eval(fit_call, envir), warning = function(w) invokeRestart("muffleWarning"))
  seconds <- proc.time()[["elapsed"]] - started
  # The mode the optimiser ended on, as it saw it: evaluating the objective
  # again after the fit can land the inner optimisation elsewhere.
  lapply(grab$runs, function(run) {
    env <- run$objective$env
    par <- run$result$par
    mode <- run$mode
    value <- as.numeric(run$value)
    problem <- run$problem
    list(history = problem$tmb_data$history, tmb_data = problem$tmb_data, parameters = problem$parameters,
         map = problem$map, random = problem$random, par = par, mode = mode, objective = value,
         mode_matches_par = max(abs(mode[-env$random] - par)),
         result_objective = as.numeric(run$result$objective),
         convergence = run$result$convergence, message = run$result$message, check = run$result$check,
         n_major = problem$tmb_data$n_major, n_states = length(problem$tmb_data$y0), seconds = seconds)
  })
}

mm_make_objective <- function(candidate, data = candidate$tmb_data) {
  TMB::MakeADFun(data, candidate$parameters, map = candidate$map, random = candidate$random, DLL = "alfak2",
                 silent = TRUE)
}
