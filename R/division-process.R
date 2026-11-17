# Internal single-daughter kernel for division dynamics. Its row sums may be
# below one; no survival or omitted-state mass is redistributed to other nodes.
division_transition_data <- function(graph, beta = graph$beta,
                                     probability_cutoff = 1e-14) {
  op <- alfak2_transition_operator_cpp(graph$karyotypes, beta, "exact",
                                      division_process = TRUE,
                                      probability_cutoff = probability_cutoff)
  cn <- sort(unique(as.integer(graph$karyotypes)))
  viability <- vapply(cn, function(parent) {
    sum(vapply(seq.int(graph$min_cn, graph$max_cn), function(child)
      alfak2_pij_cpp(parent, child, beta), numeric(1)))
  }, numeric(1))
  names(viability) <- as.character(cn)
  viable <- apply(graph$karyotypes, 1L, function(x) prod(viability[as.character(x)]))
  op$viable <- as.numeric(viable)
  op$biological_loss <- pmax(0, 1 - op$viable)
  op$computational_outflow <- pmax(0, op$viable - op$retained)
  op$leakage <- op$computational_outflow
  op
}

# Numerically independent R implementation for forward validation and prediction.
division_forward_r <- function(pi0, birth_rate, operator, dt, steps = max(4L, ceiling(dt / 0.1))) {
  if (length(pi0) != length(birth_rate) || any(!is.finite(birth_rate)) || any(birth_rate < 0))
    stop("Birth rates must be finite, non-negative, and aligned with pi0.", call. = FALSE)
  validate_scalar(dt, "dt", lower = 0)
  steps <- validate_grf_integer(steps, "steps", lower = 1L)
  from <- operator$from0 + 1L
  to <- operator$to0 + 1L
  same <- from == to
  n <- length(pi0)
  diagonal <- -birth_rate
  diagonal[from[same]] <- diagonal[from[same]] + 2 * operator$weight[same] * birth_rate[from[same]]
  shift <- max(diagonal)
  h <- dt / steps
  half <- exp(0.5 * h * (diagonal - shift))
  mass <- as.numeric(pi0) / sum(pi0)
  log_total <- 0
  off <- which(!same)
  for (s in seq_len(steps)) {
    mass <- mass * half
    z <- sum(mass); log_total <- log_total + log(z) + 0.5 * h * shift; mass <- mass / z
    term <- accum <- mass
    for (order in 1:4) {
      values <- 2 * operator$weight[off] * birth_rate[from[off]] * term[from[off]] * h / order
      next_term <- numeric(n)
      if (length(values)) {
        grouped <- rowsum(values, to[off], reorder = FALSE)
        next_term[as.integer(rownames(grouped))] <- grouped[, 1L]
      }
      accum <- accum + next_term
      term <- next_term
    }
    z <- sum(accum); log_total <- log_total + log(z); mass <- accum / z
    mass <- mass * half
    z <- sum(mass); log_total <- log_total + log(z) + 0.5 * h * shift; mass <- mass / z
  }
  list(pi1 = mass, log_total_growth = log_total)
}

division_forward_accuracy <- function(posterior, graph, dt, tolerance = 1e-5) {
  validate_scalar(tolerance, "forward_tolerance", lower = .Machine$double.eps)
  op <- division_transition_data(graph)
  pi0 <- posterior$posterior_predictive$pi0
  rates <- exp(posterior$latent_fitness_mean)
  steps <- posterior$diagnostics$forward_steps
  coarse <- division_forward_r(pi0, rates, op, dt, steps)
  fine <- division_forward_r(pi0, rates, op, dt, 2L * steps)
  frequency_error <- max(abs(coarse$pi1 - fine$pi1))
  growth_error <- abs(coarse$log_total_growth - fine$log_total_growth)
  backend_error <- max(abs(coarse$pi1 - posterior$posterior_predictive$pi1))
  list(status = if (backend_error > 1e-8) "backend_disagreement" else
         if (max(frequency_error, growth_error) <= tolerance) "passed_step_halving" else "refinement_required",
       tolerance = tolerance, steps = steps, max_frequency_change = frequency_error,
       cpp_r_frequency_agreement = backend_error,
       log_total_growth_change = growth_error,
       initial_weighted_graph_outflow = sum(pi0 * op$computational_outflow),
       maximum_single_division_graph_outflow = max(op$computational_outflow),
       maximum_single_division_biological_loss = max(op$biological_loss),
       maximum_fitted_birth_rate = max(rates))
}
