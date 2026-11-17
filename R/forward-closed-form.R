# Closed-form forward model for division dynamics on a sparse karyotype domain.
#
# Mean dynamics of the division process are dn/dt = A n with
# A = (2 P' - I) diag(b): P is the single-daughter kernel and b the birth
# rates. Write A = diag(d) + M, where d_i = (2 P_ii - 1) b_i is the net growth
# rate of state i (losses to every other state included) and M_ji = 2 P_ij b_i
# is the flux from i to its neighbour j. Expanding in M with exp(diag(d) t)
# kept exact (a Dyson series) gives, to first order,
#
#   n_j(T) = exp(d_j T) n_j(0) + sum_i M_ji phi1(d_j, d_i; T) n_i(0),
#   phi1(a, b; T) = (exp(a T) - exp(b T)) / (a - b).
#
# The expansion parameter is the per-edge flux M T, which is small for
# realistic missegregation rates, while differences in d times T may be large:
# a deleterious neighbour is held at mutation-selection balance with its
# parent, and that is exactly what phi1 describes. A state depends only on
# itself and its neighbours, so the model stays sparse on large domains.
#
# These R functions are the reference implementation. They work in log space
# and never form exp(d T).

# One-step edges of the exact single-daughter kernel inside a domain.
# `karyotypes` is an integer matrix (one row per state). Simultaneous changes
# of two or more chromosomes (probability of order beta^2) are not edges;
# their mass, like mass leaving the domain, is part of the loss in `self`.
sparse_division_kernel <- function(karyotypes, beta) {
  karyotypes <- as.matrix(karyotypes)
  storage.mode(karyotypes) <- "integer"
  n <- nrow(karyotypes)
  n_chr <- ncol(karyotypes)
  max_cn <- max(karyotypes) + 1L
  # chromosome-level probabilities: stay, gain one copy, lose one copy
  copies <- 0:max_cn
  stay <- vapply(copies, function(cn) alfak2_pij_cpp(cn, cn, beta), numeric(1))
  gain <- vapply(copies, function(cn) alfak2_pij_cpp(cn, cn + 1L, beta), numeric(1))
  lose <- vapply(copies, function(cn) if (cn > 0L) alfak2_pij_cpp(cn, cn - 1L, beta) else 0, numeric(1))
  log_stay <- matrix(log(stay)[karyotypes + 1L], nrow = n)
  log_self <- rowSums(log_stay)
  labels <- format_karyotypes(karyotypes)
  index_of <- stats::setNames(seq_len(n), labels)
  from <- to <- integer()
  log_prob <- numeric()
  for (chr in seq_len(n_chr)) {
    for (direction in c(-1L, 1L)) {
      moved <- karyotypes
      moved[, chr] <- moved[, chr] + direction
      child <- unname(index_of[format_karyotypes(moved)])
      keep <- which(!is.na(child))
      if (!length(keep)) next
      step <- if (direction > 0L) gain[karyotypes[keep, chr] + 1L] else lose[karyotypes[keep, chr] + 1L]
      positive <- step > 0
      keep <- keep[positive]
      from <- c(from, keep)
      to <- c(to, child[!is.na(child)][positive])
      log_prob <- c(log_prob, log(step[positive]) + log_self[keep] - log_stay[keep, chr])
    }
  }
  list(n = n, from = from, to = to, log_prob = log_prob, log_self = log_self)
}

# log of phi1(a, b; T) = integral_0^T exp(a (T - s) + b s) ds, stable for any
# a, b including a == b.
log_phi1 <- function(a, b, dt) {
  gap <- abs(a - b) * dt
  ratio <- ifelse(gap < 1e-6, 1 - gap / 2, -expm1(-gap) / pmax(gap, .Machine$double.xmin))
  pmax(a, b) * dt + log(dt) + log(ratio)
}

# Normalised second divided difference: with u <= v non-negative,
# g(u, v) = integral over the simplex of exp(-l2 u - l3 v), so that
# phi2(x1 >= x2 >= x3; T) = exp(x1 T) T^2 g((x1 - x2) T, (x1 - x3) T).
dyson_g2 <- function(u, v) {
  h <- function(z) ifelse(z < 1e-6, 1 - z / 2 + z^2 / 6, -expm1(-z) / pmax(z, .Machine$double.xmin))
  # -h'(z), the limit of (h(u) - h(v)) / (v - u) as v -> u
  minus_dh <- function(z) ifelse(z < 1e-4, 0.5 - z / 3 + z^2 / 8, (1 - exp(-z) * (1 + z)) / pmax(z, .Machine$double.xmin)^2)
  close <- (v - u) < 1e-6 * pmax(1, v)
  ifelse(close, minus_dh((u + v) / 2), (h(u) - h(v)) / pmax(v - u, .Machine$double.xmin))
}

# log of phi2(a, b, c; T), the time-ordered double integral of a two-step path.
log_phi2 <- function(a, b, c, dt) {
  top <- pmax(a, b, c)
  low <- pmin(a, b, c)
  middle <- a + b + c - top - low
  top * dt + 2 * log(dt) + log(dyson_g2((top - middle) * dt, (top - low) * dt))
}

row_log_sum_exp <- function(value, group, n) {
  peak <- rep(-Inf, n)
  peak_by_group <- tapply(value, group, max)
  peak[as.integer(names(peak_by_group))] <- peak_by_group
  total <- numeric(n)
  sum_by_group <- tapply(exp(value - peak[group]), group, sum)
  total[as.integer(names(sum_by_group))] <- sum_by_group
  ifelse(is.finite(peak), peak + log(total), -Inf)
}

# log n(T), unnormalised, from log n(0) and log birth rates.
# `order = 1` keeps one-step paths; `order = 2` adds every two-step path.
closed_form_forward_r <- function(log_n0, log_birth, kernel, dt, order = 1L) {
  n <- kernel$n
  birth <- exp(log_birth)
  net <- (2 * exp(kernel$log_self) - 1) * birth
  log_flux <- log(2) + kernel$log_prob + log_birth[kernel$from]
  terms <- c(log_n0 + net * dt,
             log_flux + log_n0[kernel$from] + log_phi1(net[kernel$to], net[kernel$from], dt))
  groups <- c(seq_len(n), kernel$to)
  if (order >= 2L) {
    # two-step paths i -> k -> j: join edges on the intermediate state k
    first <- split(seq_along(kernel$from), kernel$to)   # edges ending in k
    second <- split(seq_along(kernel$from), kernel$from) # edges starting in k
    shared <- intersect(names(first), names(second))
    path <- do.call(rbind, lapply(shared, function(k) {
      as.matrix(expand.grid(e1 = first[[k]], e2 = second[[k]]))
    }))
    i <- kernel$from[path[, "e1"]]
    k <- kernel$to[path[, "e1"]]
    j <- kernel$to[path[, "e2"]]
    terms <- c(terms, log_flux[path[, "e1"]] + log_flux[path[, "e2"]] + log_n0[i] +
                 log_phi2(net[j], net[k], net[i], dt))
    groups <- c(groups, j)
  }
  row_log_sum_exp(terms, groups, n)
}

# Reference solution exp(A T) n(0) on the same domain by uniformisation, with
# rates shifted so that every term is non-negative. Returns log n(T) up to the
# same additive constant convention as `closed_form_forward_r` (none dropped).
exact_forward_reference <- function(log_n0, log_birth, kernel, dt, tolerance = 1e-13) {
  n <- kernel$n
  birth <- exp(log_birth)
  net <- (2 * exp(kernel$log_self) - 1) * birth
  shift <- max(net)
  rate <- max(shift - net) + .Machine$double.eps
  flux <- 2 * exp(kernel$log_prob) * birth[kernel$from]
  # Uniformised kernel U = I + (A - shift I) / rate has no negative entry.
  operator <- Matrix::sparseMatrix(
    i = c(seq_len(n), kernel$to), j = c(seq_len(n), kernel$from),
    x = c(1 + (net - shift) / rate, flux / rate), dims = c(n, n)
  )
  offset <- max(log_n0)
  term <- exp(log_n0 - offset)
  mean_jumps <- rate * dt
  n_terms <- as.integer(ceiling(mean_jumps + 12 * sqrt(mean_jumps) + 30))
  weight <- stats::dpois(0:n_terms, mean_jumps)
  total <- weight[1L] * term
  for (jump in seq_len(n_terms)) {
    term <- as.numeric(operator %*% term)
    total <- total + weight[jump + 1L] * term
    if (jump > mean_jumps && weight[jump + 1L] * sum(term) <= tolerance * sum(total)) break
  }
  log(total) + offset + shift * dt
}
