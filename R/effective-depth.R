validate_effective_depth_mode <- function(input_depth, effective_depth_mode) {
  input_depth <- match.arg(input_depth, c("raw", "effective"))
  effective_depth_mode <- match.arg(effective_depth_mode, c("min", "cap", "fixed"))
  list(input_depth = input_depth, effective_depth_mode = effective_depth_mode)
}

round_expected_counts_to_total <- function(expected, total) {
  total <- as.integer(round(total))
  if (!is.finite(total) || total < 1L) {
    stop("Effective depth totals must round to positive integers.", call. = FALSE)
  }
  expected[!is.finite(expected) | expected < 0] <- 0
  floored <- floor(expected)
  remainder <- total - sum(floored)
  if (remainder > 0L) {
    frac <- expected - floored
    ord <- order(frac, decreasing = TRUE)
    floored[ord[seq_len(min(remainder, length(ord)))]] <-
      floored[ord[seq_len(min(remainder, length(ord)))]] + 1L
  }
  as.integer(floored)
}

resolve_effective_depth <- function(raw_depth, effective_depth, effective_depth_mode) {
  if (any(!is.finite(raw_depth)) || any(raw_depth <= 0)) {
    stop("`input_depth = \"effective\"` requires every timepoint to have positive total counts.", call. = FALSE)
  }
  if (effective_depth_mode == "min") {
    out <- rep(min(raw_depth), length(raw_depth))
  } else {
    validate_scalar(effective_depth, "effective_depth", lower = 1)
    out <- if (effective_depth_mode == "fixed") {
      rep(effective_depth, length(raw_depth))
    } else {
      pmin(raw_depth, effective_depth)
    }
  }
  as.integer(pmax(1L, round(out)))
}

apply_effective_depth_counts <- function(counts,
                                         effective_depth = NULL,
                                         effective_depth_mode = c("min", "cap", "fixed")) {
  effective_depth_mode <- match.arg(effective_depth_mode)
  observation_weights <- attr(counts, "observation_weights", exact = TRUE)
  soft_minobs <- attr(counts, "soft_minobs", exact = TRUE)
  counts <- validate_count_matrix(counts)
  raw_depth <- colSums(counts)
  target_depth <- resolve_effective_depth(raw_depth, effective_depth, effective_depth_mode)
  freq <- sweep(counts, 2L, raw_depth, "/")
  out <- matrix(0L, nrow = nrow(counts), ncol = ncol(counts), dimnames = dimnames(counts))
  for (j in seq_len(ncol(counts))) {
    out[, j] <- round_expected_counts_to_total(freq[, j] * target_depth[j], target_depth[j])
  }
  keep <- rowSums(out) > 0L
  if (!any(keep)) {
    stop("Effective-depth preprocessing removed all karyotypes.", call. = FALSE)
  }
  out <- out[keep, , drop = FALSE]
  attr(out, "effective_depth_info") <- list(
    input_depth = "effective",
    effective_depth_mode = effective_depth_mode,
    requested_effective_depth = effective_depth,
    raw_depth = as.numeric(raw_depth),
    effective_depth = as.numeric(target_depth),
    retained_karyotypes = nrow(out),
    dropped_karyotypes = sum(!keep)
  )
  if (!is.null(observation_weights)) {
    attr(out, "observation_weights") <- subset_observation_weights(observation_weights, rownames(out))
  }
  if (!is.null(soft_minobs)) attr(out, "soft_minobs") <- soft_minobs
  out
}

prepare_counts_for_input_depth <- function(counts,
                                           dt,
                                           beta,
                                           input_depth = c("raw", "effective"),
                                           effective_depth = NULL,
                                           effective_depth_mode = c("min", "cap", "fixed")) {
  modes <- validate_effective_depth_mode(input_depth, effective_depth_mode)
  input_depth <- modes$input_depth
  effective_depth_mode <- modes$effective_depth_mode
  if (inherits(counts, "alfak2_data")) {
    data <- counts
    if (input_depth == "raw") return(data)
    if (!is.null(data$metadata$observation_weights)) {
      attr(data$counts, "observation_weights") <- data$metadata$observation_weights
    }
    if (!is.null(data$metadata$soft_minobs)) {
      attr(data$counts, "soft_minobs") <- data$metadata$soft_minobs
    }
    eff_counts <- apply_effective_depth_counts(data$counts, effective_depth, effective_depth_mode)
    info <- attr(eff_counts, "effective_depth_info")
    metadata <- data$metadata
    metadata$input_depth <- info
    return(prepare_alfak2_data(eff_counts, dt = data$dt, beta = data$beta, metadata = metadata))
  }
  if (input_depth == "raw") {
    return(prepare_alfak2_data(counts, dt = dt, beta = beta))
  }
  eff_counts <- apply_effective_depth_counts(counts, effective_depth, effective_depth_mode)
  info <- attr(eff_counts, "effective_depth_info")
  prepare_alfak2_data(eff_counts, dt = dt, beta = beta, metadata = list(input_depth = info))
}

resolve_fit_observation_controls <- function(input_depth, observation_model, dm_concentration) {
  if (is.null(observation_model)) {
    observation_model <- if (identical(input_depth, "effective")) {
      "dirichlet_multinomial"
    } else {
      "multinomial"
    }
  }
  observation_model <- match_observation_model(observation_model)
  if (is.null(dm_concentration)) {
    dm_concentration <- if (identical(input_depth, "effective") &&
                            identical(observation_model, "dirichlet_multinomial")) {
      50
    } else {
      200
    }
  }
  validate_scalar(dm_concentration, "dm_concentration", lower = .Machine$double.eps)
  list(observation_model = observation_model, dm_concentration = dm_concentration)
}
