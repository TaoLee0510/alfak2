# Draw from the actual model mixture, retaining correlations before taking any
# fitness reference contrast. Components share one fitness interpretation.
draw_alfak2_posterior_components <- function(components, target_labels, n_draws = 8000L,
                                             seed = 711204L, default_labels = NULL) {
  if (!length(components)) stop("No posterior components")
  weights <- vapply(components, function(x) as.numeric(x$weight), numeric(1))
  if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) stop("Invalid mixture weights")
  forward <- unique(vapply(components, function(x) as.character(x$forward_model), character(1)))
  if (length(forward) != 1L) stop("Mixture components have incompatible fitness interpretations")
  set.seed(seed)
  selected <- sample.int(length(components), n_draws, replace = TRUE, prob = weights / sum(weights))
  draws <- matrix(NA_real_, length(target_labels), n_draws, dimnames = list(target_labels, NULL))
  for (i in seq_along(components)) {
    columns <- which(selected == i)
    if (!length(columns)) next
    component <- components[[i]]
    division <- identical(component$forward_model, "division_process")
    mean <- if (division) component$latent_mean else component$mean
    covariance <- if (division) component$latent_covariance else component$covariance
    if (is.null(mean) || is.null(covariance)) stop("Mixture component is missing its full latent distribution")
    covariance <- as.matrix(covariance)
    labels <- names(mean)
    if (is.null(labels)) labels <- rownames(covariance)
    if (is.null(labels)) labels <- component$labels
    if (is.null(labels)) labels <- default_labels
    ii <- match(target_labels, labels)
    if (anyNA(ii)) stop("Mixture component does not cover every target node")
    mean <- as.numeric(mean)[ii]
    covariance <- covariance[ii, ii, drop = FALSE]
    eig <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
    if (any(!is.finite(eig$values)) || min(eig$values) < -1e-7) stop("Invalid mixture component covariance")
    transform <- sweep(eig$vectors, 2L, sqrt(pmax(0, eig$values)), "*")
    z <- matrix(stats::rnorm(length(mean) * length(columns)), nrow = length(mean))
    value <- sweep(transform %*% z, 1L, mean, "+")
    if (division) value <- exp(value)
    draws[, columns] <- value
  }
  if (any(!is.finite(draws))) stop("Non-finite posterior mixture draws")
  draws
}
