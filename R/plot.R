#' Plot alfak2 posterior summaries
#'
#' @param x An `alfak2_fit` or `alfak2_local_fit`.
#' @param layer `"global"` or `"local"` for full fits.
#' @param ... Reserved.
#'
#' @return A `ggplot` object.
#' @export
plot_alfak2 <- function(x, layer = c("global", "local"), ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package `ggplot2` is required for plotting.", call. = FALSE)
  }
  layer <- match.arg(layer)
  s <- summarize_alfak2(x, layer = layer)
  s$node_order <- seq_len(nrow(s))
  ggplot2::ggplot(s, ggplot2::aes(.data$node_order, .data$fitness_mean, color = .data$support_tier)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, color = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$conf_low, ymax = .data$conf_high),
                           width = 0, alpha = 0.45) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::labs(x = "Graph node", y = "Posterior fitness", color = "Support tier") +
    ggplot2::theme_bw()
}

#' @export
plot.alfak2_fit <- function(x, ...) {
  plot_alfak2(x, ...)
}

#' @export
plot.alfak2_local_fit <- function(x, ...) {
  plot_alfak2(x, layer = "local", ...)
}
