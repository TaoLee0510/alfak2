new_alfak2_data <- function(counts, dt, beta, metadata = list()) {
  structure(
    list(
      counts = counts,
      labels = rownames(counts),
      dt = dt,
      beta = beta,
      metadata = metadata
    ),
    class = "alfak2_data"
  )
}

new_alfak2_graph <- function(x) {
  class(x) <- "alfak2_graph"
  x
}

new_alfak2_local_fit <- function(x) {
  class(x) <- "alfak2_local_fit"
  x
}

new_alfak2_global_fit <- function(x) {
  class(x) <- c("alfak2_global_fit", setdiff(class(x), "alfak2_global_fit"))
  x
}

new_alfak2_fit <- function(x) {
  class(x) <- "alfak2_fit"
  x
}

#' @export
print.alfak2_data <- function(x, ...) {
  cat("<alfak2_data>\n")
  cat("  karyotypes:", nrow(x$counts), "\n")
  cat("  timepoints:", paste(colnames(x$counts), collapse = ", "), "\n")
  cat("  dt:", x$dt, " beta:", x$beta, "\n")
  invisible(x)
}

#' @export
print.alfak2_graph <- function(x, ...) {
  cat("<alfak2_graph>\n")
  cat("  nodes:", length(x$labels), "\n")
  cat("  chromosomes:", x$n_chr, "\n")
  cat("  shell_depth:", x$shell_depth, "\n")
  invisible(x)
}

#' @export
print.alfak2_fit <- function(x, ...) {
  cat("<alfak2_fit>\n")
  cat("  local nodes:", nrow(x$local$summary), "\n")
  cat("  global nodes:", nrow(x$global$summary), "\n")
  cat("  convergence:", x$local$diagnostics$convergence, "\n")
  if (!is.null(x$xval$R2R)) {
    cat("  xval R2R:", x$xval$R2R, "\n")
  }
  invisible(x)
}
