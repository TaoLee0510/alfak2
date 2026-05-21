#' Summarize alfak2 fits
#'
#' @param object An `alfak2_fit`, `alfak2_local_fit`, or graph posterior list.
#' @param layer For full fits, `"global"` or `"local"`.
#' @param ... Reserved.
#'
#' @return Data frame of posterior node summaries.
#' @export
summarize_alfak2 <- function(object, layer = c("global", "local"), ...) {
  layer <- match.arg(layer)
  if (inherits(object, "alfak2_fit")) {
    if (layer == "global") return(object$global$summary)
    return(object$local$summary)
  }
  if (inherits(object, "alfak2_local_fit")) return(object$summary)
  if (is.list(object) && !is.null(object$summary)) return(object$summary)
  stop("Unsupported object type for summarize_alfak2().", call. = FALSE)
}

#' @export
summary.alfak2_fit <- function(object, ...) {
  summarize_alfak2(object, ...)
}

#' @export
summary.alfak2_local_fit <- function(object, ...) {
  summarize_alfak2(object, ...)
}

#' @export
summary.alfak2_global_fit <- function(object, ...) {
  summarize_alfak2(object, ...)
}
