list_alfak2_input_modes <- function() {
  c("full", "minobs_matched", "soft_minobs")
}

prepare_alfak2_input <- function(mode, sim, landscape, config = list()) {
  counts <- sim$counts
  totals <- rowSums(counts)
  if (identical(mode, "full")) {
    return(list(counts = counts, anchor_count_reference = NULL, input_depth = "raw"))
  }
  if (identical(mode, "minobs_matched")) {
    min_total <- config$min_total %||% 10
    keep <- totals >= min_total
    if (sum(keep) < 3L) {
      keep[order(totals, decreasing = TRUE)[seq_len(min(3L, length(totals)))]] <- TRUE
    }
    return(list(counts = counts[keep, , drop = FALSE], anchor_count_reference = NULL, input_depth = "raw"))
  }
  if (identical(mode, "soft_minobs")) {
    return(list(counts = counts, anchor_count_reference = config$anchor_count_reference %||% 10, input_depth = "raw"))
  }
  stop("Unknown alfak2 input mode: ", mode, call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
