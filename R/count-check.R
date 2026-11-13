#' Check single-cell karyotype counts before fitting
#'
#' Reports, for each sample (column) of a count matrix, how much of it a
#' sparse joint fit can build on. Nothing is rejected or changed.
#'
#' The sparse joint model anchors on karyotypes seen in several cells and
#' fits each sample's cells that lie within one step of them. A sample in
#' which almost every cell has its own karyotype, far from every other cell
#' of the same sample, gives it little to anchor on. That pattern can be real
#' diversity, but it is also what copy-number calls that vary from cell to
#' cell produce. When a timepoint has no cell next to a major karyotype,
#' [fit_sparse_joint()] stops with an error of class
#' `alfak2_domain_coverage_error`; a fit that succeeds reports the share of
#' each sample inside its domain in `diagnostics$domain_coverage`.
#'
#' For orientation: in the single-cell samples of Salehi et al. (2021) used by
#' the ALFA-K benchmark, `nearest_distance` was 0 or 1 for most samples and
#' about 8 for two samples in which every cell had a different karyotype.
#'
#' @param counts Count matrix with karyotype row names and one column per
#'   sample, or an `alfak2_data` object.
#' @param major_count Number of cells from which a karyotype counts as
#'   repeated, as `major_count` of [fit_sparse_joint()].
#'
#' @return A data frame with one row per sample: `sample`, `cells`,
#'   `karyotypes`, `repeated_fraction` (share of cells on karyotypes with at
#'   least `major_count` cells in that sample), `nearest_distance` (median over
#'   cells of the distance, summed over chromosomes, to the nearest other cell
#'   of the same sample; 0 for a cell whose karyotype is seen more than once),
#'   `nullisomy_fraction` (share of cells with a chromosome at copy number 0)
#'   and `shared` (karyotypes the sample shares with any other sample; `NA`
#'   for a single sample).
#' @export
check_sparse_joint_counts <- function(counts, major_count = 3) {
  if (inherits(counts, "alfak2_data")) counts <- counts$counts
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || !is.numeric(counts) || ncol(counts) < 1L) {
    stop("`counts` must be a numeric matrix with one column per sample.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
    stop("`counts` must have non-empty karyotype row names.", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0) || any(abs(counts - round(counts)) > 1e-8)) {
    stop("`counts` must contain non-negative integer counts.", call. = FALSE)
  }
  validate_scalar(major_count, "major_count", lower = 1)
  samples <- if (is.null(colnames(counts))) paste0("sample", seq_len(ncol(counts))) else colnames(counts)
  karyotypes <- parse_karyotypes(rownames(counts))
  rows <- lapply(seq_len(ncol(counts)), function(s) {
    cells <- counts[, s]
    seen <- cells > 0
    n <- sum(cells)
    others <- if (ncol(counts) > 1L) rowSums(counts[, -s, drop = FALSE] > 0) > 0 else NULL
    data.frame(
      sample = samples[[s]], cells = n, karyotypes = sum(seen),
      repeated_fraction = if (n > 0) sum(cells[cells >= major_count]) / n else NA_real_,
      nearest_distance = if (n > 1) count_check_nearest_distance(karyotypes[seen, , drop = FALSE], cells[seen]) else NA_real_,
      nullisomy_fraction = if (n > 0) sum(cells[rowSums(karyotypes == 0L) > 0]) / n else NA_real_,
      shared = if (is.null(others)) NA_integer_ else sum(seen & others),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Median over cells of the distance to the nearest other cell of the same
# sample: 0 for a karyotype with more than one cell, otherwise the distance
# to the nearest other karyotype seen. Distances are computed in blocks of
# rows to bound memory.
count_check_nearest_distance <- function(karyotypes, cells) {
  n <- nrow(karyotypes)
  nearest <- rep(Inf, n)
  for (block in split(seq_len(n), ceiling(seq_len(n) / 200))) {
    distance <- matrix(0, length(block), n)
    for (chr in seq_len(ncol(karyotypes))) {
      distance <- distance + abs(outer(karyotypes[block, chr], karyotypes[, chr], "-"))
    }
    distance[cbind(seq_along(block), block)] <- Inf
    nearest[block] <- apply(distance, 1L, min)
  }
  nearest[cells > 1] <- 0
  stats::median(rep(nearest, cells))
}
