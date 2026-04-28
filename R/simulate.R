#' Simulate sparse two-timepoint counts from a truth landscape
#'
#' @param landscape Object returned by `simulate_l1_gp_landscape()`, or another
#'   landscape-like list containing `karyotypes`, `labels`, and `fitness`.
#' @param beta,dt Dynamics parameters.
#' @param n0,n1 Sampling depths.
#' @param detection_threshold Minimum total observed count retained.
#' @param dropout_prob Random node dropout probability after sampling.
#' @param seed Reproducibility seed.
#' @param init_concentration Concentration around diploid-like center.
#' @param overdispersion Dirichlet concentration; `0` means multinomial.
#'
#' @return A list containing a two-column count matrix and sparsity metadata.
#' @export
simulate_sparse_counts <- function(landscape,
                                   beta = 0.01,
                                   dt = 1,
                                   n0 = 500,
                                   n1 = 500,
                                   detection_threshold = 1,
                                   dropout_prob = 0,
                                   seed = 1,
                                   init_concentration = 1.5,
                                   overdispersion = 0) {
  alfak2_simulate_counts_cpp(
    karyotypes = landscape$karyotypes,
    labels = landscape$labels,
    fitness = landscape$fitness,
    beta = beta,
    dt = dt,
    n0 = as.integer(n0),
    n1 = as.integer(n1),
    detection_threshold = as.integer(detection_threshold),
    dropout_prob = dropout_prob,
    seed = as.integer(seed),
    init_concentration = init_concentration,
    overdispersion = overdispersion
  )
}
