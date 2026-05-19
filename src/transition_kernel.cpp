// [[Rcpp::depends(RcppEigen)]]
#include "alfak2_core.h"

using namespace Rcpp;

// Chromosome-level transition probability used for diagnostics and simulation.
// [[Rcpp::export]]
double alfak2_pij_cpp(int parent_cn, int child_cn, double beta) {
  if (beta < 0.0 || beta > 1.0) Rcpp::stop("`beta` must be in [0, 1].");
  return alfak2::chromosome_transition_probability(parent_cn, child_cn, beta);
}

// [[Rcpp::export]]
Rcpp::List alfak2_transition_operator_cpp(Rcpp::IntegerMatrix karyotypes,
                                          double beta = 0.00005,
                                          std::string transition_kernel = "exact") {
  transition_kernel = alfak2::normalize_transition_kernel(transition_kernel);
  Rcpp::CharacterVector labels = alfak2::matrix_to_labels(karyotypes);
  int n = labels.size();
  std::vector< std::vector<int> > nodes(n);
  std::unordered_map<std::string, int> id;
  for (int i = 0; i < n; ++i) {
    nodes[i].resize(karyotypes.ncol());
    for (int j = 0; j < karyotypes.ncol(); ++j) nodes[i][j] = karyotypes(i, j);
    id[alfak2::make_key(nodes[i])] = i;
  }
  std::vector<alfak2::Edge> raw_edges = alfak2::one_step_edges(nodes, id, beta, transition_kernel);
  std::vector<double> self_weight = alfak2::state_self_weights(nodes, beta, transition_kernel);
  std::vector<int> from, to;
  std::vector<double> weight;
  alfak2::add_row_stochastic_transition(raw_edges, n, from, to, weight, self_weight);
  return Rcpp::List::create(
    Rcpp::Named("from0") = Rcpp::wrap(from),
    Rcpp::Named("to0") = Rcpp::wrap(to),
    Rcpp::Named("weight") = Rcpp::wrap(weight),
    Rcpp::Named("transition_kernel") = transition_kernel
  );
}
