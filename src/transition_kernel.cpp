// [[Rcpp::depends(RcppEigen)]]
#include "alfak2_core.h"

using namespace Rcpp;

static double choose_double(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  if (k == 0 || k == n) return 1.0;
  k = std::min(k, n - k);
  double out = 1.0;
  for (int i = 1; i <= k; ++i) out *= static_cast<double>(n - k + i) / i;
  return out;
}

// Chromosome-level transition probability used for diagnostics and simulation.
// [[Rcpp::export]]
double alfak2_pij_cpp(int parent_cn, int child_cn, double beta) {
  if (parent_cn < 0 || child_cn < 0) return 0.0;
  if (beta < 0.0 || beta > 1.0) Rcpp::stop("`beta` must be in [0, 1].");
  int diff = std::abs(parent_cn - child_cn);
  if (diff > parent_cn) return 0.0;
  double total = 0.0;
  for (int m = diff; m <= parent_cn; m += 2) {
    int lost_to_child = (m + parent_cn - child_cn) / 2;
    if (lost_to_child < 0 || lost_to_child > m) continue;
    total += choose_double(parent_cn, m) *
      std::pow(beta, m) *
      std::pow(1.0 - beta, parent_cn - m) *
      std::pow(0.5, m) *
      choose_double(m, lost_to_child);
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::List alfak2_transition_operator_cpp(Rcpp::IntegerMatrix karyotypes,
                                          double beta = 0.00005) {
  Rcpp::CharacterVector labels = alfak2::matrix_to_labels(karyotypes);
  int n = labels.size();
  std::vector< std::vector<int> > nodes(n);
  std::unordered_map<std::string, int> id;
  for (int i = 0; i < n; ++i) {
    nodes[i].resize(karyotypes.ncol());
    for (int j = 0; j < karyotypes.ncol(); ++j) nodes[i][j] = karyotypes(i, j);
    id[alfak2::make_key(nodes[i])] = i;
  }
  std::vector<alfak2::Edge> raw_edges = alfak2::one_step_edges(nodes, id, beta);
  std::vector<int> from, to;
  std::vector<double> weight;
  alfak2::add_row_stochastic_transition(raw_edges, n, from, to, weight);
  return Rcpp::List::create(
    Rcpp::Named("from0") = Rcpp::wrap(from),
    Rcpp::Named("to0") = Rcpp::wrap(to),
    Rcpp::Named("weight") = Rcpp::wrap(weight)
  );
}
