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
                                          std::string transition_kernel = "exact",
                                          bool division_process = false,
                                          double probability_cutoff = 1e-14) {
  transition_kernel = alfak2::normalize_transition_kernel(transition_kernel);
  if (beta < 0.0 || beta > 1.0 || probability_cutoff < 0.0 || probability_cutoff >= 1.0)
    Rcpp::stop("Invalid division transition probability controls.");
  Rcpp::CharacterVector labels = alfak2::matrix_to_labels(karyotypes);
  int n = labels.size();
  std::vector< std::vector<int> > nodes(n);
  std::unordered_map<std::string, int> id;
  for (int i = 0; i < n; ++i) {
    nodes[i].resize(karyotypes.ncol());
    for (int j = 0; j < karyotypes.ncol(); ++j) nodes[i][j] = karyotypes(i, j);
    id[alfak2::make_key(nodes[i])] = i;
  }
  // A division kernel is a single-daughter probability, not a graph random
  // walk. Keep all retained transitions, including simultaneous changes, and
  // never renormalize away death or movement outside the computational graph.
  if (division_process) {
    if (transition_kernel != "exact")
      Rcpp::stop("Division dynamics require the exact single-daughter kernel.");
    std::vector<int> from, to;
    std::vector<double> weight;
    Rcpp::NumericVector retained(n);
    int max_cn = 0;
    for (const auto& x : nodes) for (int v : x) max_cn = std::max(max_cn, v);
    std::vector<std::vector<double>> lookup(max_cn + 1, std::vector<double>(max_cn + 1));
    for (int i = 0; i <= max_cn; ++i) for (int j = 0; j <= max_cn; ++j)
      lookup[i][j] = alfak2::chromosome_transition_probability(i, j, beta);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double q = 1.0;
        for (int c = 0; c < karyotypes.ncol(); ++c) {
          q *= lookup[nodes[i][c]][nodes[j][c]];
          if (q <= probability_cutoff) break;
        }
        if (q > probability_cutoff) {
          from.push_back(i); to.push_back(j); weight.push_back(q);
          retained[i] += q;
        }
      }
    }
    return Rcpp::List::create(
      Rcpp::Named("from0") = Rcpp::wrap(from),
      Rcpp::Named("to0") = Rcpp::wrap(to),
      Rcpp::Named("weight") = Rcpp::wrap(weight),
      Rcpp::Named("retained") = retained,
      Rcpp::Named("transition_kernel") = "exact",
      Rcpp::Named("probability_cutoff") = probability_cutoff
    );
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
