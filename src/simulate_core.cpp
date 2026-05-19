// [[Rcpp::depends(RcppEigen)]]
#include "alfak2_core.h"

using namespace Rcpp;

static Rcpp::IntegerVector multinomial_draw(const Rcpp::NumericVector& prob,
                                            int size,
                                            std::mt19937& rng,
                                            double overdispersion) {
  int n = prob.size();
  std::vector<double> p(n);
  if (overdispersion > 0.0) {
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      double alpha = std::max(1e-10, overdispersion * prob[i]);
      std::gamma_distribution<double> gamma(alpha, 1.0);
      p[i] = gamma(rng);
      sum += p[i];
    }
    for (int i = 0; i < n; ++i) p[i] /= sum;
  } else {
    for (int i = 0; i < n; ++i) p[i] = prob[i];
  }
  std::discrete_distribution<int> draw(p.begin(), p.end());
  Rcpp::IntegerVector out(n);
  for (int k = 0; k < size; ++k) out[draw(rng)]++;
  return out;
}

static Rcpp::List simulate_counts_impl(Rcpp::IntegerMatrix karyotypes,
                                       Rcpp::CharacterVector labels,
                                       Rcpp::NumericVector fitness,
                                       double beta,
                                       double dt,
                                       int n0,
                                       int n1,
                                       int detection_threshold,
                                       double dropout_prob,
                                       int seed,
                                       double init_concentration,
                                       double overdispersion) {
  int n = karyotypes.nrow();
  int p = karyotypes.ncol();
  if (fitness.size() != n || labels.size() != n) Rcpp::stop("Landscape vectors have inconsistent lengths.");
  std::mt19937 rng(seed);

  Rcpp::NumericVector pi0(n);
  double z0 = 0.0;
  for (int i = 0; i < n; ++i) {
    double d2 = 0.0;
    for (int c = 0; c < p; ++c) {
      double d = karyotypes(i, c) - 2.0;
      d2 += d * d;
    }
    pi0[i] = std::exp(-init_concentration * d2);
    z0 += pi0[i];
  }
  for (int i = 0; i < n; ++i) pi0[i] /= z0;

  std::vector< std::vector<int> > nodes(n, std::vector<int>(p));
  std::unordered_map<std::string, int> id;
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < p; ++c) nodes[i][c] = karyotypes(i, c);
    id[alfak2::make_key(nodes[i])] = i;
  }
  std::vector<alfak2::Edge> raw_edges = alfak2::one_step_edges(nodes, id, beta);
  std::vector<int> from, to;
  std::vector<double> weight;
  alfak2::add_row_stochastic_transition(raw_edges, n, from, to, weight);

  std::vector<double> growth(n);
  double maxg = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) maxg = std::max(maxg, dt * fitness[i]);
  for (int i = 0; i < n; ++i) growth[i] = pi0[i] * std::exp(dt * fitness[i] - maxg);
  Rcpp::NumericVector pi1(n);
  for (size_t e = 0; e < weight.size(); ++e) pi1[to[e]] += weight[e] * growth[from[e]];
  double z1 = std::accumulate(pi1.begin(), pi1.end(), 0.0);
  for (int i = 0; i < n; ++i) pi1[i] /= z1;

  Rcpp::IntegerVector y0 = multinomial_draw(pi0, n0, rng, overdispersion);
  Rcpp::IntegerVector y1 = multinomial_draw(pi1, n1, rng, overdispersion);

  std::uniform_real_distribution<double> unif(0.0, 1.0);
  std::vector<int> keep;
  keep.reserve(n);
  int best = 0, best_count = -1;
  for (int i = 0; i < n; ++i) {
    int total = y0[i] + y1[i];
    if (total > best_count) {
      best_count = total;
      best = i;
    }
    if (total >= detection_threshold && unif(rng) >= dropout_prob) keep.push_back(i);
  }
  if (keep.empty()) keep.push_back(best);

  Rcpp::IntegerMatrix counts(keep.size(), 2);
  Rcpp::CharacterVector obs_labels(keep.size());
  Rcpp::NumericVector truth_keep(keep.size());
  for (size_t k = 0; k < keep.size(); ++k) {
    int i = keep[k];
    counts(k, 0) = y0[i];
    counts(k, 1) = y1[i];
    obs_labels[k] = labels[i];
    truth_keep[k] = fitness[i];
  }
  Rcpp::CharacterVector cn = Rcpp::CharacterVector::create("t0", "t1");
  Rcpp::colnames(counts) = cn;
  Rcpp::rownames(counts) = obs_labels;

  double entropy = 0.0;
  for (int i = 0; i < n; ++i) {
    double p0 = (y0[i] + y1[i]) / static_cast<double>(n0 + n1);
    if (p0 > 0.0) entropy -= p0 * std::log(p0);
  }
  return Rcpp::List::create(
    Rcpp::Named("counts") = counts,
    Rcpp::Named("observed_labels") = obs_labels,
    Rcpp::Named("truth_observed") = truth_keep,
    Rcpp::Named("sparsity") = Rcpp::List::create(
      Rcpp::Named("observed_nodes") = static_cast<int>(keep.size()),
      Rcpp::Named("fraction_observed") = keep.size() / static_cast<double>(n),
      Rcpp::Named("effective_support_entropy") = entropy,
      Rcpp::Named("sampling_depth") = n0 + n1,
      Rcpp::Named("detection_threshold") = detection_threshold,
      Rcpp::Named("dropout_prob") = dropout_prob
    )
  );
}

// [[Rcpp::export]]
Rcpp::List alfak2_simulate_counts_cpp(Rcpp::IntegerMatrix karyotypes,
                                      Rcpp::CharacterVector labels,
                                      Rcpp::NumericVector fitness,
                                      double beta = 0.00005,
                                      double dt = 1.0,
                                      int n0 = 500,
                                      int n1 = 500,
                                      int detection_threshold = 1,
                                      double dropout_prob = 0.0,
                                      int seed = 1,
                                      double init_concentration = 1.5,
                                      double overdispersion = 0.0) {
  return simulate_counts_impl(karyotypes, labels, fitness, beta, dt, n0, n1,
                              detection_threshold, dropout_prob, seed,
                              init_concentration, overdispersion);
}
