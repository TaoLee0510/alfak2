#ifndef ALFAK2_CORE_H
#define ALFAK2_CORE_H

#include <RcppEigen.h>
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace alfak2 {

inline std::string make_key(const std::vector<int>& x) {
  std::string out;
  out.reserve(x.size() * 3);
  for (size_t i = 0; i < x.size(); ++i) {
    if (i) out.push_back('.');
    out += std::to_string(x[i]);
  }
  return out;
}

inline std::vector<int> parse_label(const std::string& s) {
  std::vector<int> out;
  std::string token;
  std::stringstream ss(s);
  while (std::getline(ss, token, '.')) {
    if (token.empty()) Rcpp::stop("Malformed karyotype label: empty copy-number field.");
    char* end = nullptr;
    long v = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0') {
      Rcpp::stop("Malformed karyotype label '%s': copy numbers must be integers.", s);
    }
    if (v < 0 || v > std::numeric_limits<int>::max()) {
      Rcpp::stop("Malformed karyotype label '%s': copy number outside integer range.", s);
    }
    out.push_back(static_cast<int>(v));
  }
  if (out.empty()) Rcpp::stop("Malformed karyotype label: empty string.");
  return out;
}

inline Rcpp::IntegerMatrix labels_to_matrix(const Rcpp::CharacterVector& labels) {
  int n = labels.size();
  if (n == 0) Rcpp::stop("At least one karyotype label is required.");
  std::vector<int> first = parse_label(Rcpp::as<std::string>(labels[0]));
  int p = first.size();
  Rcpp::IntegerMatrix mat(n, p);
  for (int j = 0; j < p; ++j) mat(0, j) = first[j];
  std::unordered_set<std::string> seen;
  seen.insert(make_key(first));
  for (int i = 1; i < n; ++i) {
    std::vector<int> x = parse_label(Rcpp::as<std::string>(labels[i]));
    if (static_cast<int>(x.size()) != p) {
      Rcpp::stop("All karyotypes must have the same chromosome dimension.");
    }
    std::string key = make_key(x);
    if (!seen.insert(key).second) Rcpp::stop("Duplicate karyotype label: %s", key);
    for (int j = 0; j < p; ++j) mat(i, j) = x[j];
  }
  return mat;
}

inline Rcpp::CharacterVector matrix_to_labels(const Rcpp::IntegerMatrix& mat) {
  int n = mat.nrow(), p = mat.ncol();
  Rcpp::CharacterVector labels(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int> x(p);
    for (int j = 0; j < p; ++j) x[j] = mat(i, j);
    labels[i] = make_key(x);
  }
  return labels;
}

inline double choose_double(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  if (k == 0 || k == n) return 1.0;
  k = std::min(k, n - k);
  double out = 1.0;
  for (int i = 1; i <= k; ++i) {
    out *= static_cast<double>(n - k + i) / i;
  }
  return out;
}

inline std::string normalize_transition_kernel(const std::string& transition_kernel) {
  if (transition_kernel == "exact" || transition_kernel == "linear") return transition_kernel;
  Rcpp::stop("`transition_kernel` must be \"exact\" or \"linear\".");
}

inline double chromosome_transition_probability(int parent_cn, int child_cn, double beta) {
  if (parent_cn < 0 || child_cn < 0) return 0.0;
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

inline double linear_step_weight(int parent_cn, int child_cn, double beta) {
  int d = std::abs(parent_cn - child_cn);
  if (d == 0) return 1.0;
  if (d != 1) return 0.0;
  double burden = std::max(1, parent_cn);
  return std::max(0.0, beta) * burden;
}

inline double state_transition_weight(const std::vector<int>& parent,
                                      const std::vector<int>& child,
                                      double beta,
                                      const std::string& transition_kernel) {
  if (transition_kernel == "linear") {
    int changed = -1;
    for (size_t c = 0; c < parent.size(); ++c) {
      if (parent[c] != child[c]) {
        if (changed >= 0) return 0.0;
        changed = static_cast<int>(c);
      }
    }
    if (changed < 0) return 1.0;
    return linear_step_weight(parent[changed], child[changed], beta);
  }
  double out = 1.0;
  for (size_t c = 0; c < parent.size(); ++c) {
    out *= chromosome_transition_probability(parent[c], child[c], beta);
  }
  return out;
}

inline std::vector<double> state_self_weights(const std::vector< std::vector<int> >& nodes,
                                              double beta,
                                              const std::string& transition_kernel) {
  std::vector<double> out(nodes.size(), 1.0);
  if (transition_kernel == "linear") return out;
  for (size_t i = 0; i < nodes.size(); ++i) {
    out[i] = state_transition_weight(nodes[i], nodes[i], beta, transition_kernel);
  }
  return out;
}

inline int ploidy_band(const std::vector<int>& x) {
  double mean = 0.0;
  for (int v : x) mean += v;
  mean /= std::max<size_t>(1, x.size());
  if (mean < 1.75) return 0;
  if (mean < 2.75) return 1;
  return 2;
}

struct Edge {
  int from;
  int to;
  int chr;
  int direction;
  double weight;
};

inline std::vector<Edge> one_step_edges(const std::vector< std::vector<int> >& nodes,
                                        const std::unordered_map<std::string, int>& id,
                                        double beta,
                                        std::string transition_kernel = "exact") {
  transition_kernel = normalize_transition_kernel(transition_kernel);
  std::vector<Edge> edges;
  int n = nodes.size();
  if (n == 0) return edges;
  int p = nodes[0].size();
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < p; ++c) {
      for (int dir : {-1, 1}) {
        std::vector<int> child = nodes[i];
        child[c] += dir;
        auto it = id.find(make_key(child));
        if (it == id.end()) continue;
        int j = it->second;
        double w = state_transition_weight(nodes[i], child, beta, transition_kernel);
        if (w > 0) edges.push_back({i, j, c, dir, w});
      }
    }
  }
  return edges;
}

inline void add_row_stochastic_transition(const std::vector<Edge>& raw_edges,
                                          int n,
                                          std::vector<int>& from,
                                          std::vector<int>& to,
                                          std::vector<double>& w,
                                          const std::vector<double>& self_weight = std::vector<double>()) {
  std::vector<double> row_sum(n, 0.0);
  for (const auto& e : raw_edges) row_sum[e.from] += e.weight;
  from.clear(); to.clear(); w.clear();
  from.reserve(raw_edges.size() + n);
  to.reserve(raw_edges.size() + n);
  w.reserve(raw_edges.size() + n);
  for (const auto& e : raw_edges) {
    double self = self_weight.size() == static_cast<size_t>(n) ? self_weight[e.from] : 1.0;
    double denom = self + row_sum[e.from];
    if (denom <= 0.0 || !std::isfinite(denom)) denom = 1.0;
    from.push_back(e.from);
    to.push_back(e.to);
    w.push_back(e.weight / denom);
  }
  for (int i = 0; i < n; ++i) {
    double self = self_weight.size() == static_cast<size_t>(n) ? self_weight[i] : 1.0;
    double denom = self + row_sum[i];
    if (denom <= 0.0 || !std::isfinite(denom)) {
      self = 1.0;
      denom = 1.0;
    }
    from.push_back(i);
    to.push_back(i);
    w.push_back(self / denom);
  }
}

inline double mean_copy_number_row(const Rcpp::IntegerMatrix& mat, int i) {
  double out = 0.0;
  for (int j = 0; j < mat.ncol(); ++j) out += mat(i, j);
  return out / std::max(1, mat.ncol());
}

} // namespace alfak2

#endif
