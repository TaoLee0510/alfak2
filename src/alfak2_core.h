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

inline double chromosome_step_weight(int parent_cn, int child_cn, double beta) {
  int d = std::abs(parent_cn - child_cn);
  if (d == 0) return 1.0;
  if (d != 1) return 0.0;
  double burden = std::max(1, parent_cn);
  return std::max(0.0, beta) * burden;
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
                                        double beta) {
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
        double w = chromosome_step_weight(nodes[i][c], child[c], beta);
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
                                          std::vector<double>& w) {
  std::vector<double> row_sum(n, 0.0);
  for (const auto& e : raw_edges) row_sum[e.from] += e.weight;
  from.clear(); to.clear(); w.clear();
  from.reserve(raw_edges.size() + n);
  to.reserve(raw_edges.size() + n);
  w.reserve(raw_edges.size() + n);
  for (const auto& e : raw_edges) {
    double denom = 1.0 + row_sum[e.from];
    from.push_back(e.from);
    to.push_back(e.to);
    w.push_back(e.weight / denom);
  }
  for (int i = 0; i < n; ++i) {
    from.push_back(i);
    to.push_back(i);
    w.push_back(1.0 / (1.0 + row_sum[i]));
  }
}

inline double mean_copy_number_row(const Rcpp::IntegerMatrix& mat, int i) {
  double out = 0.0;
  for (int j = 0; j < mat.ncol(); ++j) out += mat(i, j);
  return out / std::max(1, mat.ncol());
}

} // namespace alfak2

#endif
