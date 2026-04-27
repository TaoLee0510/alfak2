// [[Rcpp::depends(RcppEigen)]]
#include "alfak2_core.h"

using Eigen::SparseMatrix;
using Eigen::Triplet;
using Eigen::VectorXd;

struct GraphSolveResult {
  VectorXd mean;
  VectorXd variance;
  bool ok;
  std::string status;
};

static void add_penalty(std::vector< Triplet<double> >& triplets,
                        const std::vector< std::pair<int, double> >& coeff,
                        double scale) {
  if (scale == 0.0) return;
  for (size_t a = 0; a < coeff.size(); ++a) {
    for (size_t b = 0; b < coeff.size(); ++b) {
      triplets.emplace_back(coeff[a].first, coeff[b].first,
                            scale * coeff[a].second * coeff[b].second);
    }
  }
}

static SparseMatrix<double> assemble_precision(const Rcpp::IntegerMatrix& karyotypes,
                                               const Rcpp::IntegerVector& edge_from,
                                               const Rcpp::IntegerVector& edge_to,
                                               const Rcpp::NumericVector& edge_weight,
                                               const std::vector<int>& anchor_index0,
                                               const Rcpp::NumericVector& anchor_var,
                                               double lambda_l,
                                               double lambda_e,
                                               double sigma_obs,
                                               double eps) {
  int n = karyotypes.nrow();
  int p = karyotypes.ncol();
  std::vector< Triplet<double> > triplets;
  triplets.reserve(edge_from.size() * 4 + n * (p + 1) * 9);

  for (int e = 0; e < edge_from.size(); ++e) {
    int i = edge_from[e] - 1;
    int j = edge_to[e] - 1;
    if (i == j) continue;
    double w = lambda_l * std::max(0.0, edge_weight[e]);
    if (w == 0.0) continue;
    triplets.emplace_back(i, i, w);
    triplets.emplace_back(j, j, w);
    triplets.emplace_back(i, j, -w);
    triplets.emplace_back(j, i, -w);
  }

  std::unordered_map<std::string, int> id;
  id.reserve(n * 2);
  for (int i = 0; i < n; ++i) {
    std::vector<int> x(p);
    for (int c = 0; c < p; ++c) x[c] = karyotypes(i, c);
    id[alfak2::make_key(x)] = i;
  }

  for (int i = 0; i < n; ++i) {
    std::vector<int> x(p);
    for (int c = 0; c < p; ++c) x[c] = karyotypes(i, c);
    for (int c = 0; c < p; ++c) {
      std::vector<int> lo = x, hi = x;
      lo[c] -= 1;
      hi[c] += 1;
      auto it_lo = id.find(alfak2::make_key(lo));
      auto it_hi = id.find(alfak2::make_key(hi));
      if (it_lo != id.end() && it_hi != id.end()) {
        add_penalty(triplets,
                    {{it_lo->second, 1.0}, {i, -2.0}, {it_hi->second, 1.0}},
                    lambda_e);
      }
    }
    for (int c1 = 0; c1 < p; ++c1) {
      for (int c2 = c1 + 1; c2 < p; ++c2) {
        std::vector<int> xc = x, xd = x, xcd = x;
        xc[c1] += 1;
        xd[c2] += 1;
        xcd[c1] += 1;
        xcd[c2] += 1;
        auto it_c = id.find(alfak2::make_key(xc));
        auto it_d = id.find(alfak2::make_key(xd));
        auto it_cd = id.find(alfak2::make_key(xcd));
        if (it_c != id.end() && it_d != id.end() && it_cd != id.end()) {
          add_penalty(triplets,
                      {{i, 1.0}, {it_c->second, -1.0},
                       {it_d->second, -1.0}, {it_cd->second, 1.0}},
                      0.5 * lambda_e);
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) triplets.emplace_back(i, i, eps);
  for (size_t a = 0; a < anchor_index0.size(); ++a) {
    double v = std::max(1e-10, anchor_var[a] + sigma_obs * sigma_obs);
    triplets.emplace_back(anchor_index0[a], anchor_index0[a], 1.0 / v);
  }

  SparseMatrix<double> q(n, n);
  q.setFromTriplets(triplets.begin(), triplets.end());
  q.makeCompressed();
  return q;
}

static GraphSolveResult solve_graph(const Rcpp::IntegerMatrix& karyotypes,
                                    const Rcpp::IntegerVector& edge_from,
                                    const Rcpp::IntegerVector& edge_to,
                                    const Rcpp::NumericVector& edge_weight,
                                    const std::vector<int>& anchor_index0,
                                    const Rcpp::NumericVector& anchor_mean,
                                    const Rcpp::NumericVector& anchor_var,
                                    double lambda_l,
                                    double lambda_e,
                                    double sigma_obs,
                                    double eps,
                                    bool compute_variance) {
  int n = karyotypes.nrow();
  SparseMatrix<double> q = assemble_precision(karyotypes, edge_from, edge_to, edge_weight,
                                              anchor_index0, anchor_var,
                                              lambda_l, lambda_e, sigma_obs, eps);
  VectorXd rhs = VectorXd::Zero(n);
  for (size_t a = 0; a < anchor_index0.size(); ++a) {
    double v = std::max(1e-10, anchor_var[a] + sigma_obs * sigma_obs);
    rhs(anchor_index0[a]) += anchor_mean[a] / v;
  }

  Eigen::SimplicialLDLT< SparseMatrix<double> > solver;
  solver.compute(q);
  if (solver.info() != Eigen::Success) {
    return {VectorXd::Constant(n, NA_REAL), VectorXd::Constant(n, NA_REAL),
            false, "SimplicialLDLT factorization failed"};
  }
  VectorXd mean = solver.solve(rhs);
  if (solver.info() != Eigen::Success) {
    return {VectorXd::Constant(n, NA_REAL), VectorXd::Constant(n, NA_REAL),
            false, "Sparse posterior mean solve failed"};
  }

  VectorXd var = VectorXd::Constant(n, NA_REAL);
  if (compute_variance) {
    var.setZero();
    for (int i = 0; i < n; ++i) {
      VectorXd e = VectorXd::Zero(n);
      e(i) = 1.0;
      VectorXd z = solver.solve(e);
      var(i) = std::max(0.0, z(i));
    }
  }
  return {mean, var, true, "ok"};
}

static double cv_score_grid(const Rcpp::IntegerMatrix& karyotypes,
                            const Rcpp::IntegerVector& edge_from,
                            const Rcpp::IntegerVector& edge_to,
                            const Rcpp::NumericVector& edge_weight,
                            const std::vector<int>& anchor_index0,
                            const Rcpp::NumericVector& anchor_mean,
                            const Rcpp::NumericVector& anchor_var,
                            double lambda_l,
                            double lambda_e,
                            double sigma_obs,
                            double eps) {
  int m = anchor_index0.size();
  if (m < 3) return 0.0;
  double rss = 0.0;
  for (int hold = 0; hold < m; ++hold) {
    std::vector<int> idx;
    Rcpp::NumericVector mean(m - 1), var(m - 1);
    idx.reserve(m - 1);
    int k = 0;
    for (int a = 0; a < m; ++a) {
      if (a == hold) continue;
      idx.push_back(anchor_index0[a]);
      mean[k] = anchor_mean[a];
      var[k] = anchor_var[a];
      ++k;
    }
    GraphSolveResult res = solve_graph(karyotypes, edge_from, edge_to, edge_weight,
                                       idx, mean, var, lambda_l, lambda_e,
                                       sigma_obs, eps, false);
    if (!res.ok) return std::numeric_limits<double>::infinity();
    double pred = res.mean(anchor_index0[hold]);
    double err = pred - anchor_mean[hold];
    rss += err * err / std::max(1e-10, anchor_var[hold] + sigma_obs * sigma_obs);
  }
  return rss / m;
}

// [[Rcpp::export]]
Rcpp::List alfak2_graph_posterior_cpp(Rcpp::IntegerMatrix karyotypes,
                                      Rcpp::IntegerVector edge_from,
                                      Rcpp::IntegerVector edge_to,
                                      Rcpp::NumericVector edge_weight,
                                      Rcpp::IntegerVector anchor_index,
                                      Rcpp::NumericVector anchor_mean,
                                      Rcpp::NumericVector anchor_var,
                                      Rcpp::NumericVector lambda_l_grid,
                                      Rcpp::NumericVector lambda_e_grid,
                                      Rcpp::NumericVector sigma_obs_grid,
                                      double eps = 1e-5) {
  if (anchor_index.size() == 0) Rcpp::stop("At least one local anchor is required.");
  if (anchor_mean.size() != anchor_index.size() || anchor_var.size() != anchor_index.size()) {
    Rcpp::stop("Anchor index, mean, and variance vectors must have equal length.");
  }
  std::vector<int> anchor0(anchor_index.size());
  for (int i = 0; i < anchor_index.size(); ++i) {
    anchor0[i] = anchor_index[i] - 1;
    if (anchor0[i] < 0 || anchor0[i] >= karyotypes.nrow()) {
      Rcpp::stop("Anchor index outside graph node range.");
    }
  }

  if (lambda_l_grid.size() == 0) lambda_l_grid = Rcpp::NumericVector::create(1.0);
  if (lambda_e_grid.size() == 0) lambda_e_grid = Rcpp::NumericVector::create(0.25);
  if (sigma_obs_grid.size() == 0) sigma_obs_grid = Rcpp::NumericVector::create(0.05);

  double best_score = std::numeric_limits<double>::infinity();
  double best_l = lambda_l_grid[0], best_e = lambda_e_grid[0], best_s = sigma_obs_grid[0];
  Rcpp::DataFrame grid = Rcpp::DataFrame::create();
  std::vector<double> gl, ge, gs, score;
  for (double l : lambda_l_grid) {
    for (double e : lambda_e_grid) {
      for (double s : sigma_obs_grid) {
        double sc = cv_score_grid(karyotypes, edge_from, edge_to, edge_weight,
                                  anchor0, anchor_mean, anchor_var,
                                  l, e, s, eps);
        gl.push_back(l); ge.push_back(e); gs.push_back(s); score.push_back(sc);
        if (sc < best_score) {
          best_score = sc;
          best_l = l; best_e = e; best_s = s;
        }
      }
    }
  }

  GraphSolveResult final = solve_graph(karyotypes, edge_from, edge_to, edge_weight,
                                       anchor0, anchor_mean, anchor_var,
                                       best_l, best_e, best_s, eps, true);
  if (!final.ok) Rcpp::stop(final.status);

  Rcpp::NumericVector mean(final.mean.size()), sd(final.variance.size());
  for (int i = 0; i < final.mean.size(); ++i) {
    mean[i] = final.mean(i);
    sd[i] = std::sqrt(std::max(0.0, final.variance(i)));
  }
  return Rcpp::List::create(
    Rcpp::Named("mean") = mean,
    Rcpp::Named("sd") = sd,
    Rcpp::Named("lambda_l") = best_l,
    Rcpp::Named("lambda_e") = best_e,
    Rcpp::Named("sigma_obs") = best_s,
    Rcpp::Named("cv_score") = best_score,
    Rcpp::Named("grid") = Rcpp::DataFrame::create(
      Rcpp::Named("lambda_l") = gl,
      Rcpp::Named("lambda_e") = ge,
      Rcpp::Named("sigma_obs") = gs,
      Rcpp::Named("score") = score
    ),
    Rcpp::Named("factorization_status") = final.status
  );
}
