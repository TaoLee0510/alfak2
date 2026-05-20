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

struct CvScoreResult {
  double score;
  int evaluated;
  int skipped;
};

static std::vector<int> graph_components(int n,
                                         const Rcpp::IntegerVector& edge_from,
                                         const Rcpp::IntegerVector& edge_to) {
  std::vector< std::vector<int> > adj(n);
  for (int e = 0; e < edge_from.size(); ++e) {
    int i = edge_from[e] - 1;
    int j = edge_to[e] - 1;
    if (i < 0 || i >= n || j < 0 || j >= n || i == j) continue;
    adj[i].push_back(j);
    adj[j].push_back(i);
  }
  std::vector<int> component(n, -1);
  int cid = 0;
  for (int start = 0; start < n; ++start) {
    if (component[start] >= 0) continue;
    std::deque<int> queue;
    queue.push_back(start);
    component[start] = cid;
    while (!queue.empty()) {
      int cur = queue.front();
      queue.pop_front();
      for (int nb : adj[cur]) {
        if (component[nb] >= 0) continue;
        component[nb] = cid;
        queue.push_back(nb);
      }
    }
    ++cid;
  }
  return component;
}

static std::vector<int> anchor_component_counts(const std::vector<int>& anchor_index0,
                                                const std::vector<int>& component) {
  int n_component = 0;
  for (int cid : component) n_component = std::max(n_component, cid + 1);
  std::vector<int> counts(n_component, 0);
  for (int idx : anchor_index0) {
    int cid = component[idx];
    if (cid >= 0) counts[cid] += 1;
  }
  return counts;
}

static int cv_evaluable_anchors(const std::vector<int>& anchor_index0,
                                const std::vector<int>& component,
                                const std::vector<int>& component_anchor_count) {
  int out = 0;
  for (int idx : anchor_index0) {
    int cid = component[idx];
    if (cid >= 0 && component_anchor_count[cid] >= 2) ++out;
  }
  return out;
}

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

static CvScoreResult cv_score_grid(const Rcpp::IntegerMatrix& karyotypes,
                                   const Rcpp::IntegerVector& edge_from,
                                   const Rcpp::IntegerVector& edge_to,
                                   const Rcpp::NumericVector& edge_weight,
                                   const std::vector<int>& anchor_index0,
                                   const Rcpp::NumericVector& anchor_mean,
                                   const Rcpp::NumericVector& anchor_var,
                                   const std::vector<int>& component,
                                   const std::vector<int>& component_anchor_count,
                                   double lambda_l,
                                   double lambda_e,
                                   double sigma_obs,
                                   double eps) {
  int m = anchor_index0.size();
  if (m < 3) return {NA_REAL, 0, m};
  double rss = 0.0;
  int evaluated = 0;
  int skipped = 0;
  for (int hold = 0; hold < m; ++hold) {
    int cid = component[anchor_index0[hold]];
    if (cid < 0 || component_anchor_count[cid] < 2) {
      ++skipped;
      continue;
    }
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
    if (!res.ok) return {std::numeric_limits<double>::infinity(), evaluated, skipped};
    double pred = res.mean(anchor_index0[hold]);
    double err = pred - anchor_mean[hold];
    rss += err * err / std::max(1e-10, anchor_var[hold] + sigma_obs * sigma_obs);
    ++evaluated;
  }
  if (evaluated == 0) return {NA_REAL, evaluated, skipped};
  return {rss / evaluated, evaluated, skipped};
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
                                      double eps = 1e-5,
                                      bool compute_variance = true) {
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

  if (lambda_l_grid.size() == 0 || lambda_e_grid.size() == 0 || sigma_obs_grid.size() == 0) {
    Rcpp::stop("Hyperparameter grids must be non-empty.");
  }

  std::vector<int> component = graph_components(karyotypes.nrow(), edge_from, edge_to);
  std::vector<int> component_anchor_count = anchor_component_counts(anchor0, component);
  int cv_evaluable = cv_evaluable_anchors(anchor0, component, component_anchor_count);
  bool fixed_hyperparameters = lambda_l_grid.size() == 1 &&
    lambda_e_grid.size() == 1 &&
    sigma_obs_grid.size() == 1;
  bool tune_by_cv = !fixed_hyperparameters && anchor0.size() >= 3 && cv_evaluable >= 3;
  int cv_skipped = anchor0.size() - cv_evaluable;
  int cv_report_evaluated = tune_by_cv ? cv_evaluable : 0;
  int cv_report_skipped = tune_by_cv ? cv_skipped : 0;
  std::string cv_status;
  if (fixed_hyperparameters) {
    cv_status = "fixed_hyperparameters";
  } else if (anchor0.size() < 3) {
    cv_status = "insufficient_anchors";
  } else if (cv_evaluable < 3) {
    cv_status = "insufficient_component_anchors";
  } else if (cv_skipped > 0) {
    cv_status = "partial_components";
  } else {
    cv_status = "ok";
  }
  if (!fixed_hyperparameters && !tune_by_cv) {
    Rcpp::stop("Graph hyperparameter tuning cannot run: " + cv_status);
  }
  double best_score = tune_by_cv ? std::numeric_limits<double>::infinity() : NA_REAL;
  double best_l = lambda_l_grid[0];
  double best_e = lambda_e_grid[0];
  double best_s = sigma_obs_grid[0];
  Rcpp::DataFrame grid = Rcpp::DataFrame::create();
  std::vector<double> gl, ge, gs, score;
  for (double l : lambda_l_grid) {
    for (double e : lambda_e_grid) {
      for (double s : sigma_obs_grid) {
        CvScoreResult cv = tune_by_cv ?
          cv_score_grid(karyotypes, edge_from, edge_to, edge_weight,
                        anchor0, anchor_mean, anchor_var,
                        component, component_anchor_count,
                        l, e, s, eps) :
          CvScoreResult{NA_REAL, 0, cv_report_skipped};
        double sc = cv.score;
        gl.push_back(l); ge.push_back(e); gs.push_back(s); score.push_back(sc);
        if (tune_by_cv && R_finite(sc) && sc < best_score) {
          best_score = sc;
          best_l = l; best_e = e; best_s = s;
        }
      }
    }
  }

  if (tune_by_cv && !R_finite(best_score)) {
    Rcpp::stop("Graph hyperparameter tuning failed: all CV scores are non-finite.");
  }

  GraphSolveResult final = solve_graph(karyotypes, edge_from, edge_to, edge_weight,
                                       anchor0, anchor_mean, anchor_var,
                                       best_l, best_e, best_s, eps, compute_variance);
  if (!final.ok) Rcpp::stop(final.status);

  Rcpp::NumericVector mean(final.mean.size()), sd(final.variance.size());
  for (int i = 0; i < final.mean.size(); ++i) {
    mean[i] = final.mean(i);
    sd[i] = compute_variance ? std::sqrt(std::max(0.0, final.variance(i))) : NA_REAL;
  }
  return Rcpp::List::create(
    Rcpp::Named("mean") = mean,
    Rcpp::Named("sd") = sd,
    Rcpp::Named("lambda_l") = best_l,
    Rcpp::Named("lambda_e") = best_e,
    Rcpp::Named("sigma_obs") = best_s,
    Rcpp::Named("cv_score") = best_score,
    Rcpp::Named("cv_status") = cv_status,
    Rcpp::Named("cv_evaluated") = cv_report_evaluated,
    Rcpp::Named("cv_skipped") = cv_report_skipped,
    Rcpp::Named("grid") = Rcpp::DataFrame::create(
      Rcpp::Named("lambda_l") = gl,
      Rcpp::Named("lambda_e") = ge,
      Rcpp::Named("sigma_obs") = gs,
      Rcpp::Named("score") = score
    ),
    Rcpp::Named("factorization_status") = final.status,
    Rcpp::Named("compute_variance") = compute_variance
  );
}

static std::vector< std::vector<int> > build_undirected_adjacency(int n,
                                                                  const Rcpp::IntegerVector& edge_from,
                                                                  const Rcpp::IntegerVector& edge_to) {
  std::vector< std::vector<int> > adj(n);
  for (int e = 0; e < edge_from.size(); ++e) {
    int i = edge_from[e] - 1;
    int j = edge_to[e] - 1;
    if (i < 0 || i >= n || j < 0 || j >= n || i == j) continue;
    adj[i].push_back(j);
    adj[j].push_back(i);
  }
  for (int i = 0; i < n; ++i) {
    std::sort(adj[i].begin(), adj[i].end());
    adj[i].erase(std::unique(adj[i].begin(), adj[i].end()), adj[i].end());
  }
  return adj;
}

static void fill_bfs_distances(const std::vector< std::vector<int> >& adj,
                               int source,
                               int column,
                               Eigen::MatrixXd& distances) {
  int n = adj.size();
  double inf = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) distances(i, column) = inf;
  if (source < 0 || source >= n) return;
  std::vector<int> queue(n);
  int head = 0, tail = 0;
  queue[tail++] = source;
  distances(source, column) = 0.0;
  while (head < tail) {
    int cur = queue[head++];
    double next_dist = distances(cur, column) + 1.0;
    for (int nb : adj[cur]) {
      if (std::isfinite(distances(nb, column))) continue;
      distances(nb, column) = next_dist;
      queue[tail++] = nb;
    }
  }
}

static double native_exp_kernel(double distance, double range, double sigma2) {
  if (!std::isfinite(distance)) return 0.0;
  return sigma2 * std::exp(-distance / range);
}

static double median_positive_anchor_distance(const Eigen::MatrixXd& distances,
                                              const std::vector<int>& anchor0) {
  std::vector<double> values;
  int m = anchor0.size();
  values.reserve(static_cast<size_t>(m) * static_cast<size_t>(std::max(0, m - 1)));
  for (int row = 0; row < m; ++row) {
    int idx = anchor0[row];
    if (idx < 0 || idx >= distances.rows()) continue;
    for (int col = 0; col < m; ++col) {
      double d = distances(idx, col);
      if (std::isfinite(d) && d > 0.0) values.push_back(d);
    }
  }
  if (values.empty()) return NA_REAL;
  size_t mid = values.size() / 2;
  std::nth_element(values.begin(), values.begin() + mid, values.end());
  double med = values[mid];
  if (values.size() % 2 == 0) {
    std::nth_element(values.begin(), values.begin() + mid - 1, values.end());
    med = 0.5 * (med + values[mid - 1]);
  }
  return med;
}

static double row_mean_finite(const Rcpp::NumericMatrix& x, int row) {
  double sum = 0.0;
  int n = 0;
  for (int j = 0; j < x.ncol(); ++j) {
    double v = x(row, j);
    if (!R_finite(v)) continue;
    sum += v;
    ++n;
  }
  if (n == 0) return NA_REAL;
  return sum / static_cast<double>(n);
}

static double row_variance_finite(const Rcpp::NumericMatrix& x, int row) {
  double mu = row_mean_finite(x, row);
  if (!R_finite(mu)) return NA_REAL;
  double ss = 0.0;
  int n = 0;
  for (int j = 0; j < x.ncol(); ++j) {
    double v = x(row, j);
    if (!R_finite(v)) continue;
    double r = v - mu;
    ss += r * r;
    ++n;
  }
  if (n < 2) return NA_REAL;
  return ss / static_cast<double>(n - 1);
}

// [[Rcpp::export]]
Rcpp::List alfak2_native_kriging_cpp(int n_nodes,
                                     Rcpp::IntegerVector edge_from,
                                     Rcpp::IntegerVector edge_to,
                                     Rcpp::IntegerVector anchor_index,
                                     Rcpp::NumericMatrix anchor_values,
                                     Rcpp::NumericVector anchor_sd,
                                     double range = NA_REAL,
                                     double sigma2 = NA_REAL,
                                     double nugget = 1e-4,
                                     bool compute_variance = true) {
  if (n_nodes <= 0) Rcpp::stop("`n_nodes` must be positive.");
  int m = anchor_index.size();
  if (m < 2) Rcpp::stop("Native Kriging requires at least two anchors.");
  if (anchor_values.ncol() != m) Rcpp::stop("`anchor_values` must have one column per anchor.");
  if (anchor_values.nrow() < 1) Rcpp::stop("`anchor_values` must contain at least one row.");
  if (anchor_sd.size() != m) Rcpp::stop("`anchor_sd` must have one value per anchor.");
  if (!R_finite(nugget) || nugget <= 0.0) nugget = 1e-4;

  std::vector<int> anchor0(m);
  for (int j = 0; j < m; ++j) {
    anchor0[j] = anchor_index[j] - 1;
    if (anchor0[j] < 0 || anchor0[j] >= n_nodes) {
      Rcpp::stop("Anchor index outside graph node range.");
    }
    if (!R_finite(anchor_values(0, j))) {
      Rcpp::stop("The first row of `anchor_values` must contain finite point estimates.");
    }
  }

  std::vector< std::vector<int> > adj = build_undirected_adjacency(n_nodes, edge_from, edge_to);
  Eigen::MatrixXd distances(n_nodes, m);
  for (int j = 0; j < m; ++j) {
    if (j % 16 == 0) Rcpp::checkUserInterrupt();
    fill_bfs_distances(adj, anchor0[j], j, distances);
  }

  double auto_range = median_positive_anchor_distance(distances, anchor0);
  if (!R_finite(range) || range <= 0.0) range = R_finite(auto_range) && auto_range > 0.0 ? auto_range : 1.0;
  if (!R_finite(sigma2) || sigma2 <= 0.0) sigma2 = row_variance_finite(anchor_values, 0);
  if (!R_finite(sigma2) || sigma2 <= 0.0) sigma2 = 1e-4;

  Eigen::MatrixXd kpo(n_nodes, m);
  for (int i = 0; i < n_nodes; ++i) {
    for (int j = 0; j < m; ++j) {
      kpo(i, j) = native_exp_kernel(distances(i, j), range, sigma2);
    }
  }

  Eigen::MatrixXd koo(m, m);
  for (int r = 0; r < m; ++r) {
    int node = anchor0[r];
    for (int c = 0; c < m; ++c) {
      koo(r, c) = native_exp_kernel(distances(node, c), range, sigma2);
    }
  }
  for (int j = 0; j < m; ++j) {
    double sd = anchor_sd[j];
    double obs_var = (R_finite(sd) && sd > 0.0) ? sd * sd : nugget;
    koo(j, j) += std::max(obs_var, nugget);
  }

  Eigen::MatrixXd koo_solve = koo;
  Eigen::LDLT<Eigen::MatrixXd> solver;
  solver.compute(koo_solve);
  std::string solve_status = "ldlt";
  if (solver.info() != Eigen::Success) {
    double jitter = std::max(1e-8, nugget * 10.0);
    bool ok = false;
    for (int attempt = 1; attempt <= 5; ++attempt) {
      koo_solve = koo;
      koo_solve.diagonal().array() += jitter;
      solver.compute(koo_solve);
      if (solver.info() == Eigen::Success) {
        ok = true;
        solve_status = "ldlt_jittered";
        break;
      }
      jitter *= 10.0;
    }
    if (!ok) Rcpp::stop("Native Kriging LDLT factorization failed.");
  }

  int n_draw = anchor_values.nrow();
  Rcpp::NumericMatrix predictions(n_draw, n_nodes);
  for (int r = 0; r < n_draw; ++r) {
    double mu = row_mean_finite(anchor_values, r);
    if (!R_finite(mu)) mu = row_mean_finite(anchor_values, 0);
    Eigen::VectorXd y_centered(m);
    for (int j = 0; j < m; ++j) {
      double v = anchor_values(r, j);
      if (!R_finite(v)) v = anchor_values(0, j);
      y_centered(j) = v - mu;
    }
    Eigen::VectorXd alpha = solver.solve(y_centered);
    Eigen::VectorXd pred = Eigen::VectorXd::Constant(n_nodes, mu) + kpo * alpha;
    for (int i = 0; i < n_nodes; ++i) predictions(r, i) = pred(i);
  }

  Rcpp::NumericVector mean(n_nodes), sd(n_nodes);
  for (int i = 0; i < n_nodes; ++i) mean[i] = predictions(0, i);
  if (compute_variance) {
    Eigen::MatrixXd solved = solver.solve(kpo.transpose());
    for (int i = 0; i < n_nodes; ++i) {
      double v = sigma2 + nugget - kpo.row(i).dot(solved.col(i));
      if (!std::isfinite(v)) v = sigma2 + nugget;
      sd[i] = std::sqrt(std::max(v, nugget));
    }
  } else {
    for (int i = 0; i < n_nodes; ++i) sd[i] = NA_REAL;
  }

  return Rcpp::List::create(
    Rcpp::Named("mean") = mean,
    Rcpp::Named("sd") = sd,
    Rcpp::Named("predictions") = predictions,
    Rcpp::Named("range") = range,
    Rcpp::Named("sigma2") = sigma2,
    Rcpp::Named("nugget") = nugget,
    Rcpp::Named("n_anchors") = m,
    Rcpp::Named("auto_range") = auto_range,
    Rcpp::Named("solve_status") = solve_status,
    Rcpp::Named("engine") = "cpp",
    Rcpp::Named("compute_variance") = compute_variance
  );
}
