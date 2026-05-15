// [[Rcpp::depends(RcppEigen)]]
#include "alfak2_core.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerMatrix alfak2_parse_karyotypes_cpp(Rcpp::CharacterVector labels) {
  return alfak2::labels_to_matrix(labels);
}

// [[Rcpp::export]]
Rcpp::CharacterVector alfak2_stringify_karyotypes_cpp(Rcpp::IntegerMatrix karyotypes) {
  return alfak2::matrix_to_labels(karyotypes);
}

// [[Rcpp::export]]
Rcpp::List alfak2_build_graph_cpp(Rcpp::CharacterVector labels,
                                  Rcpp::IntegerVector y0,
                                  Rcpp::IntegerVector y1,
                                  double beta = 0.00005,
                                  std::string transition_kernel = "exact",
                                  int shell_depth = 2,
                                  int min_cn = 0,
                                  int max_cn = 5,
                                  int max_nodes = 5000) {
  if (shell_depth < 0) Rcpp::stop("`shell_depth` must be non-negative.");
  if (min_cn > max_cn) Rcpp::stop("`min_cn` must be <= `max_cn`.");
  transition_kernel = alfak2::normalize_transition_kernel(transition_kernel);
  Rcpp::IntegerMatrix observed_mat = alfak2::labels_to_matrix(labels);
  int n_obs = observed_mat.nrow();
  int p = observed_mat.ncol();
  if (y0.size() != n_obs || y1.size() != n_obs) {
    Rcpp::stop("`y0` and `y1` must have one entry per karyotype label.");
  }

  std::vector< std::vector<int> > nodes;
  nodes.reserve(n_obs * std::max(1, shell_depth + 1));
  std::unordered_map<std::string, int> id;
  std::vector<int> support_distance;
  std::deque<int> queue;

  for (int i = 0; i < n_obs; ++i) {
    std::vector<int> x(p);
    for (int j = 0; j < p; ++j) {
      x[j] = observed_mat(i, j);
      if (x[j] < min_cn || x[j] > max_cn) {
        Rcpp::stop("Observed karyotype copy numbers must be within [min_cn, max_cn].");
      }
    }
    std::string key = alfak2::make_key(x);
    int new_id = nodes.size();
    id[key] = new_id;
    nodes.push_back(x);
    bool is_seed = (y0[i] + y1[i]) > 0;
    support_distance.push_back(is_seed ? 0 : shell_depth + 1);
    if (is_seed) queue.push_back(new_id);
  }

  std::vector<int> bfs_depth(nodes.size(), 0);
  while (!queue.empty()) {
    int cur = queue.front();
    queue.pop_front();
    int d = bfs_depth[cur];
    if (d >= shell_depth) continue;
    for (int c = 0; c < p; ++c) {
      for (int dir : {-1, 1}) {
        std::vector<int> nb = nodes[cur];
        nb[c] += dir;
        if (nb[c] < min_cn || nb[c] > max_cn) continue;
        std::string key = alfak2::make_key(nb);
        auto it = id.find(key);
        if (it == id.end()) {
          if (static_cast<int>(nodes.size()) >= max_nodes) {
            Rcpp::stop("Graph expansion exceeded `max_nodes`; reduce shell_depth or bounds.");
          }
          int nid = nodes.size();
          id[key] = nid;
          nodes.push_back(nb);
          support_distance.push_back(d + 1);
          bfs_depth.push_back(d + 1);
          queue.push_back(nid);
        } else if (support_distance[it->second] > d + 1) {
          support_distance[it->second] = d + 1;
        }
      }
    }
  }

  int n = nodes.size();
  Rcpp::IntegerMatrix node_mat(n, p);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) node_mat(i, j) = nodes[i][j];
  }
  Rcpp::CharacterVector node_labels = alfak2::matrix_to_labels(node_mat);

  std::vector<alfak2::Edge> raw_edges = alfak2::one_step_edges(nodes, id, beta, transition_kernel);
  std::vector<double> self_weight = alfak2::state_self_weights(nodes, beta, transition_kernel);
  std::vector<int> tr_from, tr_to;
  std::vector<double> tr_weight;
  alfak2::add_row_stochastic_transition(raw_edges, n, tr_from, tr_to, tr_weight, self_weight);

  Rcpp::IntegerVector edge_from(raw_edges.size()), edge_to(raw_edges.size()), edge_chr(raw_edges.size());
  Rcpp::IntegerVector edge_direction(raw_edges.size());
  Rcpp::NumericVector edge_weight(raw_edges.size());
  for (size_t e = 0; e < raw_edges.size(); ++e) {
    edge_from[e] = raw_edges[e].from + 1;
    edge_to[e] = raw_edges[e].to + 1;
    edge_chr[e] = raw_edges[e].chr + 1;
    edge_direction[e] = raw_edges[e].direction;
    edge_weight[e] = raw_edges[e].weight;
  }

  std::unordered_map<std::string, int> ctx_id;
  std::vector<std::string> ctx_labels;
  std::vector<int> ctx_group;
  std::vector<int> parent_from, parent_to, parent_ctx;
  std::vector<double> parent_weight;
  for (const auto& e : raw_edges) {
    if (support_distance[e.to] <= 0) continue;
    if (support_distance[e.from] >= support_distance[e.to]) continue;
    int band = alfak2::ploidy_band(nodes[e.from]);
    std::string ctx_key = std::string(e.direction > 0 ? "gain" : "loss") +
      "_chr" + std::to_string(e.chr + 1) + "_band" + std::to_string(band);
    auto it = ctx_id.find(ctx_key);
    int cid;
    if (it == ctx_id.end()) {
      cid = ctx_labels.size();
      ctx_id[ctx_key] = cid;
      ctx_labels.push_back(ctx_key);
      ctx_group.push_back(e.direction > 0 ? 1 : 0);
    } else {
      cid = it->second;
    }
    parent_from.push_back(e.from);
    parent_to.push_back(e.to);
    parent_ctx.push_back(cid);
    parent_weight.push_back(e.weight);
  }
  if (ctx_labels.empty()) {
    ctx_labels.push_back("fallback");
    ctx_group.push_back(0);
  }

  Rcpp::IntegerVector support(n);
  Rcpp::CharacterVector support_tier(n);
  for (int i = 0; i < n; ++i) {
    int d = support_distance[i];
    support[i] = d;
    if (d == 0) support_tier[i] = "directly_informed";
    else if (d == 1) support_tier[i] = "local_borrowed";
    else if (d == 2) support_tier[i] = "weakly_supported";
    else support_tier[i] = "graph_borrowed";
  }

  Rcpp::IntegerVector obs_index(n_obs);
  for (int i = 0; i < n_obs; ++i) {
    std::vector<int> x(p);
    for (int j = 0; j < p; ++j) x[j] = observed_mat(i, j);
    obs_index[i] = id[alfak2::make_key(x)] + 1;
  }

  return Rcpp::List::create(
    Rcpp::Named("labels") = node_labels,
    Rcpp::Named("karyotypes") = node_mat,
    Rcpp::Named("support_distance") = support,
    Rcpp::Named("support_tier") = support_tier,
    Rcpp::Named("observed_index") = obs_index,
    Rcpp::Named("edge_from") = edge_from,
    Rcpp::Named("edge_to") = edge_to,
    Rcpp::Named("edge_chr") = edge_chr,
    Rcpp::Named("edge_direction") = edge_direction,
    Rcpp::Named("edge_weight") = edge_weight,
    Rcpp::Named("transition_from0") = Rcpp::wrap(tr_from),
    Rcpp::Named("transition_to0") = Rcpp::wrap(tr_to),
    Rcpp::Named("transition_weight") = Rcpp::wrap(tr_weight),
    Rcpp::Named("parent_from0") = Rcpp::wrap(parent_from),
    Rcpp::Named("parent_to0") = Rcpp::wrap(parent_to),
    Rcpp::Named("parent_weight") = Rcpp::wrap(parent_weight),
    Rcpp::Named("parent_context0") = Rcpp::wrap(parent_ctx),
    Rcpp::Named("context_label") = Rcpp::wrap(ctx_labels),
    Rcpp::Named("context_group0") = Rcpp::wrap(ctx_group),
    Rcpp::Named("n_chr") = p,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("transition_kernel") = transition_kernel,
    Rcpp::Named("shell_depth") = shell_depth,
    Rcpp::Named("min_cn") = min_cn,
    Rcpp::Named("max_cn") = max_cn
  );
}
