#ifndef ALFAK2_SPARSE_JOINT_SERIES_MODEL_H
#define ALFAK2_SPARSE_JOINT_SERIES_MODEL_H

// Sparse joint count model for one population sampled at two or more times.
//
// Everything is as in `sparse_joint_model.h` (landscape prior, net rates,
// initial abundances and the balance prior are the same code), except that
// the forward solution from the first sample is evaluated at every later
// sample time tau_k:
//
//   n_j(tau_k) = exp(d_j tau_k) n_j(0) + sum_i M_ji phi1(d_j, d_i; tau_k) n_i(0),
//
// and each later sample has its own depth offset theta_k. Every sample is
// thereby compared with the first-order solution started at the first
// sample, as a fit to the first and last samples compares the last one; a
// later sample does not restart the dynamics. The Hessian keeps the sparsity
// of the two-sample model: a state still depends only on itself and its
// neighbours. With two samples the objective equals that of
// `sparse_joint_model.h`. `sparse_joint_model.h` is left unchanged so that
// two-sample fits stay exactly as they were; R uses this template for more
// than two samples.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sparse_joint_series_model(objective_function<Type>* obj) {
  DATA_MATRIX(y);                     // counts, one column per sample, first sample first
  DATA_VECTOR(tau);                   // time of each later sample since the first
  DATA_IVECTOR(edge_from);            // one-step edges inside the domain
  DATA_IVECTOR(edge_to);
  DATA_VECTOR(edge_log_prob);         // log P_ij of the single-daughter kernel
  DATA_VECTOR(log_self);              // log P_ii
  DATA_MATRIX(design);                // trend design; first column is the intercept
  DATA_IVECTOR(coef_group);           // 0 intercept, 1 dosage, 2 curvature, 3 interaction, 4 shell
  DATA_INTEGER(reference);            // state whose net rate centres theta
  DATA_VECTOR(zeta_prior_mean);
  DATA_VECTOR(zeta_prior_sd);
  DATA_SCALAR(intercept_mean);        // log baseline birth rate
  DATA_SCALAR(intercept_sd);
  DATA_SCALAR(shell_sd);
  DATA_VECTOR(log_scale_prior_mean);  // dosage, curvature, interaction, residual
  DATA_VECTOR(log_scale_prior_sd);
  DATA_INTEGER(n_major);              // majors are the first n_major states
  DATA_INTEGER(balance);              // 1: non-major states start at mutation-selection balance
  DATA_SCALAR(history);               // time over which that balance built up
  DATA_SCALAR(balance_sd);            // log-scale spread around it
  DATA_VECTOR(balance_offset);        // zero, except to evaluate the expected information
  DATA_VECTOR(weight);                // 1: a state's counts enter the likelihood; 0: held out
  DATA_MATRIX(sample_weight);         // per state and later sample: 1 counts, 0 left out

  PARAMETER_VECTOR(zeta);             // log expected count at the first sample
  PARAMETER_VECTOR(g);                // log birth rate
  PARAMETER_VECTOR(coef);             // trend coefficients
  PARAMETER_VECTOR(theta);            // per later sample: log depth ratio, centred on the reference state
  PARAMETER_VECTOR(log_scale);

  int n = y.rows();
  int n_later = tau.size();
  Type nll = 0;

  // Landscape prior: trend plus an independent residual.
  vector<Type> scale = exp(log_scale);
  for (int k = 0; k < log_scale.size(); ++k)
    nll -= dnorm(log_scale(k), log_scale_prior_mean(k), log_scale_prior_sd(k), true);
  for (int k = 0; k < coef.size(); ++k) {
    if (coef_group(k) == 0) nll -= dnorm(coef(k), intercept_mean, intercept_sd, true);
    else if (coef_group(k) == 4) nll -= dnorm(coef(k), Type(0), shell_sd, true);
    else nll -= dnorm(coef(k), Type(0), scale(coef_group(k) - 1), true);
  }
  vector<Type> trend = design * coef;
  for (int i = 0; i < n; ++i) nll -= dnorm(g(i), trend(i), scale(3), true);

  vector<Type> net(n);
  for (int i = 0; i < n; ++i) net(i) = (Type(2) * exp(log_self(i)) - Type(1)) * exp(g(i));

  // Initial abundances, as in sparse_joint_model.h.
  vector<Type> balance_residual(n);
  balance_residual.setZero();
  if (balance == 1) {
    vector<Type> log_start(n);
    vector<int> fed(n);
    fed.setZero();
    for (int e = 0; e < edge_from.size(); ++e) {
      int p = edge_from(e), j = edge_to(e);
      if (p >= n_major || j < n_major) continue;
      Type term = log(Type(2)) + edge_log_prob(e) + g(p) + zeta(p) - net(p) * history +
        sparse_joint_log_phi1(net(j), net(p), history);
      log_start(j) = fed(j) ? logspace_add(log_start(j), term) : term;
      fed(j) = 1;
    }
    for (int i = 0; i < n; ++i) {
      if (fed(i)) {
        balance_residual(i) = zeta(i) - log_start(i);
        nll -= dnorm(balance_residual(i) - balance_offset(i), Type(0), balance_sd, true);
      } else {
        nll -= dnorm(zeta(i), zeta_prior_mean(i), zeta_prior_sd(i), true);
      }
    }
  } else {
    for (int i = 0; i < n; ++i) nll -= dnorm(zeta(i), zeta_prior_mean(i), zeta_prior_sd(i), true);
  }

  // First sample.
  for (int i = 0; i < n; ++i) {
    if (weight(i) == Type(0)) continue;
    nll -= weight(i) * dpois(y(i, 0), exp(zeta(i)), true);
  }

  // Later samples: the forward solution from the first sample at tau_k.
  matrix<Type> log_intensity(n, n_later);
  for (int k = 0; k < n_later; ++k) {
    vector<Type> log_n(n);
    for (int i = 0; i < n; ++i) log_n(i) = zeta(i) + net(i) * tau(k);
    for (int e = 0; e < edge_from.size(); ++e) {
      int i = edge_from(e), j = edge_to(e);
      Type term = log(Type(2)) + edge_log_prob(e) + g(i) + zeta(i) +
        sparse_joint_log_phi1(net(j), net(i), tau(k));
      log_n(j) = logspace_add(log_n(j), term);
    }
    Type centre = net(reference) * tau(k);
    for (int i = 0; i < n; ++i) {
      log_intensity(i, k) = log_n(i) - centre + theta(k);
      if (weight(i) == Type(0) || sample_weight(i, k) == Type(0)) continue;
      nll -= weight(i) * sample_weight(i, k) * dpois(y(i, k + 1), exp(log_intensity(i, k)), true);
    }
  }

  REPORT(net);
  REPORT(log_intensity);
  REPORT(trend);
  REPORT(balance_residual);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
