#ifndef ALFAK2_SPARSE_JOINT_STATE_SPACE_MODEL_H
#define ALFAK2_SPARSE_JOINT_STATE_SPACE_MODEL_H

// Sparse joint count model for one population sampled three or more times,
// as a state-space model (`series = "state_space"`).
//
// Every state has a latent log expected count at every time point. From one
// time point to the next the population follows the first-order closed form of
// `sparse_joint_model.h` over the interval, started from the latent
// abundances at the earlier time point, plus a Gaussian deviation on the log
// scale: drift through passaging bottlenecks and whatever else the dynamics
// leave out. Its variance is drift * interval / frequency, with the
// frequency of the state about that interval taken from the counts, and
// rises to at most `max_drift_var`. The drift is one scale parameter. The
// counts at each time point enter once, and the landscape (birth rates and trend) is
// shared. With no drift the latent abundances follow the closed form from
// time point to time point; with a lot, each interval is fitted on its own.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sparse_joint_state_space_model(objective_function<Type>* obj) {
  DATA_MATRIX(y);                     // counts, states x time points, first time point first
  DATA_VECTOR(interval);              // time from each time point to the next
  DATA_IVECTOR(edge_from);            // one-step edges inside the domain
  DATA_IVECTOR(edge_to);
  DATA_VECTOR(edge_log_prob);         // log P_ij of the single-daughter kernel
  DATA_VECTOR(log_self);              // log P_ii
  DATA_MATRIX(design);                // trend design; first column is the intercept
  DATA_IVECTOR(coef_group);           // 0 intercept, 1 dosage, 2 curvature, 3 interaction, 4 shell
  DATA_INTEGER(reference);            // state whose net rate centres theta
  DATA_VECTOR(zeta_prior_mean);       // first time point
  DATA_VECTOR(zeta_prior_sd);
  DATA_SCALAR(intercept_mean);
  DATA_SCALAR(intercept_sd);
  DATA_SCALAR(shell_sd);
  DATA_VECTOR(log_scale_prior_mean);  // dosage, curvature, interaction, residual, drift
  DATA_VECTOR(log_scale_prior_sd);
  DATA_INTEGER(n_major);              // majors are the first n_major states
  DATA_INTEGER(balance);              // 1: non-major states start at mutation-selection balance
  DATA_SCALAR(history);               // time over which that balance built up
  DATA_SCALAR(balance_sd);
  DATA_VECTOR(balance_offset);        // zero, except to evaluate the expected information
  DATA_VECTOR(weight);                // 1: a state's counts enter the likelihood; 0: held out
  DATA_MATRIX(drift_weight);          // per state and interval: 1 / frequency
  DATA_SCALAR(max_drift_var);         // largest variance of a deviation
  DATA_MATRIX(transition_offset);     // zero, except to evaluate the expected information

  PARAMETER_VECTOR(zeta);             // log expected count at the first time point
  PARAMETER_MATRIX(zeta_later);       // log expected count at each later time point
  PARAMETER_VECTOR(g);                // log birth rate
  PARAMETER_VECTOR(coef);             // trend coefficients
  PARAMETER_VECTOR(theta);            // per interval: log depth ratio, centred on the reference state
  PARAMETER_VECTOR(log_scale);        // the last is the log drift

  int n = y.rows();
  int n_interval = interval.size();
  Type nll = 0;

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
  Type drift = scale(4);

  vector<Type> net(n);
  for (int i = 0; i < n; ++i) net(i) = (Type(2) * exp(log_self(i)) - Type(1)) * exp(g(i));

  // Initial abundances at the first time point, as in sparse_joint_model.h.
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
  for (int i = 0; i < n; ++i) {
    if (weight(i) != Type(0)) nll -= weight(i) * dpois(y(i, 0), exp(zeta(i)), true);
  }

  // From each time point to the next.
  matrix<Type> transition_residual(n, n_interval);
  for (int k = 0; k < n_interval; ++k) {
    vector<Type> from(n);
    for (int i = 0; i < n; ++i) from(i) = k == 0 ? zeta(i) : zeta_later(i, k - 1);
    vector<Type> log_n(n);
    for (int i = 0; i < n; ++i) log_n(i) = from(i) + net(i) * interval(k);
    for (int e = 0; e < edge_from.size(); ++e) {
      int i = edge_from(e), j = edge_to(e);
      Type term = log(Type(2)) + edge_log_prob(e) + g(i) + from(i) +
        sparse_joint_log_phi1(net(j), net(i), interval(k));
      log_n(j) = logspace_add(log_n(j), term);
    }
    Type centre = net(reference) * interval(k);
    for (int i = 0; i < n; ++i) {
      Type mean = log_n(i) - centre + theta(k);
      // drift * interval / frequency, capped smoothly at max_drift_var.
      Type v = drift * interval(k) * drift_weight(i, k);
      v = v * max_drift_var / (v + max_drift_var);
      transition_residual(i, k) = zeta_later(i, k) - mean;
      nll -= dnorm(transition_residual(i, k) - transition_offset(i, k), Type(0), sqrt(v), true);
      if (weight(i) != Type(0)) nll -= weight(i) * dpois(y(i, k + 1), exp(zeta_later(i, k)), true);
    }
  }

  REPORT(net);
  REPORT(trend);
  REPORT(balance_residual);
  REPORT(transition_residual);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
