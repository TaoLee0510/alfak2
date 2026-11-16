#ifndef ALFAK2_SPARSE_JOINT_PAIRWISE_MODEL_H
#define ALFAK2_SPARSE_JOINT_PAIRWISE_MODEL_H

// Sparse joint count model for one population sampled three or more times,
// fitted pair by pair (`series = "pairwise"`).
//
// Each pair of consecutive time points is modelled as two time points are in
// `sparse_joint_model.h`: its own domain (the majors of the pair and their
// one-step neighbours), its own initial abundances at the earlier time point,
// the first-order closed form over its interval, and its own depth. All
// pairs share one landscape: the birth rates and their trend. A time point
// between the first and the last is the end of one pair and the start of
// the next, so the dynamics restart from every time point as it was observed.
//
// States are those of the union of the pairs' domains; `member` says which
// take part in which pair. The initial abundances form a states-by-pairs
// matrix, of which R maps the entries outside a pair's domain away.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sparse_joint_pairwise_model(objective_function<Type>* obj) {
  DATA_MATRIX(y_start);               // counts at each pair's earlier time point, states x pairs
  DATA_MATRIX(y_end);                 // counts at its later time point
  DATA_VECTOR(interval);              // time between the two time points of each pair
  DATA_MATRIX(member);                // 1: the state is in that pair's domain
  DATA_MATRIX(pair_major);            // 1: a major of that pair
  DATA_IVECTOR(edge_from);            // one-step edges between states
  DATA_IVECTOR(edge_to);
  DATA_VECTOR(edge_log_prob);         // log P_ij of the single-daughter kernel
  DATA_VECTOR(log_self);              // log P_ii
  DATA_MATRIX(design);                // trend design; first column is the intercept
  DATA_IVECTOR(coef_group);           // 0 intercept, 1 dosage, 2 curvature, 3 interaction, 4 shell
  DATA_IVECTOR(reference);            // per pair: state whose net rate centres theta
  DATA_MATRIX(zeta_prior_mean);       // per state and pair
  DATA_MATRIX(zeta_prior_sd);
  DATA_SCALAR(intercept_mean);
  DATA_SCALAR(intercept_sd);
  DATA_SCALAR(shell_sd);
  DATA_VECTOR(log_scale_prior_mean);  // dosage, curvature, interaction, residual
  DATA_VECTOR(log_scale_prior_sd);
  DATA_INTEGER(balance);              // 1: non-major states start at mutation-selection balance
  DATA_VECTOR(pair_history);          // per pair: time over which that balance built up
  DATA_SCALAR(balance_sd);
  DATA_MATRIX(balance_offset);        // zero, except to evaluate the expected information
  DATA_VECTOR(weight);                // 1: a state's counts enter the likelihood; 0: held out

  PARAMETER_MATRIX(zeta);             // log expected count at each pair's earlier time point
  PARAMETER_VECTOR(g);                // log birth rate
  PARAMETER_VECTOR(coef);             // trend coefficients
  PARAMETER_VECTOR(theta);            // per pair: log depth ratio, centred on its reference state
  PARAMETER_VECTOR(log_scale);

  int n = y_start.rows();
  int n_pair = y_start.cols();
  Type nll = 0;

  // Landscape prior: trend plus an independent residual, as for two time points.
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

  matrix<Type> balance_residual(n, n_pair);
  balance_residual.setZero();
  matrix<Type> log_intensity(n, n_pair);
  log_intensity.setZero();
  for (int k = 0; k < n_pair; ++k) {
    // Initial abundances of the pair: the non-majors of the pair fed by its
    // majors over the pair's history, as in sparse_joint_model.h.
    vector<Type> log_start(n);
    log_start.setZero();
    vector<int> fed(n);
    fed.setZero();
    if (balance == 1) {
      for (int e = 0; e < edge_from.size(); ++e) {
        int p = edge_from(e), j = edge_to(e);
        if (member(p, k) == Type(0) || member(j, k) == Type(0)) continue;
        if (pair_major(p, k) == Type(0) || pair_major(j, k) == Type(1)) continue;
        Type term = log(Type(2)) + edge_log_prob(e) + g(p) + zeta(p, k) - net(p) * pair_history(k) +
          sparse_joint_log_phi1(net(j), net(p), pair_history(k));
        log_start(j) = fed(j) ? logspace_add(log_start(j), term) : term;
        fed(j) = 1;
      }
    }
    for (int i = 0; i < n; ++i) {
      if (member(i, k) == Type(0)) continue;
      if (fed(i)) {
        balance_residual(i, k) = zeta(i, k) - log_start(i);
        nll -= dnorm(balance_residual(i, k) - balance_offset(i, k), Type(0), balance_sd, true);
      } else {
        nll -= dnorm(zeta(i, k), zeta_prior_mean(i, k), zeta_prior_sd(i, k), true);
      }
      if (weight(i) != Type(0)) nll -= weight(i) * dpois(y_start(i, k), exp(zeta(i, k)), true);
    }

    // The later time point: the closed form over the pair's interval, within its
    // domain.
    vector<Type> log_n(n);
    for (int i = 0; i < n; ++i) log_n(i) = zeta(i, k) + net(i) * interval(k);
    for (int e = 0; e < edge_from.size(); ++e) {
      int i = edge_from(e), j = edge_to(e);
      if (member(i, k) == Type(0) || member(j, k) == Type(0)) continue;
      Type term = log(Type(2)) + edge_log_prob(e) + g(i) + zeta(i, k) +
        sparse_joint_log_phi1(net(j), net(i), interval(k));
      log_n(j) = logspace_add(log_n(j), term);
    }
    Type centre = net(reference(k)) * interval(k);
    for (int i = 0; i < n; ++i) {
      if (member(i, k) == Type(0)) continue;
      log_intensity(i, k) = log_n(i) - centre + theta(k);
      if (weight(i) != Type(0)) nll -= weight(i) * dpois(y_end(i, k), exp(log_intensity(i, k)), true);
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
