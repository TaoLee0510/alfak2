#ifndef ALFAK2_SPARSE_JOINT_COHORT_MODEL_H
#define ALFAK2_SPARSE_JOINT_COHORT_MODEL_H

// Sparse joint count model for several two-timepoint samples on one
// landscape.
//
// Every sample keeps its own domain, initial abundances and depth, exactly as
// in `sparse_joint_model.h`. What the samples share is the landscape: a
// "node" is a karyotype, a "state" is a node in one sample's domain, and the
// log birth rate of a state is the node's shared value, optionally with a
// sample-specific deviation. A karyotype that is abundant in one sample and
// an unobserved neighbour in another is then informed by both.

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sparse_joint_cohort_model(objective_function<Type>* obj) {
  DATA_VECTOR(y0);                    // counts, states of all samples concatenated
  DATA_VECTOR(y1);
  DATA_IVECTOR(state_sample);         // sample of each state
  DATA_IVECTOR(state_node);           // landscape node of each state
  DATA_IVECTOR(state_major);          // 1 if the state is major in its sample
  DATA_VECTOR(dt);                    // per sample
  DATA_IVECTOR(edge_from);            // one-step edges between states of a sample
  DATA_IVECTOR(edge_to);
  DATA_VECTOR(edge_log_prob);
  DATA_VECTOR(log_self);              // per state
  DATA_MATRIX(design);                // per node; first column is the intercept
  DATA_IVECTOR(coef_group);           // 0 intercept, 1 dosage, 2 curvature, 3 interaction, 4 shell
  DATA_IVECTOR(reference);            // per sample: state whose net rate centres theta1
  DATA_VECTOR(zeta_prior_mean);
  DATA_VECTOR(zeta_prior_sd);
  DATA_SCALAR(intercept_mean);
  DATA_SCALAR(intercept_sd);
  DATA_SCALAR(shell_sd);
  DATA_VECTOR(log_scale_prior_mean);  // dosage, curvature, interaction, residual, deviation
  DATA_VECTOR(log_scale_prior_sd);
  DATA_INTEGER(deviations);           // 1: states deviate from their node; 0: they equal it
  DATA_INTEGER(balance);
  DATA_VECTOR(history);               // per sample
  DATA_SCALAR(balance_sd);
  DATA_VECTOR(balance_offset);        // zero, except to evaluate the expected information

  PARAMETER_VECTOR(zeta);             // per state: log expected count at the first timepoint
  PARAMETER_VECTOR(m);                // per node: shared log birth rate
  PARAMETER_VECTOR(g);                // per state: log birth rate, used when `deviations`
  PARAMETER_VECTOR(coef);
  PARAMETER_VECTOR(theta1);           // per sample
  PARAMETER_VECTOR(log_scale);

  int n_state = y0.size();
  int n_node = m.size();
  Type nll = 0;

  vector<Type> scale = exp(log_scale);
  for (int k = 0; k < log_scale.size(); ++k) {
    if (k == 4 && deviations == 0) continue;  // no deviation scale under full sharing
    nll -= dnorm(log_scale(k), log_scale_prior_mean(k), log_scale_prior_sd(k), true);
  }
  for (int k = 0; k < coef.size(); ++k) {
    if (coef_group(k) == 0) nll -= dnorm(coef(k), intercept_mean, intercept_sd, true);
    else if (coef_group(k) == 4) nll -= dnorm(coef(k), Type(0), shell_sd, true);
    else nll -= dnorm(coef(k), Type(0), scale(coef_group(k) - 1), true);
  }
  vector<Type> trend = design * coef;
  for (int v = 0; v < n_node; ++v) nll -= dnorm(m(v), trend(v), scale(3), true);

  vector<Type> log_birth(n_state);
  for (int i = 0; i < n_state; ++i) {
    if (deviations == 1) {
      nll -= dnorm(g(i), m(state_node(i)), scale(4), true);
      log_birth(i) = g(i);
    } else {
      log_birth(i) = m(state_node(i));
    }
  }
  vector<Type> net(n_state);
  for (int i = 0; i < n_state; ++i)
    net(i) = (Type(2) * exp(log_self(i)) - Type(1)) * exp(log_birth(i));

  // Initial abundances: free, or at mutation-selection balance with the
  // sample's own majors (see sparse_joint_model.h).
  vector<Type> balance_residual(n_state);
  balance_residual.setZero();
  if (balance == 1) {
    vector<Type> log_start(n_state);
    vector<int> fed(n_state);
    fed.setZero();
    for (int e = 0; e < edge_from.size(); ++e) {
      int p = edge_from(e), j = edge_to(e);
      if (state_major(p) == 0 || state_major(j) == 1) continue;
      Type h = history(state_sample(j));
      Type term = log(Type(2)) + edge_log_prob(e) + log_birth(p) + zeta(p) - net(p) * h +
        sparse_joint_log_phi1(net(j), net(p), h);
      log_start(j) = fed(j) ? logspace_add(log_start(j), term) : term;
      fed(j) = 1;
    }
    for (int i = 0; i < n_state; ++i) {
      if (fed(i)) {
        balance_residual(i) = zeta(i) - log_start(i);
        nll -= dnorm(balance_residual(i) - balance_offset(i), Type(0), balance_sd, true);
      } else {
        nll -= dnorm(zeta(i), zeta_prior_mean(i), zeta_prior_sd(i), true);
      }
    }
  } else {
    for (int i = 0; i < n_state; ++i) nll -= dnorm(zeta(i), zeta_prior_mean(i), zeta_prior_sd(i), true);
  }

  vector<Type> log_n1(n_state);
  for (int i = 0; i < n_state; ++i) log_n1(i) = zeta(i) + net(i) * dt(state_sample(i));
  for (int e = 0; e < edge_from.size(); ++e) {
    int i = edge_from(e), j = edge_to(e);
    Type term = log(Type(2)) + edge_log_prob(e) + log_birth(i) + zeta(i) +
      sparse_joint_log_phi1(net(j), net(i), dt(state_sample(j)));
    log_n1(j) = logspace_add(log_n1(j), term);
  }

  vector<Type> log_intensity1(n_state);
  for (int i = 0; i < n_state; ++i) {
    int s = state_sample(i);
    log_intensity1(i) = log_n1(i) - net(reference(s)) * dt(s) + theta1(s);
    nll -= dpois(y0(i), exp(zeta(i)), true);
    nll -= dpois(y1(i), exp(log_intensity1(i)), true);
  }

  REPORT(net);
  REPORT(log_intensity1);
  REPORT(trend);
  REPORT(balance_residual);
  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
