#ifndef ALFAK2_SPARSE_JOINT_MODEL_H
#define ALFAK2_SPARSE_JOINT_MODEL_H

// Joint count model on a sparse karyotype domain.
//
// Division dynamics are dn/dt = (diag(d) + M) n, with d_i = (2 P_ii - 1) b_i
// the net growth rate of state i and M_ji = 2 P_ij b_i the mutational flux
// from i to a one-step neighbour j. The forward model keeps exp(diag(d) t)
// exact and M to first order:
//
//   n_j(T) = exp(d_j T) n_j(0) + sum_i M_ji phi1(d_j, d_i; T) n_i(0).
//
// A state therefore depends only on itself and its neighbours, which keeps
// the Hessian of the latent variables sparse on domains of thousands of
// states. The two multinomial samples are written in their Poisson form for
// the same reason: a softmax over the domain would couple every pair of
// states. `R/forward-closed-form.R` holds the reference implementation.

// log phi1(a, b; T) = log integral_0^T exp(a (T - s) + b s) ds. The absolute
// rate difference is smoothed at zero so that the function is differentiable
// everywhere; the smoothing also keeps x >= 1e-6, where 1 - exp(-x) loses at
// most ten digits to cancellation.
template<class Type>
Type sparse_joint_log_phi1(Type a, Type b, Type dt) {
  Type gap = (a - b) * dt;
  Type x = sqrt(gap * gap + Type(1e-12));
  return Type(0.5) * (a + b) * dt + Type(0.5) * x + log(dt) + log((Type(1) - exp(-x)) / x);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type sparse_joint_model(objective_function<Type>* obj) {
  DATA_VECTOR(y0);                    // counts at the first timepoint
  DATA_VECTOR(y1);                    // counts at the second timepoint
  DATA_SCALAR(dt);
  DATA_IVECTOR(edge_from);            // one-step edges inside the domain
  DATA_IVECTOR(edge_to);
  DATA_VECTOR(edge_log_prob);         // log P_ij of the single-daughter kernel
  DATA_VECTOR(log_self);              // log P_ii
  DATA_MATRIX(design);                // trend design; first column is the intercept
  DATA_IVECTOR(coef_group);           // 0 intercept, 1 dosage, 2 curvature, 3 interaction, 4 shell
  DATA_INTEGER(reference);            // state whose net rate centres theta1
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

  PARAMETER_VECTOR(zeta);             // log expected count at the first timepoint
  PARAMETER_VECTOR(g);                // log birth rate
  PARAMETER_VECTOR(coef);             // trend coefficients
  PARAMETER(theta1);                  // log depth ratio, centred on the reference state
  PARAMETER_VECTOR(log_scale);

  int n = y0.size();
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

  // Dynamics, in log space.
  vector<Type> net(n);
  for (int i = 0; i < n; ++i) net(i) = (Type(2) * exp(log_self(i)) - Type(1)) * exp(g(i));

  // Initial abundances. Free by default. With `balance`, a non-major state
  // starts at what the majors next to it would have fed it over `history`,
  // had each major grown at its present rate throughout:
  //
  //   n_j(0) = sum_p M_jp exp(-d_p h) phi1(d_j, d_p; h) n_p(0),
  //
  // which tends to M_jp n_p(0) / (d_p - d_j) for a long history and a less
  // fit neighbour. Its count in the first sample then informs its fitness.
  // `balance_offset` shifts the residual of that prior. It is zero in a fit.
  // Set to the residual at the mode, it leaves the Gauss-Newton part of the
  // prior's curvature, which is positive semi-definite, as replacing counts
  // by their expectations does for the likelihood.
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

  vector<Type> log_n1(n);
  for (int i = 0; i < n; ++i) log_n1(i) = zeta(i) + net(i) * dt;
  for (int e = 0; e < edge_from.size(); ++e) {
    int i = edge_from(e), j = edge_to(e);
    Type term = log(Type(2)) + edge_log_prob(e) + g(i) + zeta(i) +
      sparse_joint_log_phi1(net(j), net(i), dt);
    log_n1(j) = logspace_add(log_n1(j), term);
  }

  // Counts. theta1 absorbs the depth ratio and the common growth factor.
  Type centre = net(reference) * dt;
  vector<Type> log_intensity1 = log_n1 - centre + theta1;
  for (int i = 0; i < n; ++i) {
    if (weight(i) == Type(0)) continue;
    nll -= weight(i) * dpois(y0(i), exp(zeta(i)), true);
    nll -= weight(i) * dpois(y1(i), exp(log_intensity1(i)), true);
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
