#include <TMB.hpp>

template<class Type>
Type softplus_floor(Type x, Type floor) {
  return floor + exp(x);
}

template<class Type>
Type halfnormal_log_density_logscale(Type log_x, Type scale) {
  Type x = exp(log_x);
  return log(Type(2.0)) - log(scale) - Type(0.5) * log(Type(2.0) * Type(M_PI)) -
    Type(0.5) * pow(x / scale, Type(2.0)) + log_x;
}

template<class Type>
Type multinomial_loglik(vector<Type> prob, vector<Type> y) {
  Type out = 0.0;
  for (int i = 0; i < prob.size(); ++i) {
    out += y(i) * log(prob(i) + Type(1e-14));
  }
  return out;
}

template<class Type>
Type weighted_multinomial_loglik(vector<Type> prob, vector<Type> y, vector<Type> weight) {
  Type out = 0.0;
  for (int i = 0; i < prob.size(); ++i) {
    out += weight(i) * y(i) * log(prob(i) + Type(1e-14));
  }
  return out;
}

template<class Type>
Type dirichlet_multinomial_loglik(vector<Type> prob, vector<Type> y, Type phi) {
  Type n = 0.0;
  for (int i = 0; i < y.size(); ++i) n += y(i);
  Type out = lgamma(phi) - lgamma(n + phi);
  for (int i = 0; i < prob.size(); ++i) {
    Type alpha = phi * (prob(i) + Type(1e-12));
    out += lgamma(y(i) + alpha) - lgamma(alpha);
  }
  return out;
}

template<class Type>
Type weighted_dirichlet_multinomial_loglik(vector<Type> prob, vector<Type> y, vector<Type> weight, Type phi) {
  Type n_eff = 0.0;
  for (int i = 0; i < y.size(); ++i) n_eff += weight(i) * y(i);
  Type out = lgamma(phi) - lgamma(n_eff + phi);
  for (int i = 0; i < prob.size(); ++i) {
    Type alpha = phi * (prob(i) + Type(1e-12));
    Type y_eff = weight(i) * y(i);
    out += lgamma(y_eff + alpha) - lgamma(alpha);
  }
  return out;
}

template<class Type>
Type weighted_dirichlet_multinomial_score_loglik(vector<Type> prob, vector<Type> y, vector<Type> weight, Type phi) {
  Type n = 0.0;
  Type weighted_n = 0.0;
  for (int i = 0; i < y.size(); ++i) {
    n += y(i);
    weighted_n += weight(i) * y(i);
  }
  Type norm_weight = weighted_n / (n + Type(1e-12));
  Type out = norm_weight * (lgamma(phi) - lgamma(n + phi));
  for (int i = 0; i < prob.size(); ++i) {
    Type alpha = phi * (prob(i) + Type(1e-12));
    out += weight(i) * (lgamma(y(i) + alpha) - lgamma(alpha));
  }
  return out;
}

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y0);
  DATA_VECTOR(y1);
  DATA_INTEGER(n_nodes);
  DATA_IVECTOR(trans_from);
  DATA_IVECTOR(trans_to);
  DATA_VECTOR(trans_weight);
  DATA_IVECTOR(support_distance);
  DATA_IVECTOR(parent_from);
  DATA_IVECTOR(parent_to);
  DATA_VECTOR(parent_weight);
  DATA_IVECTOR(parent_context);
  DATA_IVECTOR(context_group);
  DATA_SCALAR(dt);
  DATA_VECTOR(obs_weight0);
  DATA_VECTOR(obs_weight1);
  DATA_INTEGER(use_observation_weights);
  DATA_VECTOR(eta_prior_mean);
  DATA_VECTOR(eta_prior_sd);
  DATA_SCALAR(anchor_prior_scale);
  DATA_SCALAR(mu_prior_scale);
  DATA_SCALAR(scale_prior_scale);
  DATA_INTEGER(observation_model);
  DATA_INTEGER(observation_weight_mode);
  DATA_SCALAR(dm_concentration);
  DATA_INTEGER(local_parameterization);
  DATA_VECTOR(f_prior_sd_multiplier);
  DATA_VECTOR(residual_prior_sd);
  DATA_VECTOR(f_center_weights);
  DATA_SCALAR(f_centering_weight);
  DATA_SCALAR(fixed_sigma_neighbor);
  DATA_SCALAR(fixed_sigma_anchor);
  DATA_SCALAR(fixed_tau_group);
  DATA_INTEGER(sdreport_mode);

  PARAMETER_VECTOR(eta);
  PARAMETER_VECTOR(f);
  PARAMETER_VECTOR(delta_context);
  PARAMETER_VECTOR(mu_group);
  PARAMETER(log_sigma_neighbor);
  PARAMETER(log_sigma_anchor);
  PARAMETER_VECTOR(log_tau_group);

  Type nll = 0.0;
  vector<Type> f_raw(n_nodes);
  for (int i = 0; i < n_nodes; ++i) f_raw(i) = f(i);
  vector<Type> f_report(n_nodes);
  if (local_parameterization == 1) {
    for (int i = 0; i < n_nodes; ++i) f_report(i) = f(i) / dt;
  } else {
    for (int i = 0; i < n_nodes; ++i) f_report(i) = f(i);
  }

  Type max_eta = eta.maxCoeff();
  vector<Type> exp_eta(n_nodes);
  Type eta_sum = 0.0;
  for (int i = 0; i < n_nodes; ++i) {
    exp_eta(i) = exp(eta(i) - max_eta);
    eta_sum += exp_eta(i);
  }
  vector<Type> pi0(n_nodes);
  for (int i = 0; i < n_nodes; ++i) pi0(i) = exp_eta(i) / eta_sum;

  vector<Type> growth(n_nodes);
  Type max_g = local_parameterization == 1 ? f.maxCoeff() : (dt * f).maxCoeff();
  for (int i = 0; i < n_nodes; ++i) {
    Type exponent = local_parameterization == 1 ? f(i) : dt * f(i);
    growth(i) = pi0(i) * exp(exponent - max_g);
  }
  vector<Type> transported(n_nodes);
  transported.setZero();
  for (int e = 0; e < trans_weight.size(); ++e) {
    transported(trans_to(e)) += trans_weight(e) * growth(trans_from(e));
  }
  Type z1 = 0.0;
  for (int i = 0; i < n_nodes; ++i) z1 += transported(i);
  vector<Type> pi1(n_nodes);
  for (int i = 0; i < n_nodes; ++i) pi1(i) = transported(i) / z1;

  if (use_observation_weights == 1 && observation_model == 1) {
    Type phi = dm_concentration;
    if (observation_weight_mode == 1) {
      nll -= weighted_dirichlet_multinomial_score_loglik(pi0, y0, obs_weight0, phi);
      nll -= weighted_dirichlet_multinomial_score_loglik(pi1, y1, obs_weight1, phi);
    } else {
      nll -= weighted_dirichlet_multinomial_loglik(pi0, y0, obs_weight0, phi);
      nll -= weighted_dirichlet_multinomial_loglik(pi1, y1, obs_weight1, phi);
    }
  } else if (use_observation_weights == 1) {
    nll -= weighted_multinomial_loglik(pi0, y0, obs_weight0);
    nll -= weighted_multinomial_loglik(pi1, y1, obs_weight1);
  } else if (observation_model == 1) {
    Type phi = dm_concentration;
    nll -= dirichlet_multinomial_loglik(pi0, y0, phi);
    nll -= dirichlet_multinomial_loglik(pi1, y1, phi);
  } else {
    nll -= multinomial_loglik(pi0, y0);
    nll -= multinomial_loglik(pi1, y1);
  }

  if (f_centering_weight > Type(0.0)) {
    Type center_num = 0.0;
    Type center_den = 0.0;
    for (int i = 0; i < n_nodes; ++i) {
      center_num += f_center_weights(i) * f_report(i);
      center_den += f_center_weights(i);
    }
    if (center_den > Type(0.0)) {
      Type center_mean = center_num / center_den;
      nll += Type(0.5) * f_centering_weight * center_mean * center_mean;
    }
  }

  Type sigma_neighbor = fixed_sigma_neighbor > Type(0.0) ? fixed_sigma_neighbor : softplus_floor(log_sigma_neighbor, Type(1e-5));
  Type sigma_anchor = fixed_sigma_anchor > Type(0.0) ? fixed_sigma_anchor : softplus_floor(log_sigma_anchor, Type(1e-5));

  if (!(fixed_sigma_neighbor > Type(0.0))) {
    nll -= halfnormal_log_density_logscale(log_sigma_neighbor, scale_prior_scale);
  }
  if (!(fixed_sigma_anchor > Type(0.0))) {
    nll -= halfnormal_log_density_logscale(log_sigma_anchor, anchor_prior_scale);
  }
  if (!(fixed_tau_group > Type(0.0))) {
    for (int g = 0; g < log_tau_group.size(); ++g) {
      nll -= halfnormal_log_density_logscale(log_tau_group(g), scale_prior_scale);
    }
  }

  for (int g = 0; g < mu_group.size(); ++g) {
    nll -= dnorm(mu_group(g), Type(0.0), mu_prior_scale, true);
  }

  for (int c = 0; c < delta_context.size(); ++c) {
    int g = context_group(c);
    Type tau = fixed_tau_group > Type(0.0) ? fixed_tau_group : softplus_floor(log_tau_group(g), Type(1e-5));
    nll -= dnorm(delta_context(c), mu_group(g), tau, true);
  }

  vector<Type> prior_sum(n_nodes);
  vector<Type> prior_weight(n_nodes);
  prior_sum.setZero();
  prior_weight.setZero();
  for (int e = 0; e < parent_weight.size(); ++e) {
    int child = parent_to(e);
    int parent = parent_from(e);
    int ctx = parent_context(e);
    Type w = parent_weight(e);
    prior_sum(child) += w * (f_report(parent) + delta_context(ctx));
    prior_weight(child) += w;
  }

  for (int i = 0; i < n_nodes; ++i) {
    if (support_distance(i) == 0) {
      nll -= dnorm(f_report(i), Type(0.0), sigma_anchor * f_prior_sd_multiplier(i), true);
    } else {
      Type mean_i = Type(0.0);
      Type sd_i = sigma_neighbor;
      if (prior_weight(i) > Type(0.0)) {
        mean_i = prior_sum(i) / prior_weight(i);
      } else {
        sd_i = sigma_anchor * Type(2.0);
      }
      if (support_distance(i) >= 2) sd_i *= Type(1.75);
      sd_i *= f_prior_sd_multiplier(i);
      nll -= dnorm(f_report(i), mean_i, sd_i, true);
      if (residual_prior_sd(i) > Type(0.0) && prior_weight(i) > Type(0.0)) {
        nll -= dnorm(f_report(i), mean_i, residual_prior_sd(i), true);
      }
    }
    nll -= dnorm(eta(i), eta_prior_mean(i), eta_prior_sd(i), true);
  }

  if (sdreport_mode == 1) {
    int n_direct = 0;
    for (int i = 0; i < n_nodes; ++i) {
      if (support_distance(i) == 0) n_direct++;
    }
    vector<Type> f_direct(n_direct);
    int pos = 0;
    for (int i = 0; i < n_nodes; ++i) {
      if (support_distance(i) == 0) {
        f_direct(pos) = f_report(i);
        pos++;
      }
    }
    ADREPORT(f_direct);
  } else if (sdreport_mode == 2) {
    int n_direct_borrowed = 0;
    for (int i = 0; i < n_nodes; ++i) {
      if (support_distance(i) <= 1) n_direct_borrowed++;
    }
    vector<Type> f_direct_borrowed(n_direct_borrowed);
    int pos = 0;
    for (int i = 0; i < n_nodes; ++i) {
      if (support_distance(i) <= 1) {
        f_direct_borrowed(pos) = f_report(i);
        pos++;
      }
    }
    ADREPORT(f_direct_borrowed);
  } else if (sdreport_mode == 3) {
    ADREPORT(f_report);
    ADREPORT(pi0);
    ADREPORT(pi1);
  } else if (sdreport_mode == 4) {
    ADREPORT(pi0);
    ADREPORT(pi1);
  }
  REPORT(sigma_neighbor);
  REPORT(sigma_anchor);
  REPORT(pi0);
  REPORT(pi1);
  REPORT(f_report);
  REPORT(f_raw);
  REPORT(f_center_weights);
  return nll;
}
