#ifndef ALFAK2_DIVISION_FORWARD_H
#define ALFAK2_DIVISION_FORWARD_H

// Sparse, positive Strang splitting for the mean of a multitype pure-birth
// process: dn/dt = (2 P' diag(b) - diag(b)) n. Diagonal growth is exact;
// off-diagonal births use a fourth-order nonnegative exponential series.
// Normalizing at each substep is only numerical rescaling; log_total retains
// every scale factor, including genuine death and computational graph outflow.
template<class Type>
struct DivisionForwardResult {
  vector<Type> frequency;
  Type log_total;
};

template<class Type>
DivisionForwardResult<Type> division_forward(
    vector<Type> initial, vector<Type> birth_rate,
    vector<int> from, vector<int> to, vector<Type> probability,
    Type dt, int steps) {
  const int n = initial.size();
  vector<Type> diagonal(n);
  for (int i = 0; i < n; ++i) diagonal(i) = -birth_rate(i);
  for (int e = 0; e < probability.size(); ++e)
    if (from(e) == to(e)) diagonal(from(e)) += Type(2) * probability(e) * birth_rate(from(e));
  Type shift = diagonal.maxCoeff();
  Type h = dt / Type(steps);
  vector<Type> half(n);
  for (int i = 0; i < n; ++i) half(i) = exp(Type(0.5) * h * (diagonal(i) - shift));
  vector<Type> state = initial;
  const bool atomic_dense = n <= 500;
  matrix<Type> off_transition;
  if (atomic_dense) {
    off_transition.resize(n, n);
    off_transition.setZero();
    for (int e = 0; e < probability.size(); ++e)
      if (from(e) != to(e)) off_transition(to(e), from(e)) += Type(2) * probability(e);
  }
  Type log_total = 0;
  for (int s = 0; s < steps; ++s) {
    for (int i = 0; i < n; ++i) state(i) *= half(i);
    vector<Type> term = state;
    vector<Type> accum = state;
    for (int order = 1; order <= 4; ++order) {
      vector<Type> next(n);
      next.setZero();
      if (atomic_dense) {
        matrix<Type> rhs(n, 1);
        for (int i = 0; i < n; ++i) rhs(i, 0) = birth_rate(i) * term(i) * h / Type(order);
        matrix<Type> product = atomic::matmul(off_transition, rhs);
        for (int i = 0; i < n; ++i) next(i) = product(i, 0);
      } else {
        for (int e = 0; e < probability.size(); ++e) {
          if (from(e) != to(e))
            next(to(e)) += Type(2) * probability(e) * birth_rate(from(e)) * term(from(e)) * h / Type(order);
        }
      }
      accum += next;
      term = next;
    }
    state = accum;
    for (int i = 0; i < n; ++i) state(i) *= half(i);
    log_total += h * shift;
    // Algebraically equivalent rescaling only once per chunk avoids repeatedly
    // introducing dense AD dependencies inside the sparse exponential action.
    if ((s + 1) % 16 == 0 || s + 1 == steps) {
      Type z = state.sum();
      log_total += log(z);
      state /= z;
    }
  }
  DivisionForwardResult<Type> out;
  out.frequency = state;
  out.log_total = log_total;
  return out;
}
#endif
