# Conditional growth-scale integration candidate

This implements the conditional second step of the confirmed plan, after the
fixed-data reference reveals material approximation error in multiple cases.
It is one numerical-inference candidate, not a new biological prior. No new
accuracy calibration/test data or family-dependent prior selection is used.

The V1 Gaussian log-birth prior contains an independent intercept coefficient
a0 ~ N(0, baseline_sd^2). Condition on this coefficient, remove only its design
column/prior from the inner fit, and add a0 to the original prior offset. Keep
the original offset (including the intercept variance correction), remaining
feature priors, proper graph residual and initial-logit priors unchanged.
Integrating the original intercept prior exactly therefore recovers the same
joint model. No centering constraint or external population measurement is
introduced. Other common-scale variation in the proper residual remains in
the inner fit; this is not claimed to marginalize every weak direction.

Each inner posterior is the existing checked Laplace approximation. Its
normalized marginal objective supplies the conditional likelihood estimate.
Gaussian importance quadrature includes the original intercept prior/proposal
density ratio. The proposal mean and variance are recovered from the original
Laplace g posterior and the exact Gaussian conditional law of a0 given g.
This proposal chooses evaluation points; it does not replace the target prior.

Quadrature orders are 9, 17, 33, at most 65. Refinement is based only on
numerical convergence: maximum birth mean change <1e-3 relative to max(1,|mean|)
and birth/reference-contrast covariance Frobenius changes <5e-3. Every positive
weight node must have trusted inner covariance; failed nodes are never silently
dropped. Weights and within/between-component covariance remain explicit.
Mixture moments are transformed component by component, not through a single
moment-matched Gaussian in log space.

First validate recovery of the original Gaussian/lognormal prior moments, then
compare against all 16 already fixed 4/8-node references. This is development
of a numerical approximation using revealed mechanism cases; it is not a new
independent accuracy test. Counts-only and growth-information arms remain
separate. Inner Gaussian error may persist, especially for rare initial states.
No default API change or full-graph accuracy promotion follows from this stage.
