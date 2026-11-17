# Numerical reference attempt 2: fixed Gaussian proposal, identical target

Attempt 1 reached 160k retained draws per chain for both frozen data cases.
Rhat and bulk/tail ESS passed, but mean MCSE was .02824/.02927 posterior SD
(required <=.025). The interaction SD/covariance MCSE also failed. These are
uncertified references; no approximation-accuracy conclusion is drawn from them.

Keep every data value, likelihood, original prior, material-discrepancy rule and
sampling/convergence cap. Reuse the unchanged independent exact mean target.
The new sampler uses a fixed full-state Gaussian background Q=N(m,S) only for
proposal generation. m and C come from the pooled last 40k draws of each of the
four chains in attempt 1's frozen 80k checkpoints; S=2*C+1e-4*I. No Laplace
comparator mean/covariance, generating truth or error against truth is used.

To preserve the exact original whitened target phi(z)*L(z), elliptical slicing
uses residual log density

    log L(z) - z'z/2 + (z-m)' S^-1 (z-m)/2.

The Gaussian normalizing constants are independent of z and cancel from slice
thresholds. Proposals are m+(z-m)*cos(theta)+chol(S)*nu*sin(theta). This is not a
replacement biological prior. The background remains fixed for all retained
draws. Validate with a deliberately shifted/non-unit Q against analytic prior
moments and independent two-node quadrature, so omission of the correction
cannot pass validation.

Four chains start at dispersed N(m,2.25*S) draws, using fresh independently
seeded numerical streams at 9,800,000+1000*case+retained_stage plus the existing
chain offsets. Use 5000 warmup and the same 20k/40k/80k/160k retained cap, serial
chains within a case and at most two model processes overall. Preserve attempt 1
and its raw chains, target hashes and failed precision diagnostics. Do not pool
its retained draws into attempt 2's certified estimator.
