# Frozen numerical localization after the failed primary gate

The certified 20k-per-chain corrected reference flags one relative-mean cell
in n8_layout_prior: 0.0540785 reference SD, 7.72 mean MCSE from zero. This is
close to the prespecified .05-SD material boundary. It is a descriptive
point-estimate flag, not proof that the true approximation error exceeds .05 SD.
The second case passes the relative/centered/edge gate. Original gate and all
results are retained; the 36-fit attribution batch is not started.

For both fixed cases, tighten scalar quadrature tolerances by a factor of ten,
retaining the identical prior/likelihood and complete-support/tail checks.
Compare resulting changes with each case's certified MCMC reference.

Also condition on the original independent intercept a0=0. Preserve the original
offset mu and remove precisely sd(a0)^2*11' from K. This is an exact conditional
prior of the existing model, not an alternative biological prior. Compare the
direct conditional Laplace component with a separate MCMC reference of this
same conditional target, keeping pi0 unknown with the unchanged SD=2 prior.

Use 10k-per-chain blocked ESS as proposal-only pilot (last 5k used), then the
validated density-corrected Gaussian-reference sampler with S=2*C+1e-4*I.
Four dispersed chains, 5k warmup, 20k/40k/80k/160k cap and all original
Rhat/ESS/MCSE criteria apply. Pilot samples are not pooled into reported
reference estimates. Numerical streams use 9,710,000+1000*case for pilots and
9,700,000+1000*case+retained_stage for final conditional chains.

These conditional references locate possible remaining Laplace error at a
single declared intercept. Agreement there cannot certify all scalar nodes
or their Laplace evidence weights. No generating truth or known true initial
frequency is supplied to these numerical posterior calculations.
