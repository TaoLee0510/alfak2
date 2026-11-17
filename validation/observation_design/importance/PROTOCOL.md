# Joint-posterior correction and resumed mechanism goal

User confirmed target/goal mode on 2026-11-05, America/New_York. Work locally,
at most two model processes. Preserve every failed attempt. No production
default change, reserved calibration/test opening, HPC submission or Git
publication. Numerical improvements retain the original likelihood/prior.

## Estimator frozen before comparison

Marginalize the Gaussian landscape coefficients exactly into g~N(mu,K), with
the original offset and covariance. Whiten g=mu+L*z_g. Initial composition uses
the unchanged eta_basis=2*orthonormal Helmert and z_eta~N(0,I). The target is
phi(z)*L(y|z); no true birth, true composition or reference MCMC is used in the
deployable correction. Oracle arms explicitly condition on their supplied pi0.

Construct a proposal from the joint posterior mode and full Hessian in z.
The proposal mixture is .80 Normal(m,C), .19 t_7 with covariance C, and .01
original standard Gaussian prior. Evaluate normalized mixture density exactly.
Draw a 16,384-point pilot. If its raw importance ESS>=max(1000,20*d) and all
moments are finite, update the proposal once: m=.9*m_weighted+.1*m_mode;
C=.9*C_weighted+.1*C_mode+1e-6*I. Otherwise retain the mode proposal. This
proposal-only pilot is discarded from final estimates; no MCMC comparator
or generating truth enters adaptation. Freeze the final proposal.

Use eight independent iid batches, 4096 samples each initially, extending each
to 8192/16384/32768 only for numerical quality. Reuse all within-batch iid draws
under the same fixed proposal. Raw (unsmoothed) normalized weights determine
all posterior estimates. Pareto k is diagnostic only. Require pooled raw
weight ESS>=10000, k<=.5, mean MCSE<=.005 posterior SD and SD MCSE<=.01 posterior
SD for birth, original-reference, globally centered, edge and initial-logit
contrasts. Independent-batch estimates cross-check influence-function MCSE.
Failed/nonfinite target evaluations cannot be omitted. Retain sample/seed
ledgers and unnormalized log weights. Posterior quantiles are weighted sample
quantiles, not quantiles of a moment-matched Gaussian.

## Numerical reference

Use the two unchanged mechanism_v2 eight-node cases. Uniformly extend both
independent MCMC references to 80k retained per chain, then 160k only if mean
MCSE>.005 SD or SD MCSE>.01 SD (or existing convergence checks fail). Keep
the fixed, already validated Gaussian-reference sampler and its pilot proposal.
The old 20k gate result remains archived. New reference sampling is driven by
precision, never by agreement with the new estimator. Reference randomness is
independent of candidate sampling. Numerical streams are 10,100,000 plus
case*1000+retained-stage for references, 10,500,000 plus case*10000+batch/stage
for candidate evaluations; development generator seeds are declared separately.

Keep the original .05-reference-SD mean and [.9,1.1] SD-ratio material limits.
Use joint candidate/reference MCSE, report both separately, and label boundary
uncertainty inconclusive rather than passed. Covariance, baseline identities,
and 50/80/95% interval posterior masses/quantiles are checked and reported.
First validate analytic Gaussian prior, independent two-node quadrature,
forward action, prior/gauge and proposal density/Jacobian calculations.

## 27-node preflight and mechanism continuation

Select two original opt_single cases before corrected accuracy inspection:
highest original runtime; then most initial zero-count states among the other
cases (ties original runtime then case_id). Fresh complete-case elapsed time
includes proposal fitting, pilot, final sampling and diagnostics, with hard
300-second deadline and <=6 GiB peak RSS. Compilation is shared setup and
recorded separately. Preflight successes are reused only with identical source,
data, model and numerical fingerprints. Do not weaken gates after a failure.

Once the small-graph and preflight gates pass, recompute all 12 A-arm finite
process/unknown-pi0 cases under this correction and the already confirmed B/C/D
diagnostic arms, 12 each: exact mean data/unknown pi0; finite process/true-pi0
oracle; exact mean data/true-pi0 oracle. Total 48 cases including preflights.
Keep actual initial integer population, original t0 observations, endpoint
sampling uniforms, reference and original-control distances paired. Candidate
data are development only. Oracle results are labelled mechanism probes and
are excluded from deployable-method ranking. Any unavailable estimate or
interval prevents complete aggregation. Average biological->landscape->family
->equal families; six landscapes remain six independent landscapes.

Report all original RMSE/Shape/Edge/IS50/80/95, coverage, width, family ratios
and paired landscape differences. Numerical target agreement is distinct from
generating-truth accuracy. A corrected estimator need not cure the original
mean-process likelihood's missing branching noise. If numerical gates fail,
continue diagnosing the target-preserving computation, not changing priors,
datasets, biological likelihoods or acceptance limits.

## Numerical proposal revision after v2 failure (2026-11-05)

Both original raw-importance runs failed, retained under importance_v1/small_cases.
No biological comparison or 27-node fit was started. Replace only the proposal:
train on four independent, data-only corrected elliptical-slice chains, 1500
warmup plus 6000 draws each, initialized from N(mode, 4.5*Cmode), using reference
Gaussian N(mode,2*Cmode). These pilot chains are independent of all comparator
chains and never enter final estimates. Rotate the prior-whitened coordinates
so coordinate one is the mean log-birth direction. Fit cubic conditional means
of the other coordinates, subtract them with a triangular shear (Jacobian 1).
Fit eight equal-mass Gaussian components in this transported coordinate, by
pilot baseline octiles; each covariance is 1.2*(.95*Cbin+.05*Cglobal)+1e-7*I.
The final proposal is .80 equal Gaussian components +.19 t7 with global pilot
covariance in transported coordinates +.01 original standard normal in z.
Density and draw transforms are exact inverses with orthogonal/unit Jacobian.
Retain the original iid batch, MCSE, Pareto-k, ESS and resource gates. The
additional training cost is included in each ordinary 300-second fit. Pilot
convergence is descriptive: unbiased target weighting does not require pilot
chains to be exact posterior samples, but final diagnostics must pass.

Revision v4: retain the original single-component proposal if the initial independent 16,384-point raw pilot has ESS>=8192 and Pareto k<=.4; otherwise use the data-only transport pilot. This avoids fitting a redundant flexible transport to the exactly Gaussian prior and simple two-node target. All final numerical gates remain fixed; the pilot is discarded in both branches. v3 two-node failed k diagnostics despite high ESS, and all its outputs remain archived.

Revision v5: v4 analytic validations passed 103/103, but both small cases failed the unchanged Pareto/MCSE gates. Add coordinatewise monotone piecewise-linear marginal Gaussianization after the cubic shear, using pilot quantiles at .0001, .9999, and Phi(linspace(-3.5,3.5,43)). Extrapolate with the end slopes. Both sampling inverse and density Jacobian are exact; retain all numerical gates, reference/data/prior, sample cap and resource limits.

Revision v6: v5 small-case raw tails still failed. Use independent t7 coordinates (unit covariance) instead of Gaussian coordinates inside each transported body component, preserving its Cholesky covariance. Retain global multivariate t7 and original-prior defensive components. Product densities and draw scales are evaluated exactly. Add refinement stages of 65,536 and 131,072 iid draws per batch; no precision/ESS/Pareto or 300s/6GiB threshold is relaxed. The added stages are numerical refinement, and cannot override the complete-case resource gate.

Revision v7: v6 attained mean/SD precision but failed k (.58/.53). Increase data-only pilot to 4*(2000 warmup +30,000 retained) and replace cubic extrapolation by a natural cubic conditional-mean spline, interior knots at baseline pilot quantiles .1/.3/.5/.7/.9 and boundaries .001/.999; tails extrapolate linearly. Use component covariance 1.3*(.8*Cbin+.2*Cglobal)+1e-7I. These are proposal-only changes to cover rare conditional tails; all target and final diagnostic thresholds remain unchanged.

Revision v8: v7 pooled k passed, but several independent batches remained above .5 (.518/.524 maximum), so it did not pass. Add three orthogonal PCA rotations followed by the same monotone marginal Gaussianization to the proposal transport. Retain exact inverse and cumulative log-Jacobian, validated against finite differences and independent density formulas. Use a single product-t7 body with covariance 1.1*Cglobal+1e-7I after these layers, rather than eight conditional components. This targets remaining non-Gaussian dependence; all pilots remain data-only and final gates remain fixed.

Revision v9: v8 passed both high-precision reference gates (maximum contrast mean differences .01151/.00786 reference SD; SD ratios .99194-1.00585). Before scaling to 27 nodes, reduce memory without changing numerical formulas or randomness: compute observable groups one at a time, collect unused temporaries explicitly, release previous stage summaries, and write sample checkpoints/results uncompressed. Repeat analytic and both small-case checks under the same source before resource preflights.

Revision v10: v9 reproduced both reference gates with identical estimates and lower peak RSS. For graphs >8 nodes, save raw samples once in the exact batch checkpoint and retain its path in the summary result, instead of duplicating the entire raw sample matrix in both files. No sample is discarded; all log weights and numerical summaries remain saved. This bounds worst-case 48-fit raw storage to about 22 GiB instead of 42 GiB. All mathematical code, seeds and gates are unchanged.

Revision v11: one v10 worker failed before inference because concurrent Rcpp cache metadata writes exposed a zero-byte token file. Preserve that failed attempt. Serialize compiled-cache initialization with an owned directory lock, with a 120s explicit timeout, and prewarm the frozen source once before the two worker controllers. This is a runtime initialization fix; no mathematical computation, data, seed or gate changes.

Revision v12: v11's two 27-node preflights failed numerical certification at the 1,048,576-sample cap: rare large weights drove mean/SD MCSE above limits. Preserve all artifacts and do not start the 48-fit batch. For the higher-dimensional proposal, train 4*(2000 warmup+60,000 retained) independent data-only chains and use eight PCA/marginal layers. The transported proposal mixture becomes .50 Gaussian(Cbody)+.30 independent-coordinate t7(Cbody)+.19 multivariate t7(Cglobal)+.01 original prior. This restores efficient sampling of the approximately Gaussian central mass while keeping both forms of defensive tails. All target, sample cap, diagnostic and resource thresholds are unchanged.

Revision v13: v12's finite-difference validation at h=1e-5 crossed a downstream piecewise-linear knot in its eight-layer transport; direct step localization at 1e-6/1e-7/1e-8 converged to the analytic log-Jacobian. Change the independent derivative check to h=1e-7, preserving the 1e-8 comparison tolerance. All estimator code and numerical gates are byte-identical to v12; v12 primary small runs were numerical probes while this validation was unresolved and are retained. Validate v13 before primary gate use.

Revision v14: v13's second 27-node case passed, while the first attained all mean/SD precision limits but failed the unchanged per-batch Pareto gate (max .52425). Before marginal/PCA layers, add conditional full-covariance whitening in 16 equal-pilot-mass baseline bins. For residual coordinates after the natural-spline mean subtraction, Cbin=.98*Covbin+.02*Covglobal+1e-6I. Within each bin use its Cholesky factor, exact inverse and log-Jacobian -logdet(Lbin); the baseline coordinate is unchanged. The map is a bijection almost everywhere, with bin boundaries of measure zero; no point or tail is omitted. This addresses conditional variances/correlations that change with baseline. Keep all gates, sample cap and resource bounds unchanged. Test the derivative at the transported global center, away from fixed prior points that may cross multiple downstream spline knots; retain fixed-point inverse and independent density checks.

Revision v15: v14 passed its second preflight, but the first failed absolute birth moments (mean MCSE .005293 SD; SD MCSE .016784 SD), despite raw ESS251,813 and max k .44735. Stage localization showed the error in birth means/SDs, not initial logits. Add a proposal-only Gaussian component shifted +1.5 SD along the transported local mean-log-birth gradient. Estimate this direction by symmetric .01 perturbations of the inverse transport at its global center, then multiply by Cbody and normalize by sqrt(gradient'Cbody*gradient). Mixture is now .30 original Gaussian+.20 shifted Gaussian+.30 product t7+.19 global multivariate t7+.01 original prior. Draws and all mixture/Jacobian densities remain exact; this deliberately allocates proposal mass to birth/birth-squared integrands without changing the posterior. All gates and resource limits remain unchanged.

Revision v16: v14's first preflight also exceeded the 6GiB RSS gate (8,209,235,968 bytes) during its million-sample diagnostics; the failure remains archived. Replace full-array moment diagnostics with 16,384-row accumulation of the identical weighted means/covariances, influence-function MCSE, covariance MCSE and independent-batch checks. Calculate weighted quantiles one observable column at a time. Validate agreement with the preserved legacy formulas on the exact same samples, including means, all covariance/MCSE matrices and quantiles. Keep v15's absolute-birth tail proposal tilt. No target, samples/seed rule, diagnostic threshold or resource bound is changed. v15 passed103/103 validation checks but had no primary candidate runs before this memory fix.
