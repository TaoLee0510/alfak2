# Fixed-layout mechanism study: posterior-reference gate

The user confirmed this local study on 2026-11-04 (America/New_York).
Keep the selected opt_single composition, original fitness prior, 600-cell
typing budget, nominal-A fitness reference and original-control distance groups.
Do not change the production default or open reserved calibration/test cases.
No HPC submission or Git publication is authorized by this study.

Stage 1 freezes two eight-node CN1:2 three-chromosome analogues of opt_single.
Shell masses from 2.2.2 are .45,0,.35,.20. The two birth profiles are the
previously fixed prior-offset and interaction planning profiles at rate scale
one. The interaction formula is .16*x1*x2-.12*x2*x3+.12*x1*x2*x3, x=k-2.
Initial population is 1200, beta=.005, T=1.5. Simulate the exact correlated-
daughter Gillespie process; observe non-depleting multinomial counts 240/360.
The reference weight is the original nominal-A analogue: .85 on 2.2.2 and
.15 split equally among its three one-hop neighbors. No case or seed is
replaced for its observed counts or posterior result.

Generation seeds are 9,100,000+100*case+1; sampling seeds are
9,200,000+100*case+1. MCMC seeds are 9,300,000+1000*case+retained_stage,
with the existing distinct per-chain offsets. Validation seeds start at
9,400,000. All are numerical-reference/development namespaces.

Use the existing blocked elliptical-slice sampler and Gaussian-prior whitening,
four dispersed chains per case. Chains execute serially within each case, with
at most two case/model processes in total. 5,000 warmup; 20k retained per chain,
extend only for convergence to 40k,80k,160k. No thinning. Rank split Rhat <=1.01,
bulk/tail ESS >=1000, mean/SD MCSE <=.025 posterior SD; covariance MCSE Frobenius
<=.05 of covariance. Monitor initial-logit centered second products as well as
log birth, birth, relative, centered and edge observables. Constants are
identified explicitly rather than removed for inconvenient precision.

The old sampler's Strang forward action is replaced in this isolated reference
by an independently coded positive Taylor/uniformization action of the full
mean generator, normalizing each block. The sampler transition kernel is
unchanged. Validate the new action against the independent exact R mean solver,
the current TMB action and analytic beta-zero propagation before sampling.
Existing Gaussian-prior and two-node quadrature references validate the sampler.
Sampled posterior states must also pass independent forward/log-likelihood
checks. The MCMC target is the original Gaussian-prior, count-only model, not
the generating process and not a fitted Laplace Gaussian.

Compare the unchanged frozen adaptive baseline integrator with the certified
MCMC target: birth/relative/centered/edge mean, SD, covariance and interval
quantiles/masses. The agreed material rule is mean discrepancy >.05 reference
SD or SD ratio outside [.9,1.1], exceeding four corresponding MCSE. Apply it
to relative, centered and edge contrasts; report all birth/log-birth diagnostics
separately. A failed reference convergence or numerical target check is an
uncertified gate, not evidence that an approximation is accurate. Quantile
comparisons are fixed-data posterior checks, not repeated-sampling coverage.

If material discrepancy appears, first locate the numerical source and stop
the process-vs-initial-composition attribution batch. A tighter scalar
quadrature run may distinguish integration tolerance from remaining conditional
Laplace error. A conditional-a0 reference can be used only as a labelled
numerical-localization experiment with the exact original conditional prior.
No change to biological likelihood/prior is made to make the gate pass.

Only if this gate passes, run the approved 2x2 diagnostic on the twelve frozen
opt_single units (six landscapes, two biological replicates): existing finite
process with unknown pi0 (12 cached), exact mean-generated data with unknown
pi0 (12 new), finite process with true-pi0 oracle (12 new), and mean-generated
data with true-pi0 oracle (12 new). For each paired generator, keep the same
realized initial integer population and t0 observations, and reuse endpoint
sampling uniforms. True birth is used only for generation/scoring. True pi0
is supplied only to the explicitly labelled oracle diagnostics, which cannot
be ranked as deployable methods. Preserve reference/distances and all failures.
At most 36 new routine fits, each hard 300 seconds and <=6 GiB. Reference MCMC
uses its declared sampling cap and separately recorded runtime; reference peak
RSS must also be <=6 GiB. At most two model processes, local only.

Final inference remains descriptive development evidence. A small-graph gate
cannot certify a 27-node posterior or establish a frequentist accuracy result.
Deliver frozen input/source ledgers, complete numerical comparisons, the gate
decision, and any supported localization. If gated out, state clearly that the
36-fit attribution batch was not started.
