# Accuracy upgrade: model, validation, and acceptance status

The existing default pipeline remains `pipeline = "two_layer"`, with relative
selection/endpoint transport and corrected Laplace local inference. The new
`pipeline = "joint_count"` is explicit: it has not passed the complete
all-distance point-and-uncertainty acceptance gate. See
[`final_exact_review.md`](../validation/accuracy_upgrade/final_exact_review.md)
for the independently held-out results and their uncertainty.

## Joint model

The candidate fits the original two-timepoint count likelihood once, on the
entire finite target graph. Its log-birth landscape is

\[
g = o + X a + u,\quad a\sim N(0,D),\quad u\sim N(0,Q_u^{-1}),\quad b_i=e^{g_i}.
\]

`X` contains an intercept, chromosome dosage and quadratic effects, and shrunk
pairwise interactions. `Q_u` is a proper graph precision including an identity
term. Prior scales are dimension-normalized; the selected defaults were fixed
using development landscapes before calibration and final testing. For small
graphs, `o_i=-Var(g_i)/2` keeps marginal prior mean birth rate at one. On larger
graphs the residual variance diagonal uses a fixed-seed stochastic estimator;
the correction is explicitly marked approximate.

For fixed beta, the full single-daughter kernel `P` is retained on graph nodes,
including simultaneous chromosome changes. It is not row-normalized. The mean
division process is

\[
\frac{dn}{dt}=[2P^T\operatorname{diag}(b)-\operatorname{diag}(b)]n.
\]

Loss outside the assumed viability bounds and omission from the computational
domain are reported separately. Viability classification is conditional on
those bounds; it does not establish biological mortality from count data.
Explicit target queries avoid complete shell enumeration, but their omitted
intermediate states remain computational outflow. A larger buffer domain can
still be necessary.

Positive Strang splitting uses exact diagonal growth and a fourth-order
nonnegative off-diagonal exponential action. Numerical rescaling preserves
log total growth. Joint fits check step halving at the fitted mode, compare the
independent R and C++ forward implementations, and refine if required.
Failure to meet tolerance invalidates the reported intervals.

Actual `population_sizes=c(N0,Nb)` may supply a Gaussian observation of the log
population ratio. Sampling depths are not population sizes. Without additional
growth information, the absolute birth baseline can remain prior dependent.

## Inference and uncertainty

Latent initial abundances and landscape effects are marginalized using TMB's
Laplace approximation. Outer and conditional gradients, a positive-definite
Hessian, and scale boundaries are checked. Joint MAP remains an explicit
ablation. Fixed scales are mapped out of optimization, and untrusted local
intervals are missing with diagnostic reasons.

Under a Gaussian log-rate approximation with mean `mu` and covariance `C`, the
reported birth moments are

\[
E[b_i]\simeq e^{\mu_i+C_{ii}/2},\qquad
Cov(b_i,b_j)\simeq E[b_i]E[b_j](e^{C_{ij}}-1).
\]

These are approximate moments, not exact posterior integration. Independent
calibration also addresses discrepancies such as finite-population dynamics.
Fixed beta, time intervals, and user-fixed prior scales remain assumptions.

Optional scale averaging uses normalized Gaussian priors and conditional
Laplace evidence on the selected graph. It requires Laplace inference, one
shared DM concentration, and an unweighted normalized likelihood. Graph and
initial-abundance priors depend on observed support, so the weights are
conditional model scores, not probabilities integrating that selection
mechanism. Total covariance includes within- and between-component variation.
Raw mixture intervals use their mixture CDF; relative contrasts require joint
component draws. Dense analytic Gaussian marginalization is an opt-in feature
limited to 32 nodes after a larger resource probe failed.

The corrected two-layer graph solver uses predictive Gaussian log scores,
supports correlated local anchors, and exposes its conditional-on-local-fit CV
scope. Passing covariance does not remove the local prior a second time; this
two-layer approximation is labeled in diagnostics. Polynomial and tree
adapters now propagate input noise and distinguish training from holdout error.

## Evaluation and gate

Every estimate and truth uses one prespecified t0-count-weighted observed-state
reference contrast across the entire target domain. There is no test-truth
affine correction or separate centering at each distance. Intervals use the
same contrast's covariance or joint draws.

The workflow freezes independent landscape, population, and sampling seeds.
Development uses six landscapes with four repeats each. A separate six
landscapes fit empirical uncertainty calibration; twelve further landscapes
with four repeats each provide the final test. Exact continuous-time Gillespie
simulation removes the old tau-leap step bias from the formal evaluation.
The original tau-leap experiments and all regressions remain available.

Calibration changes only interval endpoints, using landscape-balanced signed
standardized errors. Model/scenario fingerprints and distance support must
match. Six calibration landscapes do not provide a finite-sample guarantee.
Both baseline and candidate receive their own independently fitted calibration
for fair uncertainty comparisons.

The final gate requires evidence across exact distances, complete prediction
coverage, and calibrated uncertainty; a lower pooled RMSE is insufficient.
The current final test shows lower mean point error at all tested distances,
but paired confidence intervals cross zero and calibrated far-field interval
scores do not uniformly improve. The candidate therefore remains opt-in.
The 22-chromosome query probes establish numerical execution only, not
high-dimensional statistical accuracy or practical 5,000-node inference.

## Second-round development checkpoint

Two experimental controls were fixed before the second development screen.
`interaction_support="observed_nullspace"` rotates the pair-feature basis into
directions visible at observed rows and directions vanishing there, with a
separate `unidentified_interaction_sd` for the latter. This is a support-adaptive
prior: the support is selected from positive counts at either timepoint.
Its linear design nullspace is not a likelihood-identifiability result.
Normalization and mutation flow can still inform those directions. Conditional
Laplace component weights do not integrate the support-selection mechanism.

Optional `joint_control$prior_components` supplies a list of prior-control
lists, with `model_prior` as their prior weights. It cannot be combined with
`model_scales`. Candidate A changes only the support-aware interaction prior;
candidate B also averages two fixed residual-prior components. Neither receives
generating-family labels, true initial frequencies or true population growth.
The diagnostic-only `oracle_initial_probability` interface fixes known initial
frequencies exactly, including structural zeros; oracle results are kept
separate from fair count-only algorithm comparisons.

The [18-landscape screen](../validation/accuracy_round2/development_stage_review.md)
finds no credible single winner. B lowers relative RMSE by 0.35–2.37% across
distances 0–3, but centered shape risk rises in every distance mean and two
family/shape safety cells exceed 5% harm descriptively. Its reduction in common
reference error accounts for 101.61% of total MSE improvement because shape MSE
increases. This is a development stopping decision, not an independent-test
significance claim. No new calibration or test data are opened.

The empirical calibration experiment preserves native posterior means/SDs and
selects interval-only controls using development landscapes. Missing strata or
unsupported held-out intervals cannot be dropped during selection. Selected
LOO scores are not nested OOF estimates of the complete selection procedure.
The prospective power helper requires explicit complete metric/distance vectors
and unique landscape keys; its output covers only declared core risk/shape
gates, not a certified full promotion probability. No final power computation
or test sample-size certification is claimed after this failed stage.

## Independent numerical posterior checkpoint

The [4/8-node reference study](../validation/posterior_reference_4to8/REVIEW.md)
passes all sixteen fixed-data MCMC precision checks and 99 independent
target/sampler assertions. Under its explicitly fixed initial-logit priors,
the original count-only log-Gaussian/Laplace approximation overstates relative
fitness SD by factors 1.42–3.28. Supplying population-growth information changes
the target and largely removes that discrepancy; it is a mechanism probe, not
a count-only algorithmic improvement.

A [scalar-integration prototype](../validation/baseline_quadrature/REVIEW.md)
conditions on the existing independent Gaussian intercept and integrates its
original prior, preserving the original offset, remaining fitness prior and
initial-frequency treatment. Component-wise mixture moments reduce the
reference relative-SD ratios to 0.923–1.035; two rare-state cases retain material
mean discrepancies. Gaussianizing the full mixture again would lose the
improvement. It remains an approximation, not exact posterior integration.

On twelve unchanged 27-node development samples with deployment eta priors,
relative RMSE improves at d0–d2 but shape error rises slightly at d1/d2. Raw
proper interval scores regress in smooth-d0 and rough-d2. No calibrated
all-distance comparison or new independent test is claimed. The prototype is
not connected to the default API, and numerical agreement alone does not
resolve structural-prior mismatch or missing information about remote states.

## Fixed-budget observation-design checkpoint

The [observation-design study](../validation/observation_design/REVIEW.md)
now has all sixty frozen 27-node development cases available. A prespecified
control-depth contrast includes all twelve
paired biological units on six independent landscapes: doubling typed cells
from 300 to 600 lowers observed-state mean RMSE by 1.69%, but raises mean RMSE
by 2.28% at one hop and 3.74% at two hops. These are descriptive development
results, not a formal significance or noninferiority claim.

The adaptive scalar-integration prototype adds 65 analytic assertions to the
previous 211 implementation checks and passes all eight small-graph reference
checks. It preserves the Gaussian baseline prior using a logit coordinate and
explicit density/Jacobian, with full-support panels and mandatory endpoint
refinement. These are numerical error estimates, not rigorous global bounds.
All six unavailable cases are restored in at most 266.27 seconds and 0.554 GiB
per case. Three representative 27-node compatibility checks pass; 54 converged
legacy posteriors are retained rather than claiming all sixty were recomputed.
Original failures and one diagnostic layout's legacy 305-second budget failure
remain recorded. Some legacy runs lack individual peak-memory measurements.

The complete fixed-budget comparison rejects both primary candidates for the
all-distance development objective. Against depth600, the shared-anchor layout
lowers mean RMSE by 7.29%/6.59% at one/two hops, while observed-state RMSE rises
2.02%, centered shape risk rises 8.37%, and IS95 rises 13.45%. Both primary
candidates also fail family-level safety conditions. These results describe a
near/remote tradeoff on six development landscapes; they do not establish
significance or calibrated uncertainty. No larger-domain expansion, reserved
calibration/test, default promotion or HPC submission is made. The next proposed
design study would preserve near-state information as a constraint while
optimizing initial compositions and the fixed typing budget.
