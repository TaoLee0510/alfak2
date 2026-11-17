# Near-state constrained observation design: confirmed local development study

Confirmed by the user on 2026-11-04. Local execution only. No production-default
change, HPC submission, reserved calibration/test opening or Git publication.

Keep the 27-node 1:3 three-chromosome graph, beta=.005, T=1.5, initial total 1200,
600 typed cells, original Gaussian landscape and initial-logit priors, nominal-A
fitness reference, and the existing control300 actual-distance groups. All main
fits use the frozen adaptive_source_v1, direct_source_v7 and
direct_forward_source_v2 numerical implementation. The main comparison has six
independent landscapes, two biological replicates each. Repeating a numerical
fit does not increase the independent sample size.

## Candidate grid, frozen before Fisher selection

Single culture: probabilities constant within Manhattan shells 0,1,2,3 from
2.2.2. Shell total masses run over the full simplex in increments .05. Expanded
candidates require shell2+shell3 >= .05. The original A is a planning control.
The initial culture contains 1200 cells. Initial/end typing fractions .4/.6,
.5/.5 or .6/.4 give 240/360,300/300 or 360/240 cells.

Two cultures: A remains the original A. Each culture initially has 600 cells.
Culture B has total mass .15 on the six states neighboring 1.1.1 or 3.3.3,
mass h on 2.2.2, and (.85-h)/2 on each extreme; h in {.15,.30,.45}.
A/B typing budgets are 400/200 or 450/150. Within each culture independently,
initial/end typing fractions are .4/.6,.5/.5,.6/.4. Budgets are integer-valued.

## Truth-blind planning and selection

Reuse the nine previously specified prior/smooth/interaction log-birth profiles
at rate scales .8,1,1.2 from fisher.R, not the six generating truths. Compute
expected multinomial information jointly for log birth and unknown initial
composition. Inference uses eta_basis=2*orthonormal Helmert and zeta~N(0,I),
so include unit precision for zeta and the original full landscape precision.
Marginalize zeta using the joint precision Schur complement. Also report the
no-eta-prior profile-information diagnostic for selected designs.

Evaluate both planning probability floors 1e-4 and 1e-6, renormalizing after
flooring; simulated compositions have no artificial floor or pseudo-cells.
For every planning point, floor and each of the seven originally planned
observed states, fixed-reference birth variance must be <= depth600's value
times (1+1e-6); the 1e-6 is numerical tolerance, not permitted scientific harm.
Minimize the worst ratio of group-mean fixed-reference variance at planned
distances 1 and 2 across all points and floors. Ties within 1e-10 use smaller
mean risk ratio, then lexicographic candidate ID. Select at most one feasible
expanded single-culture candidate and one feasible two-culture candidate.
Never weaken the constraint when a class has no feasible candidate. Fisher
variance is a local planning approximation, not real RMSE or coverage evidence.

Before candidate-data generation, freeze the selected specifications and every
source hash. If no class is feasible, report the planning failure and complete
the independent process-noise diagnostic; do not manufacture a winner.

## Data, fitting and empirical gates

Retain the existing twelve depth600 datasets unchanged, but freshly recompute
their posteriors with the same adaptive integrator as the new candidates.
Candidate process seeds use namespace 8,100,000 plus family*10000 +
landscape*100 + replicate*10 + culture (1/2). This couples matched culture
processes across the two selected designs. Candidate sampling uniforms use
8,200,000 plus the same family/landscape/replicate offsets, one matrix per
biological unit shared across candidate layouts. Only development truths are
reused, never supplied to inference. No seed replacement or rejection based
on realized support, estimates or scores. Preserve original reference and
baseline distances and also report native/planned distances.

At most 36 main cases, with exactly two single-child controllers at most.
Hard 300-second wall deadline; peak RSS <=6 GiB measured by wait4. Preflight
cases are part of the main comparison, and successful results are reused.
Any failure remains in the denominator and prevents a complete comparison.

Use the frozen score.R and biological->landscape->family->equal-family
aggregation. RMSE and Edge nonworse at every original distance; globally
centered Shape at least 5% better at every distance; raw IS50/80/95 nonworse;
each family/distance RMSE, Shape, Edge, IS95 ratio <=1.05. Complete predictions,
intervals and numerical/resource gates are all required. Report all six paired
landscape differences and raw coverage as descriptive development evidence.
No formal significance/noninferiority or final-test sample-size claim.

## Separate process-noise diagnostic

At each of the three fixed profiles at scale 1, simulate 100 Gillespie paths
from the exact integer original-A starting population, N=1200. Seeds are
8,300,000 + profile*1000 + replicate. Thus founder allocation is held fixed;
this isolates branching-process noise from initial-allocation variability.
For endpoint population fractions q and a future multinomial sample of 300,
decompose Cov(Y/300) = E[diag(q)-q q']/300 + Cov(q). Compare this conditional
simulation decomposition with a deterministic-mean multinomial approximation,
report means, covariances, count/relative-variance ratios and Monte Carlo
uncertainty. Diagnostic truths are not supplied to an estimator. This does
not fit 300 extra posteriors, calibrate intervals or change current likelihoods.

Deliver a protocol/source/data ledger, all planning results, frozen candidates,
complete-case/failed-case results and a written development decision. If a
candidate passes, further independent and larger-domain verification remains
a subsequent study; a pilot pass is not proof of universal accuracy.
