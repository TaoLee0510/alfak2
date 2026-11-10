# Fixed-cell-budget observation design: development simulation

This stage is the confirmed comparison of depth, timepoints and initial
composition. It is simulation only. Production functions, old archives and
reserved accuracy calibration/test seeds remain unchanged.

## Frozen designs and costs

The target is the full 3-chromosome CN1–3 domain (27 states), beta=.005,
T=1.5 and total initial population 1,200. The common fitness landscape is shared
across trajectories; each trajectory has its own unknown initial composition.

| Design | Trajectories | Observations | Typed cells | Initial cells |
|---|---:|---|---:|---:|
| control300 | 1 | t=0,T: 150 each | 300 | 1,200 |
| depth600 | 1 | t=0,T: 300 each | 600 | 1,200 |
| time600 | 1 | t=0,T/2,T: 200 each | 600 | 1,200 |
| composition600 | 2 | t=0,T: 150 each per trajectory | 600 | 600+600 |

Separate culture count and sampling occasions are reported as additional costs;
equal typed-cell count is not asserted to mean equal total experimental cost.
No population-growth observation is supplied in the primary comparison.
Readouts are non-depleting multinomial draws from latent finite populations,
matching the existing benchmark observation model.

Initial A has .85 mass at 2.2.2 and .15 split over its six one-copy neighbors.
Initial B has .425 at 1.1.1, .425 at 3.3.3, and .15 split over the union of their
six one-copy neighbors. These choices are fixed without truth/fit inspection.
Generate three independent 600-cell branching lineages A1,A2,B. A1+A2 is the
single 1,200-cell population for control/depth/time; A1 and B supply the two
composition trajectories. Superposition is exact under the assumed linear
division process. This coupling shares stochastic variation across designs.

There are six new development landscapes, two per smooth/rough/interaction
family, and two biological replicates per landscape (12 paired units). Fixed
count uniforms couple cell-budget allocations across designs. Family labels
are evaluation metadata only. Landscape seed namespace is 5,100,000 plus
1,000*family+landscape; process and sampling namespaces are 5,300,000 and
5,400,000. No seed is replaced in response to missing distances or poor fits.
Truth uses the existing family formulas with the pre-existing bounded link
exp(log(2)*tanh(g/log(2))), so birth is positive and bounded .5–2. All designs
within a landscape use exactly the same truth.

## Common model and numerical inference

One count likelihood contains every observation. Fitness has the frozen V1
proper Gaussian coefficient-plus-graph prior. The initial logit prior is
independent N(0,2^2) for every state and trajectory, identical across designs
and independent of observed support. A mean-logit gauge is integrated out with
an orthonormal contrast basis. Initial compositions are estimated from counts;
true/nominal initial probabilities are never supplied as model information.

Integrate the original independent Gaussian baseline coefficient as in the
previous numerical prototype, retaining the original mean-one offset and all
remaining priors. Positive-weight conditional nodes must all have trusted
gradient/Hessian/covariance. Orders 9,17,33,65 refine until birth means change
<1e-3 relative to max(1,|mean|), and birth/reference covariance Frobenius changes
<.005. Step-halving checks at component modes check forward accuracy; failed
numerics remain visible. No more than two actual fits run concurrently. The
6-GiB/300-second per-fit gate is checked in the smoke run before a full batch.

## Fair evaluation and decisions

Every scheme uses the same 27 target nodes and the fixed nominal-A reference
w'A birth. Reference weights are a prespecified reporting contrast, not an
initial-frequency oracle. Primary distance groups use each paired unit's
control300 actual observed support and remain fixed across designs. Also report
the planned-A stencil groups and each design's actual observed-support distance.
No truth-fitted affine correction or shell-specific recentering is applied.

Report complete prediction availability, RMSE, full-domain centered shape error,
one-copy edge-gradient error, rank/amplitude and raw proper interval scores
and 50/80/95% coverage. Both means/covariance and joint interval draws propagate
the baseline mixture. Average biological repeats within landscape, then
landscapes within family, then families equally. Independent units are landscapes.

Fisher screening is a local planning diagnostic averaged over multiple
prespecified parameter points; it is not accuracy evidence. Exact Gillespie
data plus actual joint fits determine the design comparison. A prospective 5%
shape benefit and the existing 5% family safety margin are descriptive pilot
triage thresholds, not claims of formal significance. Keep every model/family/
distance cell and its denominator. Missing strata cannot silently become passes.
The 81/64-node stages are considered only if this pilot gives a credible shape
and uncertainty improvement. No independent test N is certified by this pilot.

Fisher planning uses prior-mean, fixed smooth, and fixed interaction log-rate
profiles at rate multipliers .8,1,1.2 (nine points). Nominal initial compositions
are smoothed with a fixed 1e-4 floor only for finite-logit expansion points;
this smoothing is not supplied to the estimator or simulator. Profiling omits
the proper eta prior, so these are data-information planning quantities rather
than actual posterior precision. Doubling depth must double the profiled data
information at the identical point, checked numerically.
