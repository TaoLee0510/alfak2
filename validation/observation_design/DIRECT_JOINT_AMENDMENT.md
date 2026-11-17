# Equivalent joint-Hessian computation, fixed before accuracy comparison

The v4 batch exposes resource failures in nested TMB inference: one time-layout
case peaks at 6.18 GiB, and a shared-anchor case is externally cut at 300 seconds
with a 9.06-GiB peak. Original cases, cached nodes, resource failures and their
offline recovery/availability results are retained. No biological design,
fitness prior, initial-logit prior, data, reference, scoring rule, or accuracy
threshold is changed in response to these failures.

Here every hierarchical scale is fixed and all model parameters are random
effects in the nested implementation. There are zero outer parameters. Its
conditional Laplace approximation can therefore be computed directly: optimize
the same normalized joint density, form its full Hessian H, and use
`joint_nll + .5*logdet(H) - d/2*log(2*pi)` for marginal negative log evidence.
For linear log-fitness map J, covariance is `J H^-1 J'`. The same joint Gaussian
is transformed to birth moments and mixed over the original Gaussian intercept.
This removes unnecessary nested Laplace/AD-report differentiation, while
preserving the mathematical approximation. It does not make Laplace exact or
add information about generating truth.

`direct_joint.R` retains all existing forward step-halving/backend checks,
normalized priors, reference covariance, 9/17/33/65 quadrature orders and
refinement tolerances. A finite objective, joint gradient <=.001, positive
joint Hessian and finite birth moments are required. Optimizer return flags
are retained; stationarity and curvature determine trust, as for the original
zero-dimensional outer optimization. No scales or initial compositions are
treated as known.

Before a new full run, eight fixed small-graph combinations compare the old
and direct computations: 4/8 nodes, beta 0/.005, one/two trajectories. Latent
means/covariance, evidence and predicted counts must agree within fixed
tolerances (1e-5 for moments/predictions, 1e-6 for evidence). The prior 134
target/data checks remain applicable because the TMB joint objective is
unchanged. A fresh integrated smoke uses the previous 9.06-GiB shared-anchor
case, under the original 6-GiB/300-second gate.

Only if both checks pass may all sixty unchanged cases be recomputed, freshly,
with the direct method. No old fit cache is used to claim new runtime. External
hard deadlines apply to each fresh R child, peak RSS is recorded, and at most
two actual fits run concurrently, including any still-running v4 fits. The
same statistical development rules apply. The full v4 report is archived as
`REVIEW_v4.md`; production/default APIs and reserved calibration/test remain
unchanged. Larger domains remain conditional on all-distance statistical and
resource success.

The first direct implementation was blocked by validation (34/48 assertions):
eight failures compared a non-statistical R `logarithm` attribute, and six
showed small residual mode error, up to about 1.8e-5. Newton polishing was added
to require joint gradient <=1e-6; a subsequent gradient-matrix shape error was
also caught before any full batch. All failed source snapshots, validation
tables and logs remain. Version 3 passes 48/48 checks with no warnings/errors.
The validation compares explicit absolute errors under the stated tolerances;
it does not compare an R determinant attribute. No accuracy results were
inspected to select these repairs.

Later checks expose overflow in unnormalized mean-action intermediates and
failed high-scale conditional optimization. Independent mean checks now factor
common growth, normalize in blocks, and retain log growth. The Taylor action
normalizes every step; 8/32/64 terms are chosen using the fixed intercept node,
with unchanged forward tolerances. Node visits proceed from the center outward,
using the nearest solved mode only as an initial value; original nodes and
mixture order/weights remain. The current 134+77 implementation assertions pass.
However, version 7 still reaches the hard 300-second limit at the extreme
65-order node near 7.333 (1.15-GiB peak). Thus the full recomputation is blocked;
the current prototype is not a qualified replacement. The next adaptive-scalar
integration proposal in `NEXT_NUMERICAL_PLAN.md` has not been implemented.
