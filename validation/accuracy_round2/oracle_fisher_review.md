# Round 2: local expected Fisher information

At the three examined v1 development modes, the 27-node count likelihood has algebraic fitness rank 26, yet the prior-relative amount of learning at 150 sampled cells per timepoint is only 0.54–0.71 effective directions after profiling initial frequencies. Scaling expected sampling information to 5,000 cells raises this diagnostic to 4.09–4.50 directions. These are local information calculations, not additional datasets, fitted methods or accuracy guarantees.

## Expansion points and construction

One prespecified finite-population dataset from the first landscape and first sampling replicate of each smooth, rough and interaction family supplies a frozen v1 conditional mode. Three 27-node count-only fits reproduce the required complete mode; no population-growth observation is added. All fitness priors retain v1 values. The independent forward calculation uses the exact continuous-time mean generator A=(2 P^T−I) diag(b), a converged positive matrix-exponential series, and central differences in g=log(b). Initial softmax derivatives and their propagation into final frequencies are analytic.

For each timepoint t, standardize the probability Jacobian by sqrt(N_t) diag(1/sqrt(pi_t)). Stack the two initial-logit Jacobians into A_eta and the log-birth Jacobians into A_g (the initial-time birth block is zero). An SVD supplies an orthonormal basis Q_eta for the data-visible nuisance-column space. The projected design is R=A_g−Q_eta(Q_eta^T A_g), and the profiled expected multinomial Fisher matrix is F=R^T R. This PSD construction does not subtract nearly equal Fisher blocks or use an unstable inverse of the logit gauge.

Let K be the complete frozen Gaussian prior covariance of g, including the feature coefficients and graph residual, and let K=L L^T. The prior-whitened matrix is (R L)^T(R L), with eigenvalues lambda_j. The reported effective learned directions are sum_j lambda_j/(1+lambda_j): a local Gaussian, prior-relative variance-learning diagnostic. It is not a literal count of globally identifiable parameters. Profiling eta omits its proper prior regularization, so this is not the full fitted posterior precision.

## Information and depth planning

Both timepoint depths are scaled together from 150 to 1,000 or 5,000, keeping the expansion point and probabilities fixed. Fisher information scales linearly with depth; the learned fractions saturate nonlinearly. Fixing pi0 at the same mode removes nuisance uncertainty at the identical expansion point and must increase information in PSD order. Replacing pi0 by its true value also changes the expansion point, and is reported separately rather than called a pure same-point information gain.

| Family | Initial-frequency assumption | Learned dirs, 150 | 1,000 | 5,000 | Shape dirs, 150 | Shape dirs, 5,000 | Data rank |
|---|---|---:|---:|---:|---:|---:|---:|
| smooth | Initial frequencies profiled | 0.707 | 2.257 | 4.495 | 0.325 | 3.806 | 26 |
| smooth | Initial frequencies fixed at same mode | 0.936 | 2.825 | 5.043 | 0.509 | 4.366 | 26 |
| rough | Initial frequencies profiled | 0.591 | 1.878 | 4.094 | 0.271 | 3.535 | 26 |
| rough | Initial frequencies fixed at same mode | 0.792 | 2.377 | 4.774 | 0.417 | 4.224 | 26 |
| interaction | Initial frequencies profiled | 0.541 | 1.883 | 4.234 | 0.262 | 3.619 | 26 |
| interaction | Initial frequencies fixed at same mode | 0.756 | 2.445 | 4.997 | 0.422 | 4.390 | 26 |

| Family | True-pi0 sensitivity: learned dirs, 150 | 5,000 | Local data rank at 150 |
|---|---:|---:|---:|
| smooth | 0.837 | 5.067 | 26 |
| rough | 0.835 | 5.059 | 26 |
| interaction | 0.761 | 4.844 | 26 |

## Physical null and shape directions

With pi0 fixed, one final composition contains only n−1 independent probabilities for n free node birth rates. Consequently the data Fisher rank is at most 26 here; beta>0 does not remove this dimension bound. At beta=0, a common additive birth shift leaves frequencies unchanged, and its log-birth tangent is 1/b, not the all-ones vector. A common log-birth shift instead multiplies all rates and generally changes the evolutionary speed. Mutation can rotate the exact data-null direction into a mixture of baseline and shape effects.

For the shape spectrum, map the common additive-birth direction into prior-whitened coordinates, a=L^(-1)(1/b)/||L^(-1)(1/b)||. Choose an orthonormal basis Q_shape perpendicular to a, then diagonalize Q_shape^T F_white Q_shape. This holds that particular prior-metric absolute coordinate fixed; it is not another fit or an unconditional separation of biological shape and absolute growth. Shape and absolute Rayleigh summaries are coupled and must not be added as independent information counts.

| Family | Weakest direction: physical delta-birth alignment with a constant | Fixed-shape absolute Rayleigh at 150 |
|---|---:|---:|
| smooth | 0.503 | 0.692 |
| rough | 0.619 | 0.576 |
| interaction | 0.493 | 0.452 |

The weakest physical directions therefore mix a common baseline with shape. A positive Rayleigh value when shape is held fixed does not establish absolute-rate identification when all shape directions are free. Full node loadings of the three weakest directions are retained in the CSV.

## Numerical and mathematical checks

Halving the central-difference step changes the log-birth Jacobian by at most 8.33e-11 and the profiled Fisher matrix by at most 3.08e-10 in relative Frobenius norm. Analytic initial/final logit Jacobians agree with finite differences to 5.51e-11. The independent continuous-time mode frequencies agree with the fitted numerical solver to 9.24e-09. Known-pi0 minus profiled information has minimum eigenvalue -9.98e-16, consistent with PSD to numerical precision.

The independent mathematical regression file passes 16 assertions: analytic beta-zero derivatives and additive-birth null, nuisance gauge/PSD ordering, mutation changing the null without breaking the composition rank bound, and linear information scaling with saturating learned fractions. Machine-readable results and the test transcript are saved beside the information tables.

## Limits

This is local expected information under the specified multinomial forward model at three prior-influenced modes. It omits finite-population process variance, does not use observed-Hessian curvature or posterior multimodality, and provides no global identifiability proof. Whitening and effective learned directions depend on K; increasing a prior variance can increase this relative diagnostic without increasing data information. Depth extrapolation does not simulate new data, refit modes, remove model bias or certify future RMSE/coverage. This diagnostic must not be substituted for independent candidate validation or used to relabel pair-feature support rank as statistical identification.

Generated: 2026-11-04 01:02:00 EST
