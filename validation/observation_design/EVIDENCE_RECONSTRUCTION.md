# Deterministic evidence audit after the first recorded worker failure

The frozen v4 worker recorded an integration error for
`design_rough_L01_R1_composition_anchor600`, after 319.578 seconds. It retained
all 59 conditional nodes for quadrature orders 9, 17 and 33. One cached node
(order 33, node 30) has a NaN marginal objective despite finite saved posterior
moments and trusted gradient/Hessian diagnostics. The original failure,
elapsed time and caches are not overwritten. The 300-second gate has failed
for this case; an offline reconstruction cannot turn that into a resource pass.

At the stored mode, direct joint-objective derivatives give maximum gradient
6.03e-13, positive Hessian minimum eigenvalue 2.174, and a finite Laplace
objective 876.4842. The latent covariance reconstructed from the inverse joint
Hessian agrees to 6.94e-18. The same calculation at three finite-evidence nodes
agrees with their stored objectives and covariances. The problem is therefore
in a recorded marginal evidence value, not a nonfinite physical-birth moment.
Its timing is consistent with an elapsed-time interruption caught inside TMB,
but the saved evidence alone does not prove that cause.

`evidence_reconstruct.R` checks every cached positive-weight node at its saved
mode, without optimization or additional fitting. It reconstructs normalized
Laplace evidence as joint negative log density plus half log Hessian determinant
minus d/2 log(2*pi). It requires joint gradient <=.001, a positive Hessian,
latent covariance agreement <1e-8, and agreement <1e-6 wherever original
evidence is finite. It recovers the importance Gaussian from the stored affine
nodes, checks that recovery at every order, and applies the unchanged prior/
proposal ratio and original moment-refinement tolerances.

Missing cached nodes cannot be fabricated or dropped. If available complete
orders do not meet refinement, the case remains unavailable. Deterministically
reconstructed posteriors are stored separately from raw worker outputs. Any
exploratory risk summary must distinguish numerical availability from the
failed computing-budget gate. No truth-dependent change, new calibration/test
case, or larger-domain run is authorized by this repair.
