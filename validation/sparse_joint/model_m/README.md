# Stage M: integrating out the initial abundance of sparse states (stopped)

Part of the multi-timepoint plan (handoff of 2026-11-13). The user approved
stage M on 2026-11-14 and stopped it the same day, on the evidence below.

## Question

Stage C (`validation/sparse_joint/convergence/`) found that the scale
optimisation of `fit_sparse_joint()` fails where the Laplace approximation over
the latent variables degenerates. On the development benchmark at ff13cf4, 91
of 420 candidate fits were degenerate. Their inner Hessian was nearly singular
along the initial abundance (zeta) of non-major karyotypes with few cells.
Searching the scales again (build cv1) raised the share of accepted candidates
only from 73% to 77%.

Stage M asked whether integrating each non-major zeta out exactly, by
one-dimensional quadrature, and leaving the rest to the Laplace approximation,
would remove the degeneracy. Non-major sources would feed their neighbours at
their balance value.

## Files

- `m0.R`, `m0_report.R`, `reference.R`: step M0, the diagnosis.
  - Includes an R reference of the joint density (equal to TMB's to 1.5e-9)
    and of the per-state marginal over zeta.
  - Results are in `results/m0/` (git-ignored); the report is
    `results/m0/report.md`.
- `m1_quadrature.patch`: the M1 prototype. It adds
  `fit_sparse_joint(minor_initial = "quadrature")`, for two timepoints only.
  - New template `src/sparse_joint_quadrature_model.h`, dispatched from
    `src/sparse_joint_model.h` when the data carry `minor_quadrature = 1`.
  - Trapezoid rule on standard normal nodes, step 0.2 on +-8.
  - The posterior uses the joint model's expected information at log E[n0].
  - The rebuild keeps the zeta map.

## Findings

M0: the 91 degenerate candidates and the 39 candidates cv1 still rejected after
searching again.

- **The flat direction** was at least 90% non-major zeta in all 91 degenerate
  candidates. The most negative term of the curvature along it:
  - in 61 of them, the second-sample mixture residual: "was there, not
    sampled" against "made by a neighbour";
  - in 30, the second-order part of the balance prior (long histories).
- **The rejected candidates**: 37 of the 39 jump. A step of 1e-3 in a scale
  moves the objective by a median of 3.95.
  - At the jump, the zeta of majors also move, by up to 5.2.
  - The stage-M density changes across it too: median 0.46, at most 18.
- **Quadrature**: the prior-centred trapezoid rule is accurate to 1e-6 on the
  states of real fits. Every prior-centred rule fails when the counts lie far
  from the prior (synthetic stress test).

M1: the prototype, tried on two degenerate items. It matched the R reference
to 2e-11, but neither item stopped being degenerate.

| Item | Smallest inner eigenvalue | Flat direction now | Time |
| --- | --- | --- | --- |
| dev L05 sparse, history 360 | 1.4e-5 | zeta of two majors with no first-sample cell (0 -> 4, 0 -> 8) | 12x |
| dev L03 long, history 1080 | 1.3e-4 | g of non-majors with 0-1 then 0 cells | 8x |

Of the 91 degenerate candidates, 76 lie in items that have a major with 0-1
first-sample cells. Among healthy candidates it is 160 of 285.

## Conclusion

The degeneracy comes from karyotypes with few cells, not from the non-major
class. Majors are defined by their share in any one sample, so they can have
few cells too. Integrating one class of latent variable out moves the flat
direction to the next.

Integrating (zeta, g) jointly over every sparse karyotype, majors included,
was not tried. It would take tens of times longer, and nothing guarantees it
would remove the degeneracy.

Meanwhile the degeneracy does not cost accuracy. In C3, the fits that cv1
changed had a median RMSE ratio to ff13cf4 of 1.000 at every distance.

## Decision (user, 2026-11-14)

- Stage M is stopped.
- The prototype is kept here as a patch and removed from the package.
- The package keeps:
  - the convergence check and its labels;
  - no search again (internal option `alfak2.convergence_restarts`, off);
  - the history chosen by the lowest objective, as in ff13cf4 (internal
    option `alfak2.choose_accepted`, off).

  The fits are those of ff13cf4.
- A fit that is not accepted now names its flat direction:
  - `diagnostics$optimisation$flat_direction`: the leading latent variables,
    with karyotype, cells per sample, major and share;
  - `flat_composition`: the share by kind of latent variable;
  - the convergence warning names both.

## Applying the patch

On the package tree of build cv3
(`validation/sparse_joint/multi_timepoint/results/src/alfak2_cv3.tar.gz`):

```sh
git apply validation/sparse_joint/model_m/m1_quadrature.patch
rm -f src/*.o src/*.so   # make does not see the header change
R CMD INSTALL .
```

Checked on 2026-11-14: the patched build gives the same fit as the session's
M1 build on `simulate_sparse_joint_counts()` (history = dt). The objective is
179.53250559392 in both and the summary is identical.
