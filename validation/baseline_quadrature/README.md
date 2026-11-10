# Prior-preserving baseline quadrature prototype

The [review](REVIEW.md) separates numerical posterior accuracy from accuracy
against generating truth. This is a validation prototype; it is not connected
to the default `fit_alfak2()` API.

`quadrature.R` contains the prior-conditioning, proposal, quadrature and mixture
moment functions. `run.R` validates the original Gaussian/lognormal prior law
and evaluates the sixteen fixed 4/8-node references. `analyze.R` compares full
moments/covariance and the recovered original-intercept moments with MCMC.

`pilot.R` and `PILOT_PROTOCOL.md` preserve the original deployment initial-logit
priors on twelve previously revealed 27-node development datasets. The two
processes are invoked with `ALFAK2_PILOT_PARITY=0` and `=1`; each performs actual
conditional fits serially, so at most two fits run concurrently.
`pilot_report.R` requires all dataset/model/distance cells and complete finite
primary risk values before producing any aggregate. It records an explicit
failed development stage, including shape and raw interval-score regressions.

Run commands from the repository root:

```sh
Rscript validation/baseline_quadrature/run.R
Rscript validation/baseline_quadrature/analyze.R
ALFAK2_PILOT_PARITY=0 Rscript validation/baseline_quadrature/pilot.R
ALFAK2_PILOT_PARITY=1 Rscript validation/baseline_quadrature/pilot.R
Rscript validation/baseline_quadrature/pilot_report.R
```

Existing source/case/cache fingerprints are strict. Use the frozen recorded
sources to resume the existing results; preserve these results when exploring
any new version. `ALFAK2_REFERENCE_DLL` selects the required installed native
library. The current pilot runner records the original temporary dataset/V1
cache paths; exact consumed datasets and complete frozen V1 fit objects are
also archived in the twelve `results/pilot_part*/..._inputs.rds` files. This
allows reconstruction if temporary inputs are removed. The complete MCMC
reference cases and compact summaries are retained in the neighboring
`posterior_reference_4to8/results/` directory.

No independent accuracy calibration/test landscape was opened, and no default
was promoted. A numerical posterior improvement is not presented as satisfying
the all-distance fitness/uncertainty objective.
