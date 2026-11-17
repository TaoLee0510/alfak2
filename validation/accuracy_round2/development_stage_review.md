# Round 2 development stage decision

The 18-landscape 81-node development screen does not identify a credible single winner for the requested all-distance point, shape and uncertainty improvement. No new calibration or test dataset is generated. This is a development stopping decision, not a formal rejection of noninferiority or an independent-test result.

All models use the same unchanged counts, 81 target nodes and observed-reference contrast. Six independent old exact landscapes per family are used, one A1/S1 dataset per landscape. Old calibration/revealed-test data are explicitly development material. A/B each completed 18/18 fits with trusted TMB covariance. Frozen V1 summaries are reused with exact dataset/reference/truth checks.

## Point and shape screening

Percent changes are candidate/V1 minus one; negative is better. These are equal-family means of landscape risks. The 5% target is the prospective planning effect, not a claim that a superiority test must establish a full 5% margin.

| Candidate | Distance | RMSE change | Shape RMSE change | Edge RMSE change | Raw interval-score change |
|---|---:|---:|---:|---:|---:|
| A | 0 | -0.001% | -0.000% | -0.001% | +0.288% |
| A | 1 | -0.001% | -0.000% | -0.000% | -6.506% |
| A | 2 | -0.001% | -0.000% | -0.001% | -18.940% |
| A | 3 | -0.000% | +0.001% | -0.001% | -28.572% |
| B | 0 | -0.351% | +1.116% | -0.919% | -37.342% |
| B | 1 | -2.365% | +0.159% | +0.999% | -32.080% |
| B | 2 | -2.017% | +0.024% | +0.086% | -33.360% |
| B | 3 | -2.210% | +0.006% | +0.040% | -34.790% |

A changes point/shape risk only at numerical scale; its raw interval improvement does not establish improved landscape reconstruction. B lowers global reference-contrast error, while shape/gradient improvement is absent across several distances.

## Descriptive family safety flags

The following point/shape cells exceed the 5% family margin in their observed mean. They are retained as safety warnings; no p-value or formal noninferiority conclusion is inferred from this small development sample.

| Candidate | Family | Distance | Metric | Ratio to V1 |
|---|---|---:|---|---:|
| B | interaction | 0 | spatially_centered_rmse | 1.05243 |
| B | smooth | 0 | edge_gradient_rmse | 1.05697 |

## Reference offset versus shape

For the complete target domain, MSE equals globally centered shape MSE plus the squared mean error. No truth-fitted slope, shift or shell centering is applied to an estimate.

| Candidate | Total MSE change | Shape MSE change | Offset-squared change | Offset fraction of MSE reduction |
|---|---:|---:|---:|---:|
| A | -0.000001521 | -0.000000237 | -0.000001283 | 84.39% |
| B | -0.002258095 | +0.000036423 | -0.002294518 | 101.61% |

B's common-offset reduction exceeds its total MSE gain because its centered shape MSE increases. Wider or better calibrated intervals cannot resolve that point/shape failure.

## Calibration, power and sealed validation

No final all-distance calibrated interval comparison or joint power estimate is claimed. Ordinary landscape-LOO scores used to choose gamma/kappa are development selection scores. An actual joint OOF power vector would require nested outer-landscape exclusion for that choice, plus complete RMSE, interval-score and shape/gradient gates. Because the prerequisite point/shape stage is not credible, this additional computation is deferred. The point/shape paired loss CSV is explicitly marked oof_verified=FALSE to prevent accidental use as the complete joint OOF vector.

The primary eventual promotion remains one intersection-union conjunction at one-sided alpha 0.05, with all predeclared distance benefits and family-by-distance 5% noninferiority. Simultaneous/Holm individual-cell claims remain separate. There is no Bonferroni-based assertion that the task is impossible, and no enlargement of test N to chase significance. Provisional N=36 and the reserved cap N=60 remain uncertified until a suitable single winner and credible OOF loss vectors exist.

New domain-calibration/test seeds remain unused. The 64-node generator is prepared and checked, including honest actual-distance coverage and the separate no-information d5/d6 fixture; no candidate accuracy claim in that domain is made from generator preparation alone.

Full model/interval/family tables are in [the development comparison](development81_review.md). The [protocol](PROTOCOL.md) and machine-readable stage CSVs preserve the stop conditions. No prior weight, nullspace scale, truth generator or test seed is changed in response to this screen.

Generated: 2026-11-04 01:33:02 EST
