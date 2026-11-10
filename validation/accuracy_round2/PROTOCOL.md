# Round 2: development and sealed-validation protocol

Round-1 reports, datasets and archives remain read-only. Its revealed exact calibration/test landscapes are now development material; they cannot certify a method designed after seeing them.

## Comparators and model selection

The primary comparator is frozen V1, with joint-count Laplace inference, isotropic interactions and prior SDs baseline 0.5, main 0.07, curvature 0.04, interaction 0.04, residual 0.125 and edge 0.25. The actual current default two-layer API and an interval-free zero contrast predictor are separate controls. The original pre-upgrade baseline is historical context.

Only two development candidates are considered. A uses the observed-design interaction nullspace and unidentified-interaction SD 0.2. B adds weak/rough residual components with fixed SD pairs (0.03,0.06)/(0.3,0.6) and prior weights (0.75,0.25). Generating-family labels never enter an estimator or select a component. At most one winner may enter the new test. These prior weights/nullspace scales are not changed from a new calibration/test result.

The 27-node finite-population fixture is only an interface/development pilot. The prespecified 81-node development pilot reuses one A1/S1 sample from each old exact L03–L08 landscape: 18 landscapes, six per family. Imported copies are labeled development and retain original dataset/cache SHA256. Frozen V1 summaries are reused only after exact counts, reference weights, truth and native-distance checks. A/B use deployment eta priors; no known-pi0 or total-growth oracle is provided.

## Domains and generation

| Domain | Complete target | Truth | Primary distance | Seed namespace |
|---|---|---|---|---:|
| cn4_1to3 | 4 chromosomes, CN1–3, 81 states | Unchanged round-1 linear truth functions | Actual observed-support distance | 3,100,000 |
| cn3_1to4 | 3 chromosomes, CN1–4, 64 states | Bounded positive log link | Actual observed-support distance | 3,200,000 |

Both domains use the exact Gillespie correlated-daughter process, beta 0.005 per copy, interval 1.5, initial population 1,200 and count depths 150/150. The second domain defines `birth = exp(log(2) * tanh(g/log(2)))`, fixing birth in [0.5,2] and founder birth at 1. The smooth/rough/interaction structures, parameter ranges and seeds are fixed in `truth_and_data.R`. The first domain retains its original linear formula and positivity check; an invalid predeclared truth is a generator failure, never clipped or replaced with another seed.

Primary strata always use the count sample's actual observed support. A fixed founder-plus-shell-1 stencil is an auxiliary view only; it has maximum distance 5 in the 64-node domain. An auxiliary d5 is never advertised as actual d5. Missing actual strata remain missing; seeds are not resampled to produce far-field coverage. A separate predeclared beta-zero/founder-only support-limit fixture contains the actual d5 query 4.4.3. It tests honest unsupported predictions, not a claim of data-identified distant fitness.

## Splits, budget and unsealing

Each domain has its own numeric seed namespace. Development indices are 1–8; new calibration indices 9–14 give 18 landscapes; test indices 15–26 provisionally give 36 landscapes. Indices through 34 are reserved solely to permit a power-based choice up to 60 landscapes **before any new test is opened**. The provisional count is not a power certification.

There are two independent biological trajectories and two count samples per trajectory. Independent units are landscapes. Calibration costs 72 datasets/model/domain; the provisional test costs 144; the maximum test costs 240. No more than two actual fits run simultaneously. At the measured 81-node pilot rates, A is approximately 17–21 seconds/fit and B approximately 33–42 seconds/fit; the 64-node domain must be timed separately before a full budget is accepted.

New calibration cannot be generated before candidate/protocol freeze. New test cannot be generated before model, calibration, analysis gates, final sample size and all relevant SHA256 are frozen. In the present development-only state, `prepare_data.R` refuses every non-development request even if a provisional freeze record is supplied. Releasing that phase lock requires an explicit protocol-state change after the stage gate passes. No new calibration or test truth is generated.

## Risk, shape and uncertainty

All methods share the same target nodes and one t0-count-weighted observed-reference contrast. There is no truth-fitted affine rescaling or per-shell centering. Reported risks average within landscape, then landscapes within family, then families equally.

The hard audit set includes RMSE, whole-domain spatially centered error, one-copy edge-gradient error, rank/correlation, amplitude/variance ratio and proper interval score. Spatial centering removes one error mean over the complete domain before stratifying. Edge risk includes every unordered one-copy edge and uses maximum endpoint support distance. A gain caused only by a common reference offset is insufficient evidence of landscape-shape reconstruction. Amplitude ratios are interpreted together with shape/gradient risk; valid shrinkage is not forced to amplitude 1. Numeric rank/amplitude release thresholds must be frozen from development before test.

The primary promotion is one conjunctive claim: all predeclared distance RMSE/proper-interval-score benefits and all family-by-distance 5% noninferiority requirements must pass. Its component one-sided tests use alpha 0.05 under intersection-union logic. Individual-cell discoveries need separate simultaneous/Holm reporting. A 76/190-cell conservative Bonferroni budget is a reporting sensitivity, not a mandatory correction for the single conjunction.

Both comparator and winner use the same development-only calibration-selection rule. The selected gamma/kappa never modifies the posterior SD or point prediction; it only defines a separate empirical interval scale/CDF. No kappa>0 is presumed beneficial. New calibration data subsequently fit only frozen interval calibration, with scenario/model fingerprints and unsupported-distance checks.

## OOF planning and stop rules

Old rho/equal-variance power curves are budget proxies, not candidate-v2 power. Actual joint pass probability requires complete paired OOF loss vectors, resampling landscapes within family and inflating variance by 1.5. Shape/gradient gates remain substantive checks; aggregate RMSE alone cannot authorize testing. The 64-node far-field domain has no old-domain variance certificate.

When gamma/kappa is selected, actual OOF calibration losses require nested landscape exclusion: choose hyperparameters by inner LOO on the outer training landscapes, fit there, then predict the held landscape. Selecting on all LOO scores and reusing those same scores is development screening, not nested OOF certification. If point/shape feasibility already fails, unnecessary nested calibration computation is skipped and the stage stays closed.

One sample size is frozen before the test and analyzed once. No model update, prior-weight change, extra repeats or reserved-seed opening is used to chase significance afterward. A failed/inconclusive stage remains failed/inconclusive. Complete prediction/interval availability and absent-stratum denominators are reported separately.

## Completed development checkpoint

The fixed A/B screen on 18 old exact 81-node landscapes completed with 18/18 trusted fits per candidate. A changes point/shape risks only at numerical scale. B improves reference-contrast RMSE by about 0.35–2.36% across d0–d3 but increases spatially centered shape error in every distance mean. Its interaction-d0 shape and smooth-d0 gradient risk ratios exceed 1.05 descriptively. The common-offset reduction explains 101.61% of B's whole-domain MSE gain because shape MSE rises. These are development findings, not formal independent-test significance claims.

No single winner is frozen. Nested calibration and full joint-power computation are deferred under the predeclared point/shape feasibility stop rule. New calibration/test remain sealed; reserved test N remains uncertified. See [the stage decision](development_stage_review.md), [full comparison](development81_review.md), and `results/development81_comparison/stage_decision.csv`. Neither better raw coverage nor a wider interval overrides the landscape-shape requirement.
