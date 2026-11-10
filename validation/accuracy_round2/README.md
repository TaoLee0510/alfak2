# Accuracy round 2

This directory adds development diagnostics and protocol preparation without editing round-1 reports/archives. New calibration and test datasets remain sealed. The completed [18-landscape stage decision](development_stage_review.md) finds no credible all-distance point/shape winner; [full 81-node results](development81_review.md), [the n27 pilot](candidate_pilot_review.md), [the protocol](PROTOCOL.md) and machine-readable tables remain inspectable.

Useful entry points:

```sh
Rscript validation/accuracy_round2/diagnose_round1.R
Rscript validation/accuracy_round2/diagnose_uncertainty_ordering.R
Rscript validation/accuracy_round2/power_plan.R
Rscript validation/accuracy_round2/prepare_protocol.R
Rscript validation/accuracy_round2/test_protocol.R
Rscript validation/accuracy_round2/prepare_data.R --domain cn3_1to4 --split development --landscapes 1 --output /tmp/domain64-development
Rscript validation/accuracy_round2/import_old_development.R
Rscript validation/accuracy_round2/reuse_frozen_v1.R
Rscript validation/accuracy_round2/screen_development.R
Rscript validation/accuracy_round2/archive_development.R
```

`run_models.R` accepts `--data`, `--output`, `--models`, either `--source` plus `--dll` or `--library`, optional `--select-mode`/`--limit`, and a shared `--fit-slot-dir`. It records source/DLL, argument, dataset and scoring SHA256 and refuses a third occupied fit slot. A, B, v1_joint, current_default_two_layer and zero_point_predictor are explicit model names. Native support distance overrides process-pair common distance when available. Generating-family and oracle fit inputs are forbidden. B intervals propagate the actual posterior-component mixture before taking the shared reference contrast.

The 81-node A/B pilot uses separate output directories/processes, at most two actual fits. Frozen V1 reuse is a read-only operation with a dedicated SHA audit. Run the cheap current-default control after a slot becomes available. A zero predictor has no intervals and is never treated as an uncertainty model.

`write_pilot_report.R --results DIR1,DIR2,... --report OUTPUT.md` requires matching model dataset/distance sets, reports all native-distance/family risk cells, and retains failed or unavailable predictions. `pooled_calibration.R` and `joint_power.R` are independently maintained validation modules; do not confuse ordinary selected LOO scores with nested OOF joint-power inputs.

`screen_development.R` records point/shape stage decisions and separates the common-reference offset gain from shape MSE. Its point/shape paired losses are deliberately marked `oof_verified=FALSE`; they are not the complete nested interval-calibration OOF vector. `archive_development.R` verifies the exact fitted source/DLL and scoring fingerprints, archives 18 unchanged development datasets, 90 node tables, V1 reuse audit and source snapshots. The 13 MB DLL stays under `/private/tmp`; its SHA256 and exact path are recorded. Repository archives are under `results/development81_comparison/archive/`.

New seed indices 9–14 (calibration) and 15+ (test) are reserved in separate domain namespaces. No generation-freeze record exists yet. In the present development-only phase, `prepare_data.R` refuses every non-development generation request even if supplied a provisional freeze. Power remains deferred until a candidate has credible paired OOF losses, all substantive gates are fixed and the test budget is certified before unsealing.
