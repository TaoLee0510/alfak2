# Completed local round

2026-11-05 01:50:52 EST

No model workers remain active. Current numerical source: results/importance_source_v16. Analysis source: importance_analysis_source_v2. Report/figure source: importance_report_source_v3.

107/107 validation checks passed; both high-precision8-node reference gates passed; both27-node preflights passed.48/48 final mechanism cases passed all numerical/resource gates. Maximum complete wall200.224s; peak5.516GiB. All dataset pairings, source/runtime integrity and the two user-original files passed.413 selected result/sample/report artifacts (7,468,842,816bytes) hashed in importance_v1/final_artifact_sha256.csv. Runtime binaries archived separately.

The all-distance generating-truth improvement gate FAILED: correctedA RMSE changes+0.180/+0.110/-0.090%; coverage95=87.15/82.68/89.31% at original-control distance0/1/2. Initial-composition oracle improves point RMSE more than replacing finite process by mean propagation, but fails to improve all coverage. Do not promote oracle as deployable, open reserved tests, change production defaults, submit HPC, commit or push under this completed round.

Final review: importance/REVIEW.md. Main outputs: importance_v1/mechanism (all_fit_status.csv, all_case_metrics.csv, risk/paired differences, mechanism_summary.pdf/png). Historical failures and frozen snapshots remain intact. Future model work requires a separately scoped phase; no such model changes were executed here.
