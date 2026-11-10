# Fixed-budget observation-design development study

This is a standalone validation prototype. It does not change production
functions or the package default. The initial design, disclosed amendments and
decision rules are in `PROTOCOL.md`, `DESIGN_AMENDMENT.md`,
`RESOURCE_AMENDMENT.md` and `ANALYSIS_RULES.md`.

The completed checkpoint is `REVIEW.md`: all 60 original cases were attempted,
54 are numerically available and 6 remain unavailable. The five-design ranking
is withheld. `REVIEW_v4.md` preserves the original batch report. The latest
numerical prototype is frozen in `results/direct_source_v7` with forward helpers
in `results/direct_forward_source_v2`: 211 checks pass, but the worst-case smoke
hits the 300-second hard deadline at a 1.15-GiB peak. No full direct-method batch
was launched. `NEXT_NUMERICAL_PLAN.md` is a proposed follow-up, not completed work.

The frozen statistical computation is `results/computation_source_v4/`, with
60 cases partitioned by complete biological unit between two workers. Its
internal R timeout proved insufficient because TMB can catch the interrupt.
The first workers were stopped after 31 completed cases, then the remaining
cases continued under `results/control_source_v1/bounded_run.py --part 0/1`,
with `ALFAK2_REPO` set to the absolute repository path. Each controller owns
only one fitting R child and imposes a cumulative external deadline. Do not
launch extra workers or blindly restart completed controllers. Preserve the
existing outputs when reproducing in a fresh workspace.

The workers need the exact native alfak2 diagnostic DLL recorded in their
`source_sha256.csv`, the original prior/reference helpers and the TMB toolchain.
The standalone observation DLL is compiled into a source-hash-specific temporary
directory. Complete and conditional-node caches are keyed by source and data
hashes; an unrelated cache must never be reused. Frozen data generation and
immutable source directories are evidence, not disposable run scratch.

After both bounded controllers have completed, run these sequentially from the
repository root (substitute the absolute path in the environment variable):

```sh
Rscript validation/observation_design/prepare_analysis.R
export ALFAK2_DESIGN_OUTPUT=/absolute/path/to/alfak2/validation/observation_design/results/analysis_input
Rscript validation/observation_design/analyze.R
Rscript validation/observation_design/audit.R
Rscript validation/observation_design/summarize.R
Rscript validation/observation_design/integrity.R
```

`prepare_analysis.R` performs no new fits. It keeps raw failures/resource flags,
uses only matching complete caches, and can independently reconstruct evidence
at saved joint modes. It cannot fabricate missing nodes. `analyze.R` refuses
incomplete case sets, failed fits, nonfinite primary risks
or unavailable predictions/intervals. It never averages away a failure.
If numerical availability remains incomplete, `analyze.R` writes the complete
status/availability tables and stops before ranking; skip `audit.R` in that
case and run `summarize.R` for the explicit incomplete-case report.
`audit.R` checks the saved posterior mixture moments, fixed-reference projection,
global error decomposition and input/source integrity. `summarize.R` writes
`REVIEW.md`, descriptive alternative distance views and their denominators.
`integrity.R` checks all data/consumed sources even when numerical availability
is incomplete. None of these scripts performs additional fits or opens reserved
test seeds. See `POSTPROCESSING_AMENDMENT.md` for the timeline of the internal
analysis freeze and the distinction between recovered numerics and resources.

After reconstructing the original archives, `python3
validation/observation_design/checkpoint_report.py` writes the final checkpoint
summary including the later failed numerical probes. It does not perform fits
or turn an incomplete comparison into an accuracy ranking. Preserve all failed
source versions and resource logs. The direct batch launchers are intentionally
blocked until a smoke test passes both numerical refinement and resources.

All primary comparisons use each paired unit's original control300 support
distance, keeping the same target states across layouts. Native-distance
summaries have different support and sometimes different sample denominators;
they cannot replace the primary comparison. The independent experimental units
are six landscapes, rather than the twelve biological observations.

The archived v1/v2 smoke failures and v3 resource repair remain visible. Do not
promote an all-distance claim from a favorable average or the Fisher planning
table. Follow the frozen conjunction of shape, point, edge, interval-score and
family-safety checks, and treat it as development evidence only.
