# Numerical recovery and resource accounting, before accuracy inspection

The original statistical thresholds and primary groups remain unchanged.
Analysis rules were written at 21:02:36 UTC, after the frozen fit batch started
at 20:58:26 UTC but before comparative accuracy was inspected. Rules and initial
postprocessing sources were hash-frozen at 21:24:40 UTC. This is an internal
development analysis plan, not an external preregistration. Only numerical
status, runtime, source/data construction and resource results were reviewed
before these rules and the present numerical amendments.

External hard deadlines now control the remaining fit subprocesses. Analysis
uses two separate axes: numerical prediction availability, and compliance with
the original computing budget. Original failures are never overwritten.

After all sixty planned cases are attempted, `prepare_analysis.R` may recover
a posterior from an already complete, matching source/data cache if the child
was cut during output/scoring, or through independently verified joint-mode /
Hessian evidence reconstruction. It performs no new fit. Recovered cases retain
their raw status, runtime, resource flag and recovery label. Missing nodes or
failed refinement remain unavailable, with zero prediction/interval availability
and NA risks; they cannot be deleted or averaged away.

A complete numerical comparison is produced only if every case has a trusted
posterior. Otherwise the full five-design ranking is withheld and all failures
and availability counts are reported. Even a complete off-budget posterior
cannot qualify a design for a larger-domain stage. Per-case budget failures
remain a separate veto, without changing any point, shape, edge, interval-score
or family-safety threshold. Neither reconstruction nor this amendment opens
reserved calibration/test seeds or changes the production API.
