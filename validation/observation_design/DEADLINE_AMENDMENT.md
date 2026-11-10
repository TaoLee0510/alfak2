# External deadline control; statistical computation remains frozen v4

During the first batch, the R elapsed-time interrupt did not enforce a hard
300-second boundary. TMB can catch an interrupt during a marginal-objective
evaluation and return NaN, then successfully compute an sdreport. Recorded
cases took 305.006, 319.578 and 345.510 seconds. Their original numerical
status and budget failure remain visible. They cannot be certified as passes
under the original resource gate.

At 2026-11-04 21:32:33 UTC, only the two owned worker PIDs (19616/19617) were
suspended. There were 31 completed cases. The in-progress cases had used at
most 17.677 and 5.476 seconds since the previous result was written; these
conservative elapsed amounts are charged on resumption. The two old workers
were then terminated. No other task or HPC process was affected.

The external Python controller starts one fresh R subprocess per remaining
case, with a hard kill at the remaining cumulative budget. It runs the exact
frozen v4 model/driver, using the unchanged source/data-keyed conditional caches.
Two controllers mean at most two fitting R subprocesses. Old complete results
are replayed directly; an interrupted case gets only 300 minus its already-used
seconds. Restarting does not reset the budget. Child startup and loading time
are additionally included in the conservative external wall-clock allowance.
Per-child peak RSS is recorded through wait4 resource statistics. Outcomes,
including killed/incomplete cases, are retained in `fits_bounded_0/1`.

A separate normal-exit and hard-timeout test passed before launch: the timed
child was killed by SIGKILL in 0.208 seconds for a 0.2-second allowance. The
controller changes process supervision only, not model, data, posterior target,
reference, distance grouping, scoring or development thresholds.

Any unavailable case still prevents a complete statistical ranking. Offline
evidence reconstruction, if possible from complete cached orders, is labeled
separately and cannot erase a computing-budget failure. No larger-domain run
or reserved calibration/test is opened after a failed resource gate.
