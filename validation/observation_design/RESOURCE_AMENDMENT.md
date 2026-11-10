# Numerical resource amendment, before primary design fits

The first composition600 smoke fit converged at quadrature order 17 in 118.32s,
but `/usr/bin/time -l` recorded 7,299,710,976 bytes maximum RSS (6.798 GiB).
It therefore failed the frozen 6-GiB resource gate. No primary batch was started.
Exact source v2 and the smoke output are retained.

The inference target, data, priors, integration rule and evaluation criteria are
unchanged. The standalone forward solver now uses an eighth-order positive
Taylor action of A-min(diag(A))*I, tracking the factored common growth exactly.
The diagonal shift is a dynamic AD conditional, not a comparison fixed at tape
construction. This aims to reduce the high-step nested-AD tapes required by
second-order splitting at extreme quadrature scales. Production code is untouched.

Every fitted component still compares step halving and a separately coded
converged continuous-time mean solution to 1e-5, with backend agreement to 1e-8.
Unit comparisons against the original solver use its finer .025 step; direct
pointwise likelihood checks use an independent continuous-time calculation.
New and old complete smoke means/covariance are compared before accepting the
optimization. The same single-fit data must pass the original 6-GiB/300s gate
before at most two primary workers are launched. Solver changes are numerical
repairs, not additional fitted biological candidates or changes selected from
design accuracy results.

Conditional fits explicitly release their primary AD tapes and collect
unreferenced `sdreport` AD-report tapes before the next component. Only numeric
posterior summaries are retained. This avoids accumulating external native
allocations between quadrature nodes and does not alter the inference target.
