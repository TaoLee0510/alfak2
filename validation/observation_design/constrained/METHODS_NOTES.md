# Interpretation and supplemental checks

The Fisher constraint is a local variance constraint, not a bound on realized
fitness error or a frequentist noninferiority result. It uses a delta-method
birth covariance around each fixed parameter point, and the original fitness
and initial-composition prior curvatures. Bias, nonlinear posterior effects and
finite-population process noise are assessed only through empirical validation.

The frozen search constrained the seven nominal-A seeded states. A subsequent
truth-blind check, before precision results were inspected, also evaluates the
union of the original control300 actual observed states (15 states) and every
one of the 27 states. Both selected candidates satisfy the same variance limit
at all of these states across all nine points and both probability floors.
Worst all-state variance ratios are 0.9984983 (single) and 0.9849186 (dual).
This check changes neither the chosen candidates nor their parameters. Its
reproducible calculation is `check_support_and_data.R`, with results in
`../results/constrained_v2/actual_observed_information_diagnostic.csv`.

The fixed-reference variance constraint does not imply a reduction of globally
centered Shape variance, which is a different contrast using off-diagonal
covariances. In the saved selected-design planning diagnostics, the single
candidate's worst nominal-d0 mean Shape variance ratio is 1.025878; the dual
candidate's is 0.992931. Omitting the eta prior gives 1.024940 and 0.993465,
respectively. These are planning values, not observed error changes. The
unchanged empirical Shape gate remains part of the final conjunction. Full
ratios are in `../results/constrained_v2/planning_sensitivity_ratios.csv`.

Primary accuracy comparisons still use each biological unit's unchanged
original-control actual distances. Recomputing distances after broader seeding
would describe a different question, so native and planned distances are
auxiliary views only.

In the process-noise table, each stratum's numerator and denominator are means
of individual karyotype-frequency marginal variances, followed by a ratio. They
are not the variance of the sum of all counts within that stratum. The total
conditional readout covariance is computed as

    mean(diag(q) - q q') / 300 + sample_covariance(q).

The first term describes multinomial sampling conditional on a random process
endpoint. The second describes endpoint frequency variation from the branching
process, with initial integer population fixed. The deterministic comparator
normalizes the independent exact mean population action. Its mean fraction
need not equal the mean of stochastic endpoint fractions. Both are saved in
the per-node diagnostics. No derived frequency-variance ratio is used as a
fitness-posterior calibration factor.

The bootstrap resamples the 100 process trajectories at each fixed point, not
the six development landscapes. It describes Monte Carlo sensitivity of this
diagnostic; especially for rare distant states it does not provide a rigorous
tail guarantee. No additional posterior fits or reserved datasets are used.
