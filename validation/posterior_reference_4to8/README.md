# Posterior reference, 4–8 nodes

This fixed-data checkpoint isolates numerical posterior approximation error.
The model/data cases and precision criteria are in [PROTOCOL.md](PROTOCOL.md).
The completed numerical conclusions are written to [REVIEW.md](REVIEW.md).
It consumes no reserved accuracy calibration/test seed.

From the repository root, a fresh reproduction uses installed R packages
`Rcpp`, `TMB`, `Matrix`, `posterior`, and `testthat`, and an oracle-capable
alfak2 DLL built from this repository. Point `ALFAK2_REFERENCE_DLL` to that
installed DLL. Choose fresh output and chain directories to preserve the
recorded results:

```sh
export ALFAK2_REFERENCE_DLL="/path/to/library/alfak2/libs/alfak2.so"
export ALFAK2_REFERENCE_OUTPUT="/private/tmp/alfak2-reference-reproduction/results"
export ALFAK2_REFERENCE_RAW="/private/tmp/alfak2-reference-reproduction/chains"
export ALFAK2_REFERENCE_REPORT="/private/tmp/alfak2-reference-reproduction/REVIEW.md"
Rscript validation/posterior_reference_4to8/validate.R
Rscript validation/posterior_reference_4to8/run.R
Rscript validation/posterior_reference_4to8/analyze.R
Rscript validation/posterior_reference_4to8/report.R
```

The previous two-node quadrature input is archived in
`inputs/two_node_reference.rds`; an optional `ALFAK2_REFERENCE_TWONODE` override
is available. No input download or old temporary reference file is required.
The code refuses changed case/source fingerprints or mixed raw-chain caches.

The original primary run's exact pre-execution code is archived under
`results/computation_source/`, with all recorded source hashes verified before
archiving. Portable file-path overrides and stricter future cache guards were
added while the original sampling loop was running; they did not change that
loop, its loaded sampler/target functions, case data, seeds, or its results.
The original fresh raw-chain directory had no pre-existing cache. Use fresh
directories for reproducing with the portable scripts; their code fingerprints
intentionally differ from the original archived runner.

The live original process emitted a trailing parse error after all 16 outputs
were written because those I/O edits changed its file position. Every saved
case/chain/hash was independently checked afterward; see
`results/completion_integrity.txt`. The error did not change the parsed loop or
loaded inference functions. A separate portable-path rerun passes all 99
validation assertions. The error and its verification remain in the record.

Raw chains stay outside Git; `case_status.csv` records their exact paths and
SHA256. Compact posterior summaries, diagnostics, all comparison tables,
fixed cases, source archives and independent validation results are retained
under `results/`. Any reference failing its preregistered precision cap remains
explicitly uncertified rather than disappearing from the results.
