# alfa2_benchmark_ground_true

This benchmark reproduces the ALFA-K original ABM ground-truth generation process and runs matching `alfak2` / `alfakR` method grids on the generated truth.

## Ground Truth

The generator follows `/Users/4482173/Downloads/ALFA-K_orignal/scripts/S01_run_abm_sims.R`:

- founder karyotype: 22 chromosomes, all copy number 2
- `Nwaves = 10`
- `gen_randscape()` logic unchanged
- `times = c(0, 300)`
- `pmis = 5e-05`
- `run_abm_simulation_grf()` settings unchanged:
  - `abm_pop_size = 5e4`
  - `abm_max_pop = 2e6`
  - `abm_delta_t = 0.1`
  - `abm_culling_survival = 0.01`
  - `abm_record_interval = -1`
  - `abm_seed = 42`
  - `normalize_freq = FALSE`
- `resample_sim()` logic unchanged

The extension is `sample_depth = 1000, 200`. The original wavelengths `0.2, 0.4, 0.8, 1.6` are retained, with 5 ground-truth repeats per depth and wavelength.

## Method Grid

Per `sample_depth`, `alfak2` has 4 input settings times 9 d1/d2 extrapolation methods:

- `full`
- `soft_minobs = 5`
- `soft_minobs = 10`
- `soft_minobs = 20`

Per `sample_depth`, `alfakR` has 3 `minobs` settings times 4 `NN_prior` settings:

- `minobs = 5, 10, 20`
- `NN_prior = None, empirical, empirical_censored, empirical_censored_weighted`

Each method parameter pair is run 5 fit repeats for every ground truth.

## Commands

Prepare indexes:

```sh
Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=prepare
```

Generate one ground truth:

```sh
Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=ground-truth --ground-truth-index=1
```

Run one task:

```sh
Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=fit-task --task-index=1
```

Summarize completed tasks:

```sh
Rscript benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R --mode=summarize
```

Submit all run tasks to Slurm as one array job:

```sh
bash benchmark/alfa2_benchmark_ground_true/submit_alfa2_benchmark_ground_true_slurm.sh
```

The Slurm submitter assigns each run task 1 CPU, 256G memory, and 7 days. It writes logs to `benchmark/results/alfa2_benchmark_ground_true/slurm_logs`.

On HPC, the submitter uses:

- `ALFAK2_REPO=/share/lab_crd/lab_crd/taoli/Project/alfak2`
- `ALFAKR_REPO=/share/lab_crd/lab_crd/taoli/Project/alfakR`

The R runner requires `pkgload` and loads both packages with `pkgload::load_all()` from those source repositories. It validates the loaded namespace paths and fails instead of falling back to any `alfakR` installed in the active R module library.

Outputs default to `benchmark/results/alfa2_benchmark_ground_true`.
