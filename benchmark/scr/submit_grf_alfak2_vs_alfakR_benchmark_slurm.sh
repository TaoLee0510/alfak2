#!/usr/bin/env bash
#
# Fully expanded SLURM runner for the cross-repo GRF benchmark comparing alfak2
# against alfakR nn_prior modes. The workflow is intentionally split into three
# stages so every fit combination is one SLURM array task:
#
#   1. Prepare shared GRF/ABM inputs and task table. This stage also rebuilds
#      alfakR/alfak2 native code on the current Linux machine:
#        MODE=prepare sbatch benchmark/scr/submit_grf_alfak2_vs_alfakR_benchmark_slurm.sh
#
#   2. Submit one fit per array task after prepare finishes:
#        N_TASKS=$(($(wc -l < benchmark/results/grf_alfak2_vs_alfakR/tables/task_table.tsv) - 1))
#        sbatch --array=1-${N_TASKS}%64 benchmark/scr/submit_grf_alfak2_vs_alfakR_benchmark_slurm.sh
#
#      task_table.tsv is sorted with alfak2 first, then alfakR. To enforce that
#      order at scheduler level, submit the alfak2 index range first and submit
#      the remaining range after it finishes.
#
#   3. Summarize completed fit parts:
#        MODE=summarize sbatch benchmark/scr/submit_grf_alfak2_vs_alfakR_benchmark_slurm.sh
#
# Each fit task reads prepared input_rds/grf_rds files; it does not simulate ABM,
# rebuild inputs, or recompile native code at runtime.

#SBATCH --job-name=alfak2_grf_fit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=alfak2_grf_fit_%A_%a.out
#SBATCH --error=alfak2_grf_fit_%A_%a.err

set -euo pipefail

ALFAK2_REPO="${ALFAK2_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfak2}"
ALFAKR_REPO="${ALFAKR_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfakR}"
BASE_OUTPUT_DIR="${BASE_OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/grf_alfak2_vs_alfakR}"
SOURCE_INPUT_DIR="${SOURCE_INPUT_DIR:-}"
MODE="${MODE:-fit-task}"
MODULES="${MODULES:-R/4.4}"

METHODS="${METHODS:-none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step}"
MINOBS="${MINOBS:-5,10,20}"
N_SIM="${N_SIM:-12}"
LAMBDAS="${LAMBDAS:-0.2,0.4,0.8,1.6}"
TIME_GAPS="${TIME_GAPS:-2,4,8}"
TIME_STARTS="${TIME_STARTS:-0}"
PM="${PM:-5e-05}"
BETA_LEVELS="${BETA_LEVELS:-1e-05,5e-05,1e-04,1e-03,1e-02}"
NBOOT="${NBOOT:-45}"
GRID_N="${GRID_N:-81}"
SAMPLE_DEPTH="${SAMPLE_DEPTH:-2000}"
TIME_MAX="${TIME_MAX:-360}"
PASSAGE_INTERVAL="${PASSAGE_INTERVAL:-45}"
ALFAK2_INPUT_POLICIES="${ALFAK2_INPUT_POLICIES:-full,minobs_matched}"
ALFAK2_INPUT_DEPTH="${ALFAK2_INPUT_DEPTH:-effective}"
ALFAK2_OBSERVATION_MODEL="${ALFAK2_OBSERVATION_MODEL:-dirichlet_multinomial}"
ALFAK2_DM_CONCENTRATION="${ALFAK2_DM_CONCENTRATION:-50}"
ALFAK2_EFFECTIVE_DEPTH_MODE="${ALFAK2_EFFECTIVE_DEPTH_MODE:-min}"
ALFAK2_EFFECTIVE_DEPTH="${ALFAK2_EFFECTIVE_DEPTH:-}"
ALFAK2_LEGACY_WEIGHT="${ALFAK2_LEGACY_WEIGHT:-directly_informed}"
CORRECT_EFFLUX="${CORRECT_EFFLUX:-true}"
ALFAK2_LAMBDA_L_GRID="${ALFAK2_LAMBDA_L_GRID:-0.2}"
ALFAK2_LAMBDA_E_GRID="${ALFAK2_LAMBDA_E_GRID:-1}"
ALFAK2_SIGMA_OBS_GRID="${ALFAK2_SIGMA_OBS_GRID:-0.02}"
ALFAK2_GRAPH_EDGE_WEIGHT="${ALFAK2_GRAPH_EDGE_WEIGHT:-mutation}"
ALFAK2_ANCHOR_COUNT_REFERENCE="${ALFAK2_ANCHOR_COUNT_REFERENCE:-minobs}"
ALFAK2_ANCHOR_COUNT_POWER="${ALFAK2_ANCHOR_COUNT_POWER:-1}"
ALFAK2_LOCAL_SHELL_DEPTH="${ALFAK2_LOCAL_SHELL_DEPTH:-0}"
ALFAK2_GLOBAL_EXTRA_SHELL="${ALFAK2_GLOBAL_EXTRA_SHELL:-1}"
ALFAK2_MAX_NODES="${ALFAK2_MAX_NODES:-150000}"
FORCE_REFIT="${FORCE_REFIT:-false}"
FORCE_SIM="${FORCE_SIM:-false}"
REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE:-false}"
RECOMPILE_DLL="${RECOMPILE_DLL:-}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ -n "${MODULES:-}" ]] && type module >/dev/null 2>&1; then
  module purge
  for module_name in ${MODULES}; do
    module load "${module_name}"
  done
fi

mkdir -p "${BASE_OUTPUT_DIR}" "${BASE_OUTPUT_DIR}/slurm_logs"

cd "${ALFAK2_REPO}"

EXTRA_ARGS_ARRAY=()
if [[ -n "${EXTRA_ARGS:-}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_ARGS_ARRAY=(${EXTRA_ARGS})
fi

COMMON_ARGS=(
  "--alfak2-repo=${ALFAK2_REPO}"
  "--alfakR-repo=${ALFAKR_REPO}"
  "--output-dir=${BASE_OUTPUT_DIR}"
  "--n-cores=1"
  "--methods=${METHODS}"
  "--minobs=${MINOBS}"
  "--n-sim=${N_SIM}"
  "--lambdas=${LAMBDAS}"
  "--time-gaps=${TIME_GAPS}"
  "--time-starts=${TIME_STARTS}"
  "--pm=${PM}"
  "--beta-levels=${BETA_LEVELS}"
  "--nboot=${NBOOT}"
  "--grid-n=${GRID_N}"
  "--sample-depth=${SAMPLE_DEPTH}"
  "--time-max=${TIME_MAX}"
  "--passage-interval=${PASSAGE_INTERVAL}"
  "--correct-efflux=${CORRECT_EFFLUX}"
  "--alfak2-input-policies=${ALFAK2_INPUT_POLICIES}"
  "--alfak2-input-depth=${ALFAK2_INPUT_DEPTH}"
  "--alfak2-observation-model=${ALFAK2_OBSERVATION_MODEL}"
  "--alfak2-dm-concentration=${ALFAK2_DM_CONCENTRATION}"
  "--alfak2-effective-depth-mode=${ALFAK2_EFFECTIVE_DEPTH_MODE}"
  "--alfak2-legacy-weight=${ALFAK2_LEGACY_WEIGHT}"
  "--alfak2-lambda-l-grid=${ALFAK2_LAMBDA_L_GRID}"
  "--alfak2-lambda-e-grid=${ALFAK2_LAMBDA_E_GRID}"
  "--alfak2-sigma-obs-grid=${ALFAK2_SIGMA_OBS_GRID}"
  "--alfak2-graph-edge-weight=${ALFAK2_GRAPH_EDGE_WEIGHT}"
  "--alfak2-anchor-count-reference=${ALFAK2_ANCHOR_COUNT_REFERENCE}"
  "--alfak2-anchor-count-power=${ALFAK2_ANCHOR_COUNT_POWER}"
  "--alfak2-local-shell-depth=${ALFAK2_LOCAL_SHELL_DEPTH}"
  "--alfak2-global-extra-shell=${ALFAK2_GLOBAL_EXTRA_SHELL}"
  "--alfak2-max-nodes=${ALFAK2_MAX_NODES}"
  "--force-refit=${FORCE_REFIT}"
  "--force-sim=${FORCE_SIM}"
  "--reuse-dirty-cache=${REUSE_DIRTY_CACHE}"
)

if [[ -n "${SOURCE_INPUT_DIR}" ]]; then
  COMMON_ARGS+=("--source-input-dir=${SOURCE_INPUT_DIR}")
fi
if [[ -n "${ALFAK2_EFFECTIVE_DEPTH}" ]]; then
  COMMON_ARGS+=("--alfak2-effective-depth=${ALFAK2_EFFECTIVE_DEPTH}")
fi
if [[ -n "${RECOMPILE_DLL}" ]]; then
  COMMON_ARGS+=("--recompile-dll=${RECOMPILE_DLL}")
fi

echo "[$(date)] alfak2 vs alfakR GRF benchmark"
echo "  mode:              ${MODE}"
echo "  alfak2 repo:       ${ALFAK2_REPO}"
echo "  alfakR repo:       ${ALFAKR_REPO}"
echo "  output dir:        ${BASE_OUTPUT_DIR}"
if [[ -n "${SOURCE_INPUT_DIR}" ]]; then
  echo "  source input dir:  ${SOURCE_INPUT_DIR}"
fi
echo "  cpus per task:     ${SLURM_CPUS_PER_TASK:-1}"
if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
  echo "  mem per task:      ${SLURM_MEM_PER_NODE}M"
elif [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
  echo "  mem per task:      ${SLURM_MEM_PER_CPU}M per CPU"
else
  echo "  mem per task:      unknown"
fi
echo "  n_sim:             ${N_SIM}"
echo "  lambdas:           ${LAMBDAS}"
echo "  time gaps:         ${TIME_GAPS}"
echo "  minobs:            ${MINOBS}"
echo "  sim pm:            ${PM}"
echo "  fitted beta grid:  ${BETA_LEVELS}"
echo "  methods:           ${METHODS}"
echo "  alfak2 policies:   ${ALFAK2_INPUT_POLICIES}"
echo "  alfak2 depth:      ${ALFAK2_INPUT_DEPTH}"
echo "  alfak2 obs model:  ${ALFAK2_OBSERVATION_MODEL}"
echo "  alfak2 dm conc:    ${ALFAK2_DM_CONCENTRATION}"
if [[ -n "${ALFAK2_EFFECTIVE_DEPTH}" ]]; then
  echo "  alfak2 eff depth:  ${ALFAK2_EFFECTIVE_DEPTH}"
fi
echo "  alfak2 legacy wt:  ${ALFAK2_LEGACY_WEIGHT}"
echo "  correct efflux:    ${CORRECT_EFFLUX}"
echo "  alfak2 lambda_l:   ${ALFAK2_LAMBDA_L_GRID}"
echo "  alfak2 lambda_e:   ${ALFAK2_LAMBDA_E_GRID}"
echo "  alfak2 sigma_obs:  ${ALFAK2_SIGMA_OBS_GRID}"
echo "  alfak2 shell:      local=${ALFAK2_LOCAL_SHELL_DEPTH}, global=${ALFAK2_GLOBAL_EXTRA_SHELL}"
echo "  alfak2 max nodes:  ${ALFAK2_MAX_NODES}"
echo "  force refit:       ${FORCE_REFIT}"
echo "  force sim:         ${FORCE_SIM}"
if [[ -n "${RECOMPILE_DLL}" ]]; then
  echo "  recompile dll:     ${RECOMPILE_DLL}"
fi

case "${MODE}" in
  prepare)
    Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R \
      --mode=prepare \
      "${COMMON_ARGS[@]}" \
      "${EXTRA_ARGS_ARRAY[@]}"
    ;;

  fit-task|fit_task|fit)
    TASK_INDEX="${TASK_INDEX:-${SLURM_ARRAY_TASK_ID:-}}"
    if [[ -z "${TASK_INDEX}" ]]; then
      echo "TASK_INDEX or SLURM_ARRAY_TASK_ID is required for MODE=${MODE}." >&2
      exit 2
    fi
    echo "  task index:        ${TASK_INDEX}"
    Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R \
      --mode=fit-task \
      --task-index="${TASK_INDEX}" \
      "${COMMON_ARGS[@]}" \
      "${EXTRA_ARGS_ARRAY[@]}"
    ;;

  summarize)
    Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R \
      --mode=summarize \
      "${COMMON_ARGS[@]}" \
      "${EXTRA_ARGS_ARRAY[@]}"
    ;;

  *)
    echo "Unsupported MODE=${MODE}. Use prepare, fit-task, or summarize." >&2
    exit 2
    ;;
esac

echo "[$(date)] Done MODE=${MODE}"
