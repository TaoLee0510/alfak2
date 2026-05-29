#!/usr/bin/env bash
#
# Fully expanded SLURM runner for alfak2-only GRF parameter calibration.
#
# This workflow is separate from the alfak2-vs-alfakR performance benchmark:
# calibration searches alfak2 working parameters on sparse GRF inputs; the
# selected parameters can then be passed to the comparison benchmark.
#
# Recommended 15-core workflow:
#
#   1. Prepare calibration task table from the already prepared comparison
#      benchmark inputs in SOURCE_INPUT_DIR:
#        MODE=prepare sbatch benchmark/scr/submit_grf_alfak2_parameter_calibration_slurm.sh
#
#   2. Submit one fit per array task, capped at 15 concurrent tasks:
#        N_TASKS=$(($(wc -l < benchmark/results/grf_alfak2_calibration/tables/task_table.tsv) - 1))
#        sbatch --array=1-${N_TASKS}%15 benchmark/scr/submit_grf_alfak2_parameter_calibration_slurm.sh
#
#   3. Summarize and rank parameter sets after the array finishes:
#        MODE=summarize sbatch benchmark/scr/submit_grf_alfak2_parameter_calibration_slurm.sh
#
# Each fit task reads reused input_rds/grf_rds files; it does not simulate ABM
# or rebuild inputs at runtime. If SOURCE_INPUT_DIR is empty or unavailable,
# the R calibration script can still generate its own inputs by omitting
# --source-input-dir, but this SLURM wrapper defaults to reuse.

#SBATCH --job-name=alfak2_cal_fit
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=alfak2_cal_fit_%A_%a.out
#SBATCH --error=alfak2_cal_fit_%A_%a.err

set -euo pipefail

ALFAK2_REPO="${ALFAK2_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfak2}"
ALFAKR_REPO="${ALFAKR_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfakR}"
BASE_OUTPUT_DIR="${BASE_OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/grf_alfak2_calibration}"
SOURCE_INPUT_DIR="${SOURCE_INPUT_DIR:-${ALFAK2_REPO}/benchmark/results/grf_alfak2_vs_alfakR}"
MODE="${MODE:-fit-task}"
MODULES="${MODULES:-R/4.4}"
FORCE_REFIT="${FORCE_REFIT:-false}"
FORCE_SIM="${FORCE_SIM:-false}"
REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE:-false}"
RECOMPILE_DLL="${RECOMPILE_DLL:-}"

# Sparse input grid. Defaults are the full calibration run:
# 162 parameter sets * 3 sims * 4 lambdas * 3 time gaps * 3 minobs = 17496 fits.
MINOBS="${MINOBS:-5,10,20}"
N_SIM="${N_SIM:-3}"
LAMBDAS="${LAMBDAS:-0.2,0.4,0.8,1.6}"
TIME_GAPS="${TIME_GAPS:-2,4,8}"
TIME_STARTS="${TIME_STARTS:-0}"
PM="${PM:-5e-05}"
SAMPLE_DEPTH="${SAMPLE_DEPTH:-2000}"
TIME_MAX="${TIME_MAX:-360}"
PASSAGE_INTERVAL="${PASSAGE_INTERVAL:-45}"
INPUT_POLICY="${INPUT_POLICY:-minobs_matched}"
ALFAK2_INPUT_DEPTH="${ALFAK2_INPUT_DEPTH:-effective}"
ALFAK2_OBSERVATION_MODEL="${ALFAK2_OBSERVATION_MODEL:-dirichlet_multinomial}"
ALFAK2_EFFECTIVE_DEPTH="${ALFAK2_EFFECTIVE_DEPTH:-}"
ALFAK2_MIN_CN="${ALFAK2_MIN_CN:-0}"
ALFAK2_MAX_CN="${ALFAK2_MAX_CN:-}"
ALFAK2_MAX_NODES="${ALFAK2_MAX_NODES:-150000}"

# Core alfak2 calibration grid. Default size:
# 3 legacy weights * 2 efflux settings * 3 lambda_l * 3 lambda_e * 3 sigma_obs
# * 1 dm concentration * 1 effective-depth mode * 1 shell setting = 162 params.
LEGACY_WEIGHTS="${LEGACY_WEIGHTS:-pi0,directly_informed,uniform}"
CORRECT_EFFLUX_VALUES="${CORRECT_EFFLUX_VALUES:-true,false}"
LAMBDA_L_VALUES="${LAMBDA_L_VALUES:-0.2,1,5}"
LAMBDA_E_VALUES="${LAMBDA_E_VALUES:-0.05,0.25,1}"
SIGMA_OBS_VALUES="${SIGMA_OBS_VALUES:-0.02,0.05,0.1}"
DM_CONCENTRATIONS="${DM_CONCENTRATIONS:-50}"
EFFECTIVE_DEPTH_MODES="${EFFECTIVE_DEPTH_MODES:-min}"
LOCAL_SHELL_DEPTHS="${LOCAL_SHELL_DEPTHS:-0}"
GLOBAL_EXTRA_SHELLS="${GLOBAL_EXTRA_SHELLS:-1}"

OBJECTIVE_SCOPE="${OBJECTIVE_SCOPE:-nn}"
GRAPH_EDGE_WEIGHTS="${GRAPH_EDGE_WEIGHTS:-mutation}"
ALFAK2_ANCHOR_COUNT_REFERENCE="${ALFAK2_ANCHOR_COUNT_REFERENCE:-minobs}"
ALFAK2_ANCHOR_COUNT_POWER="${ALFAK2_ANCHOR_COUNT_POWER:-1}"

OBJECTIVE_METRIC="${OBJECTIVE_METRIC:-sparse_composite}"
DIRECT_WEIGHT="${DIRECT_WEIGHT:-0.25}"
HOLDOUT_WEIGHT="${HOLDOUT_WEIGHT:-0.15}"
HOLDOUT_FRACTION="${HOLDOUT_FRACTION:-0.25}"
HOLDOUT_MIN_DIRECT="${HOLDOUT_MIN_DIRECT:-6}"
HOLDOUT_FAILURE_PENALTY="${HOLDOUT_FAILURE_PENALTY:-1}"
BIAS_WEIGHT="${BIAS_WEIGHT:-0.10}"
SPEARMAN_WEIGHT="${SPEARMAN_WEIGHT:-0.10}"
FALSE_HIGH_WEIGHT="${FALSE_HIGH_WEIGHT:-0.10}"

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
  "--source-input-dir=${SOURCE_INPUT_DIR}"
  "--n-cores=1"
  "--minobs=${MINOBS}"
  "--n-sim=${N_SIM}"
  "--lambdas=${LAMBDAS}"
  "--time-gaps=${TIME_GAPS}"
  "--time-starts=${TIME_STARTS}"
  "--pm=${PM}"
  "--sample-depth=${SAMPLE_DEPTH}"
  "--time-max=${TIME_MAX}"
  "--passage-interval=${PASSAGE_INTERVAL}"
  "--input-policy=${INPUT_POLICY}"
  "--alfak2-input-depth=${ALFAK2_INPUT_DEPTH}"
  "--alfak2-observation-model=${ALFAK2_OBSERVATION_MODEL}"
  "--alfak2-min-cn=${ALFAK2_MIN_CN}"
  "--alfak2-max-nodes=${ALFAK2_MAX_NODES}"
  "--legacy-weights=${LEGACY_WEIGHTS}"
  "--correct-efflux-values=${CORRECT_EFFLUX_VALUES}"
  "--lambda-l-values=${LAMBDA_L_VALUES}"
  "--lambda-e-values=${LAMBDA_E_VALUES}"
  "--sigma-obs-values=${SIGMA_OBS_VALUES}"
  "--graph-edge-weights=${GRAPH_EDGE_WEIGHTS}"
  "--alfak2-anchor-count-reference=${ALFAK2_ANCHOR_COUNT_REFERENCE}"
  "--alfak2-anchor-count-power=${ALFAK2_ANCHOR_COUNT_POWER}"
  "--dm-concentrations=${DM_CONCENTRATIONS}"
  "--effective-depth-modes=${EFFECTIVE_DEPTH_MODES}"
  "--local-shell-depths=${LOCAL_SHELL_DEPTHS}"
  "--global-extra-shells=${GLOBAL_EXTRA_SHELLS}"
  "--objective-scope=${OBJECTIVE_SCOPE}"
  "--objective-metric=${OBJECTIVE_METRIC}"
  "--direct-weight=${DIRECT_WEIGHT}"
  "--holdout-weight=${HOLDOUT_WEIGHT}"
  "--holdout-fraction=${HOLDOUT_FRACTION}"
  "--holdout-min-direct=${HOLDOUT_MIN_DIRECT}"
  "--holdout-failure-penalty=${HOLDOUT_FAILURE_PENALTY}"
  "--bias-weight=${BIAS_WEIGHT}"
  "--spearman-weight=${SPEARMAN_WEIGHT}"
  "--false-high-weight=${FALSE_HIGH_WEIGHT}"
  "--force-refit=${FORCE_REFIT}"
  "--force-sim=${FORCE_SIM}"
  "--reuse-dirty-cache=${REUSE_DIRTY_CACHE}"
)

if [[ -n "${ALFAK2_EFFECTIVE_DEPTH}" ]]; then
  COMMON_ARGS+=("--alfak2-effective-depth=${ALFAK2_EFFECTIVE_DEPTH}")
fi
if [[ -n "${ALFAK2_MAX_CN}" ]]; then
  COMMON_ARGS+=("--alfak2-max-cn=${ALFAK2_MAX_CN}")
fi
if [[ -n "${RECOMPILE_DLL}" ]]; then
  COMMON_ARGS+=("--recompile-dll=${RECOMPILE_DLL}")
fi

echo "[$(date)] alfak2 GRF parameter calibration"
echo "  mode:              ${MODE}"
echo "  alfak2 repo:       ${ALFAK2_REPO}"
echo "  alfakR repo:       ${ALFAKR_REPO}"
echo "  output dir:        ${BASE_OUTPUT_DIR}"
echo "  source input dir:  ${SOURCE_INPUT_DIR}"
echo "  cpus per task:     ${SLURM_CPUS_PER_TASK:-1}"
echo "  mem per task:      8G"
echo "  n_sim:             ${N_SIM}"
echo "  lambdas:           ${LAMBDAS}"
echo "  time gaps:         ${TIME_GAPS}"
echo "  minobs:            ${MINOBS}"
echo "  fitted pm:         ${PM}"
echo "  input policy:      ${INPUT_POLICY}"
echo "  alfak2 depth:      ${ALFAK2_INPUT_DEPTH}"
echo "  alfak2 obs model:  ${ALFAK2_OBSERVATION_MODEL}"
if [[ -n "${ALFAK2_EFFECTIVE_DEPTH}" ]]; then
  echo "  alfak2 eff depth:  ${ALFAK2_EFFECTIVE_DEPTH}"
fi
echo "  alfak2 max nodes:  ${ALFAK2_MAX_NODES}"
echo "  objective:         ${OBJECTIVE_SCOPE}/${OBJECTIVE_METRIC}"
echo "  legacy weights:    ${LEGACY_WEIGHTS}"
echo "  lambda_l values:   ${LAMBDA_L_VALUES}"
echo "  lambda_e values:   ${LAMBDA_E_VALUES}"
echo "  sigma_obs values:  ${SIGMA_OBS_VALUES}"
echo "  force refit:       ${FORCE_REFIT}"
if [[ -n "${RECOMPILE_DLL}" ]]; then
  echo "  recompile dll:     ${RECOMPILE_DLL}"
fi

case "${MODE}" in
  prepare)
    Rscript benchmark/scr/run_grf_alfak2_parameter_calibration.R \
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
    Rscript benchmark/scr/run_grf_alfak2_parameter_calibration.R \
      --mode=fit-task \
      --task-index="${TASK_INDEX}" \
      "${COMMON_ARGS[@]}" \
      "${EXTRA_ARGS_ARRAY[@]}"
    ;;

  summarize)
    Rscript benchmark/scr/run_grf_alfak2_parameter_calibration.R \
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
