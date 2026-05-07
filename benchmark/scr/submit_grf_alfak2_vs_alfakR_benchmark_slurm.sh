#!/usr/bin/env bash
#
# SLURM array runner for the cross-repo GRF benchmark comparing alfak2 against
# alfakR nn_prior modes. Each array task runs one lambda x time-gap slice and
# writes an isolated slice directory. Merge slice tables with ordinary shell/R
# tooling after completion; this script intentionally does not modify package
# code or require either repo to be installed.
#
# Example:
#   sbatch --array=0-7%4 benchmark/scr/submit_grf_alfak2_vs_alfakR_benchmark_slurm.sh
#
# Optional overrides:
#   BASE_OUTPUT_DIR=/path/to/results METHOD_CORES=4 sbatch --array=0-7%4 ...

#SBATCH --job-name=alfak2_grf_cmp
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=alfak2_grf_cmp_%A_%a.out
#SBATCH --error=alfak2_grf_cmp_%A_%a.err

set -euo pipefail

ALFAK2_REPO="${ALFAK2_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfak2}"
ALFAKR_REPO="${ALFAKR_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfakR}"
BASE_OUTPUT_DIR="${BASE_OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/grf_alfak2_vs_alfakR}"

METHODS="${METHODS:-none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step}"
MINOBS="${MINOBS:-5,10,20}"
N_SIM="${N_SIM:-12}"
LAMBDAS="${LAMBDAS:-0.2,0.4,0.8,1.6}"
TIME_GAPS="${TIME_GAPS:-2,4,8}"
TIME_STARTS="${TIME_STARTS:-0}"
NBOOT="${NBOOT:-45}"
GRID_N="${GRID_N:-81}"
SAMPLE_DEPTH="${SAMPLE_DEPTH:-2000}"
TIME_MAX="${TIME_MAX:-360}"
PASSAGE_INTERVAL="${PASSAGE_INTERVAL:-45}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

if [[ -n "${MODULES:-}" ]] && type module >/dev/null 2>&1; then
  module purge
  for module_name in ${MODULES}; do
    module load "${module_name}"
  done
fi

split_csv() {
  local value="$1"
  local -n out_ref="$2"
  IFS=',' read -r -a out_ref <<< "${value}"
}

format_label() {
  local value="$1"
  value="${value//-/_m}"
  value="${value//./p}"
  echo "${value}"
}

split_csv "${LAMBDAS}" LAMBDA_VALUES
split_csv "${TIME_GAPS}" GAP_VALUES

if [[ "${#LAMBDA_VALUES[@]}" -eq 0 || "${#GAP_VALUES[@]}" -eq 0 ]]; then
  echo "LAMBDAS and TIME_GAPS must not be empty." >&2
  exit 2
fi

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
N_GAPS="${#GAP_VALUES[@]}"
N_TASKS=$(("${#LAMBDA_VALUES[@]}" * N_GAPS))

if (( TASK_ID < 0 || TASK_ID >= N_TASKS )); then
  echo "SLURM_ARRAY_TASK_ID=${TASK_ID} is outside 0..$((N_TASKS - 1)). Nothing to do."
  exit 0
fi

LAMBDA_IDX=$((TASK_ID / N_GAPS))
GAP_IDX=$((TASK_ID % N_GAPS))
LAMBDA="${LAMBDA_VALUES[$LAMBDA_IDX]}"
TIME_GAP="${GAP_VALUES[$GAP_IDX]}"
LAMBDA_LABEL="$(format_label "${LAMBDA}")"
GAP_LABEL="$(format_label "${TIME_GAP}")"

METHOD_CORES="${METHOD_CORES:-${SLURM_CPUS_PER_TASK:-4}}"
SLICE_OUTPUT_DIR="${BASE_OUTPUT_DIR}/slices/lambda_${LAMBDA_LABEL}/gap_${GAP_LABEL}"

mkdir -p "${SLICE_OUTPUT_DIR}" "${BASE_OUTPUT_DIR}/slurm_logs"

echo "[$(date)] Starting alfak2 vs alfakR GRF slice"
echo "  alfak2 repo:       ${ALFAK2_REPO}"
echo "  alfakR repo:       ${ALFAKR_REPO}"
echo "  base output:       ${BASE_OUTPUT_DIR}"
echo "  slice output:      ${SLICE_OUTPUT_DIR}"
echo "  task id:           ${TASK_ID}/${N_TASKS}"
echo "  lambda:            ${LAMBDA}"
echo "  time gap:          ${TIME_GAP}"
echo "  time starts:       ${TIME_STARTS}"
echo "  minobs:            ${MINOBS}"
echo "  methods:           ${METHODS}"
echo "  n_sim:             ${N_SIM}"
echo "  nboot:             ${NBOOT}"
echo "  sample depth:      ${SAMPLE_DEPTH}"
echo "  time max:          ${TIME_MAX}"
echo "  passage interval:  ${PASSAGE_INTERVAL}"
echo "  method cores:      ${METHOD_CORES}"

cd "${ALFAK2_REPO}"

EXTRA_ARGS_ARRAY=()
if [[ -n "${EXTRA_ARGS:-}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_ARGS_ARRAY=(${EXTRA_ARGS})
fi

Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R \
  --alfak2-repo="${ALFAK2_REPO}" \
  --alfakR-repo="${ALFAKR_REPO}" \
  --output-dir="${SLICE_OUTPUT_DIR}" \
  --n-cores="${METHOD_CORES}" \
  --methods="${METHODS}" \
  --minobs="${MINOBS}" \
  --n-sim="${N_SIM}" \
  --lambdas="${LAMBDA}" \
  --time-gaps="${TIME_GAP}" \
  --time-starts="${TIME_STARTS}" \
  --nboot="${NBOOT}" \
  --grid-n="${GRID_N}" \
  --sample-depth="${SAMPLE_DEPTH}" \
  --time-max="${TIME_MAX}" \
  --passage-interval="${PASSAGE_INTERVAL}" \
  "${EXTRA_ARGS_ARRAY[@]}"

touch "${SLICE_OUTPUT_DIR}/SLICE_DONE"
echo "[$(date)] Finished slice: lambda=${LAMBDA}, time_gap=${TIME_GAP}"
