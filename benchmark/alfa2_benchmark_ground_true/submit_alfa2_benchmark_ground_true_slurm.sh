#!/usr/bin/env bash

set -euo pipefail

SCRIPT_SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

write_export() {
  local file="$1"
  local name="$2"
  local value="${3-}"
  printf 'export %s=%q\n' "${name}" "${value}" >> "${file}"
}

load_requested_modules() {
  if [[ -n "${MODULES:-}" ]] && type module >/dev/null 2>&1; then
    module purge
    for module_name in ${MODULES}; do
      module load "${module_name}"
    done
  fi
}

run_task_stage() {
  local env_file="${1:-}"
  if [[ -z "${env_file}" || ! -f "${env_file}" ]]; then
    echo "Usage: ${SCRIPT_SELF} __run_task ENV_FILE" >&2
    exit 2
  fi

  # shellcheck disable=SC1090
  source "${env_file}"
  load_requested_modules

  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1

  local task_index="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
  cd "${ALFAK2_REPO}"

  echo "[$(date)] alfa2_benchmark_ground_true fit-task"
  echo "  task index:  ${task_index}"
  echo "  job id:      ${SLURM_JOB_ID:-NA}"
  echo "  array job:   ${SLURM_ARRAY_JOB_ID:-NA}"
  echo "  repo:        ${ALFAK2_REPO}"
  echo "  alfakR repo: ${ALFAKR_REPO}"
  echo "  output dir:  ${OUTPUT_DIR}"

  "${R_BIN}" "${RUNNER}" \
    "--mode=fit-task" \
    "--task-index=${task_index}" \
    "--output-dir=${OUTPUT_DIR}" \
    "--alfakR-repo=${ALFAKR_REPO}" \
    "--sample-depths=${SAMPLE_DEPTHS}" \
    "--wavelengths=${WAVELENGTHS}" \
    "--ground-truth-reps=${GROUND_TRUTH_REPS}" \
    "--fit-repeats=${FIT_REPEATS}" \
    "--soft-minobs=${SOFT_MINOBS}" \
    "--ntp=${NTP}" \
    "--nboot=${NBOOT}" \
    "--pmis=${PMIS}" \
    "--n0=${N0}" \
    "--nb=${NB}" \
    "--alfak2-local-shell-depth=${ALFAK2_LOCAL_SHELL_DEPTH}" \
    "--alfak2-global-extra-shell=${ALFAK2_GLOBAL_EXTRA_SHELL}" \
    "--alfak2-max-nodes=${ALFAK2_MAX_NODES}" \
    "--alfak2-eval-max=${ALFAK2_EVAL_MAX}" \
    "--alfak2-iter-max=${ALFAK2_ITER_MAX}" \
    "--alfak2-retry-eval-max=${ALFAK2_RETRY_EVAL_MAX}" \
    "--alfak2-retry-iter-max=${ALFAK2_RETRY_ITER_MAX}" \
    "--alfak2-lambda-l-grid=${ALFAK2_LAMBDA_L_GRID}" \
    "--alfak2-lambda-e-grid=${ALFAK2_LAMBDA_E_GRID}" \
    "--alfak2-sigma-obs-grid=${ALFAK2_SIGMA_OBS_GRID}" \
    "--force=${FORCE}"
}

if [[ "${1:-}" == "__run_task" ]]; then
  shift
  run_task_stage "$@"
  exit 0
fi

ALFAK2_REPO="${ALFAK2_REPO:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
ALFAKR_REPO="${ALFAKR_REPO:-$(cd "${ALFAK2_REPO}/.." && pwd)/alfakR}"
OUTPUT_DIR="${OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/alfa2_benchmark_ground_true}"
RUNNER="${RUNNER:-${ALFAK2_REPO}/benchmark/alfa2_benchmark_ground_true/run_alfa2_benchmark_ground_true.R}"
R_BIN="${R_BIN:-Rscript}"
MODULES="${MODULES:-R/4.4}"

SAMPLE_DEPTHS="${SAMPLE_DEPTHS:-1000,200}"
WAVELENGTHS="${WAVELENGTHS:-0.2,0.4,0.8,1.6}"
GROUND_TRUTH_REPS="${GROUND_TRUTH_REPS:-1:5}"
FIT_REPEATS="${FIT_REPEATS:-1:5}"
SOFT_MINOBS="${SOFT_MINOBS:-5,10,20}"
NTP="${NTP:-2}"
NBOOT="${NBOOT:-45}"
PMIS="${PMIS:-5e-05}"
N0="${N0:-2e5}"
NB="${NB:-2e7}"

ALFAK2_LOCAL_SHELL_DEPTH="${ALFAK2_LOCAL_SHELL_DEPTH:-0}"
ALFAK2_GLOBAL_EXTRA_SHELL="${ALFAK2_GLOBAL_EXTRA_SHELL:-2}"
ALFAK2_MAX_NODES="${ALFAK2_MAX_NODES:-150000}"
ALFAK2_EVAL_MAX="${ALFAK2_EVAL_MAX:-500}"
ALFAK2_ITER_MAX="${ALFAK2_ITER_MAX:-500}"
ALFAK2_RETRY_EVAL_MAX="${ALFAK2_RETRY_EVAL_MAX:-2000}"
ALFAK2_RETRY_ITER_MAX="${ALFAK2_RETRY_ITER_MAX:-2000}"
ALFAK2_LAMBDA_L_GRID="${ALFAK2_LAMBDA_L_GRID:-0.2}"
ALFAK2_LAMBDA_E_GRID="${ALFAK2_LAMBDA_E_GRID:-0.01}"
ALFAK2_SIGMA_OBS_GRID="${ALFAK2_SIGMA_OBS_GRID:-0.05}"
FORCE="${FORCE:-false}"

JOB_NAME="${JOB_NAME:-alfa2_gt_task}"
MEM="${MEM:-256G}"
CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
TIME_LIMIT="${TIME_LIMIT:-7-00:00:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"
DRY_RUN="${DRY_RUN:-false}"

mkdir -p "${OUTPUT_DIR}/slurm" "${OUTPUT_DIR}/slurm_logs"

load_requested_modules

cd "${ALFAK2_REPO}"

"${R_BIN}" "${RUNNER}" \
  "--mode=prepare" \
  "--output-dir=${OUTPUT_DIR}" \
  "--alfakR-repo=${ALFAKR_REPO}" \
  "--sample-depths=${SAMPLE_DEPTHS}" \
  "--wavelengths=${WAVELENGTHS}" \
  "--ground-truth-reps=${GROUND_TRUTH_REPS}" \
  "--fit-repeats=${FIT_REPEATS}" \
  "--soft-minobs=${SOFT_MINOBS}" \
  "--ntp=${NTP}" \
  "--nboot=${NBOOT}" \
  "--pmis=${PMIS}" \
  "--n0=${N0}" \
  "--nb=${NB}" \
  "--alfak2-local-shell-depth=${ALFAK2_LOCAL_SHELL_DEPTH}" \
  "--alfak2-global-extra-shell=${ALFAK2_GLOBAL_EXTRA_SHELL}" \
  "--alfak2-max-nodes=${ALFAK2_MAX_NODES}" \
  "--alfak2-eval-max=${ALFAK2_EVAL_MAX}" \
  "--alfak2-iter-max=${ALFAK2_ITER_MAX}" \
  "--alfak2-retry-eval-max=${ALFAK2_RETRY_EVAL_MAX}" \
  "--alfak2-retry-iter-max=${ALFAK2_RETRY_ITER_MAX}" \
  "--alfak2-lambda-l-grid=${ALFAK2_LAMBDA_L_GRID}" \
  "--alfak2-lambda-e-grid=${ALFAK2_LAMBDA_E_GRID}" \
  "--alfak2-sigma-obs-grid=${ALFAK2_SIGMA_OBS_GRID}"

RUN_INDEX="${OUTPUT_DIR}/tables/run_index.tsv"
if [[ ! -f "${RUN_INDEX}" ]]; then
  echo "Missing run index after prepare: ${RUN_INDEX}" >&2
  exit 2
fi

N_TASKS=$(($(wc -l < "${RUN_INDEX}") - 1))
if (( N_TASKS < 1 )); then
  echo "run_index.tsv contains no tasks: ${RUN_INDEX}" >&2
  exit 2
fi

ENV_FILE="${OUTPUT_DIR}/slurm/alfa2_benchmark_ground_true.env"
: > "${ENV_FILE}"
write_export "${ENV_FILE}" "ALFAK2_REPO" "${ALFAK2_REPO}"
write_export "${ENV_FILE}" "ALFAKR_REPO" "${ALFAKR_REPO}"
write_export "${ENV_FILE}" "OUTPUT_DIR" "${OUTPUT_DIR}"
write_export "${ENV_FILE}" "RUNNER" "${RUNNER}"
write_export "${ENV_FILE}" "R_BIN" "${R_BIN}"
write_export "${ENV_FILE}" "MODULES" "${MODULES}"
write_export "${ENV_FILE}" "SAMPLE_DEPTHS" "${SAMPLE_DEPTHS}"
write_export "${ENV_FILE}" "WAVELENGTHS" "${WAVELENGTHS}"
write_export "${ENV_FILE}" "GROUND_TRUTH_REPS" "${GROUND_TRUTH_REPS}"
write_export "${ENV_FILE}" "FIT_REPEATS" "${FIT_REPEATS}"
write_export "${ENV_FILE}" "SOFT_MINOBS" "${SOFT_MINOBS}"
write_export "${ENV_FILE}" "NTP" "${NTP}"
write_export "${ENV_FILE}" "NBOOT" "${NBOOT}"
write_export "${ENV_FILE}" "PMIS" "${PMIS}"
write_export "${ENV_FILE}" "N0" "${N0}"
write_export "${ENV_FILE}" "NB" "${NB}"
write_export "${ENV_FILE}" "ALFAK2_LOCAL_SHELL_DEPTH" "${ALFAK2_LOCAL_SHELL_DEPTH}"
write_export "${ENV_FILE}" "ALFAK2_GLOBAL_EXTRA_SHELL" "${ALFAK2_GLOBAL_EXTRA_SHELL}"
write_export "${ENV_FILE}" "ALFAK2_MAX_NODES" "${ALFAK2_MAX_NODES}"
write_export "${ENV_FILE}" "ALFAK2_EVAL_MAX" "${ALFAK2_EVAL_MAX}"
write_export "${ENV_FILE}" "ALFAK2_ITER_MAX" "${ALFAK2_ITER_MAX}"
write_export "${ENV_FILE}" "ALFAK2_RETRY_EVAL_MAX" "${ALFAK2_RETRY_EVAL_MAX}"
write_export "${ENV_FILE}" "ALFAK2_RETRY_ITER_MAX" "${ALFAK2_RETRY_ITER_MAX}"
write_export "${ENV_FILE}" "ALFAK2_LAMBDA_L_GRID" "${ALFAK2_LAMBDA_L_GRID}"
write_export "${ENV_FILE}" "ALFAK2_LAMBDA_E_GRID" "${ALFAK2_LAMBDA_E_GRID}"
write_export "${ENV_FILE}" "ALFAK2_SIGMA_OBS_GRID" "${ALFAK2_SIGMA_OBS_GRID}"
write_export "${ENV_FILE}" "FORCE" "${FORCE}"

ARRAY_SPEC="1-${N_TASKS}"
if [[ -n "${ARRAY_LIMIT}" ]]; then
  ARRAY_SPEC="${ARRAY_SPEC}%${ARRAY_LIMIT}"
fi

echo "[$(date)] submitting alfa2_benchmark_ground_true"
echo "  tasks:       ${N_TASKS}"
echo "  array:       ${ARRAY_SPEC}"
echo "  cpu/task:    ${CPUS_PER_TASK}"
echo "  mem/task:    ${MEM}"
echo "  time/task:   ${TIME_LIMIT}"
echo "  output dir:  ${OUTPUT_DIR}"
echo "  env file:    ${ENV_FILE}"

SBATCH_ARGS=(
  "--parsable"
  "--array=${ARRAY_SPEC}"
  "--cpus-per-task=${CPUS_PER_TASK}"
  "--mem=${MEM}"
  "--time=${TIME_LIMIT}"
  "--job-name=${JOB_NAME}"
  "--output=${OUTPUT_DIR}/slurm_logs/${JOB_NAME}_%A_%a.out"
  "--error=${OUTPUT_DIR}/slurm_logs/${JOB_NAME}_%A_%a.err"
  "--chdir=${ALFAK2_REPO}"
  "${SCRIPT_SELF}"
  "__run_task"
  "${ENV_FILE}"
)

if [[ "${DRY_RUN}" == "true" ]]; then
  printf 'DRY-RUN: sbatch'
  printf ' %q' "${SBATCH_ARGS[@]}"
  printf '\n'
else
  JOB_ID=$(sbatch "${SBATCH_ARGS[@]}" | awk -F';' '{print $1}')
  printf '%s\n' "${JOB_ID}" > "${OUTPUT_DIR}/slurm/fit_task_array.jobid"
  echo "Submitted job array: ${JOB_ID}"
fi
