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
  if [[ -z "${MODULES:-}" ]]; then
    return 0
  fi

  if ! type module >/dev/null 2>&1 && ! type ml >/dev/null 2>&1; then
    for init_script in \
      /etc/profile.d/modules.sh \
      /etc/profile.d/lmod.sh \
      /usr/share/Modules/init/bash \
      /usr/share/lmod/lmod/init/bash; do
      if [[ -f "${init_script}" ]]; then
        # shellcheck disable=SC1090
        source "${init_script}"
        break
      fi
    done
  fi

  if type module >/dev/null 2>&1; then
    module purge
    for module_name in ${MODULES}; do
      module load "${module_name}"
    done
  elif type ml >/dev/null 2>&1; then
    ml purge
    for module_name in ${MODULES}; do
      ml "${module_name}"
    done
  else
    echo "Unable to load required module(s): ${MODULES}" >&2
    echo "Neither module nor ml is available in this shell." >&2
    exit 2
  fi

  if ! command -v "${R_BIN:-Rscript}" >/dev/null 2>&1; then
    echo "Rscript is unavailable after loading module(s): ${MODULES}" >&2
    exit 2
  fi
}

run_task_stage() {
  local env_file="${1:-}"
  if [[ -z "${env_file}" || ! -f "${env_file}" ]]; then
    echo "Usage: ${SCRIPT_SELF} __run_task ENV_FILE [TASK_ID_MAP]" >&2
    exit 2
  fi
  local task_id_map="${2:-}"

  # shellcheck disable=SC1090
  source "${env_file}"
  load_requested_modules

  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1

  local array_task_index="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
  local task_index="${array_task_index}"
  if [[ -n "${task_id_map}" ]]; then
    if [[ ! -f "${task_id_map}" ]]; then
      echo "Missing task id map: ${task_id_map}" >&2
      exit 2
    fi
    task_index="$(awk -v row="${array_task_index}" 'NR == row { print $1; found = 1; exit } END { if (!found) exit 3 }' "${task_id_map}")"
  fi

  cd "${ALFAK2_REPO}"

  echo "[$(date)] alfa2_benchmark_ground_true fit-task"
  echo "  array task:  ${array_task_index}"
  echo "  task index:  ${task_index}"
  echo "  job id:      ${SLURM_JOB_ID:-NA}"
  echo "  array job:   ${SLURM_ARRAY_JOB_ID:-NA}"
  echo "  task map:    ${task_id_map:-NA}"
  echo "  repo:        ${ALFAK2_REPO}"
  echo "  alfakR repo: ${ALFAKR_REPO}"
  echo "  output dir:  ${OUTPUT_DIR}"
  echo "  modules:     ${MODULES}"
  echo "  Rscript:     $(command -v "${R_BIN}")"

  "${R_BIN}" "${RUNNER}" \
    "--mode=fit-task" \
    "--task-index=${task_index}" \
    "--output-dir=${OUTPUT_DIR}" \
    "--alfak2-repo=${ALFAK2_REPO}" \
    "--alfakR-repo=${ALFAKR_REPO}" \
    "--sample-depths=${SAMPLE_DEPTHS}" \
    "--wavelengths=${WAVELENGTHS}" \
    "--ground-truth-times=${GROUND_TRUTH_TIMES}" \
    "--passage-times=${PASSAGE_TIMES}" \
    "--ground-truth-reps=${GROUND_TRUTH_REPS}" \
    "--fit-repeats=${FIT_REPEATS}" \
    "--soft-minobs=${SOFT_MINOBS}" \
    "--ntp=${NTP}" \
    "--alfakR-dt=${ALFAKR_DT}" \
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
MODULES="${MODULES-R/4.4}"

SAMPLE_DEPTHS="${SAMPLE_DEPTHS:-1000,200}"
WAVELENGTHS="${WAVELENGTHS:-0.2,0.4,0.8,1.6}"
GROUND_TRUTH_TIMES="${GROUND_TRUTH_TIMES:-0,180}"
PASSAGE_TIMES="${PASSAGE_TIMES:-0,180}"
GROUND_TRUTH_REPS="${GROUND_TRUTH_REPS:-1:5}"
FIT_REPEATS="${FIT_REPEATS:-1:5}"
SOFT_MINOBS="${SOFT_MINOBS:-5,10,20}"
NTP="${NTP:-2}"
ALFAKR_DT="${ALFAKR_DT:-1}"
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
CPUS_PER_TASK="${CPUS_PER_TASK:-1}"
TIME_LIMIT="${TIME_LIMIT:-7-00:00:00}"
ARRAY_LIMIT="${ARRAY_LIMIT:-}"
DRY_RUN="${DRY_RUN:-false}"

MEM_DEPTH200_ALFAK2_REGULAR="${MEM_DEPTH200_ALFAK2_REGULAR:-16G}"
MEM_DEPTH200_ALFAKR="${MEM_DEPTH200_ALFAKR:-8G}"
MEM_DEPTH200_ALFAK2_GGFULL="${MEM_DEPTH200_ALFAK2_GGFULL:-128G}"
MEM_DEPTH1000_ALFAK2_REGULAR="${MEM_DEPTH1000_ALFAK2_REGULAR:-32G}"
MEM_DEPTH1000_ALFAKR="${MEM_DEPTH1000_ALFAKR:-16G}"
MEM_DEPTH1000_ALFAK2_GGFULL="${MEM_DEPTH1000_ALFAK2_GGFULL:-256G}"

mkdir -p "${OUTPUT_DIR}/slurm" "${OUTPUT_DIR}/slurm_logs"

load_requested_modules

cd "${ALFAK2_REPO}"

"${R_BIN}" "${RUNNER}" \
  "--mode=prepare" \
  "--output-dir=${OUTPUT_DIR}" \
  "--alfak2-repo=${ALFAK2_REPO}" \
  "--alfakR-repo=${ALFAKR_REPO}" \
  "--sample-depths=${SAMPLE_DEPTHS}" \
  "--wavelengths=${WAVELENGTHS}" \
  "--ground-truth-times=${GROUND_TRUTH_TIMES}" \
  "--passage-times=${PASSAGE_TIMES}" \
  "--ground-truth-reps=${GROUND_TRUTH_REPS}" \
  "--fit-repeats=${FIT_REPEATS}" \
  "--soft-minobs=${SOFT_MINOBS}" \
  "--ntp=${NTP}" \
  "--alfakR-dt=${ALFAKR_DT}" \
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
write_export "${ENV_FILE}" "GROUND_TRUTH_TIMES" "${GROUND_TRUTH_TIMES}"
write_export "${ENV_FILE}" "PASSAGE_TIMES" "${PASSAGE_TIMES}"
write_export "${ENV_FILE}" "GROUND_TRUTH_REPS" "${GROUND_TRUTH_REPS}"
write_export "${ENV_FILE}" "FIT_REPEATS" "${FIT_REPEATS}"
write_export "${ENV_FILE}" "SOFT_MINOBS" "${SOFT_MINOBS}"
write_export "${ENV_FILE}" "NTP" "${NTP}"
write_export "${ENV_FILE}" "ALFAKR_DT" "${ALFAKR_DT}"
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

write_task_map() {
  local run_index="$1"
  local sample_depth="$2"
  local class="$3"
  local task_map="$4"
  awk -v target_depth="${sample_depth}" -v class="${class}" '
    BEGIN { FS = "\t" }
    NR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      required = "task_id sample_depth package input_mode extrapolation_method"
      nreq = split(required, req, " ")
      for (i = 1; i <= nreq; i++) {
        if (!col[req[i]]) {
          printf("Missing required run_index.tsv column: %s\n", req[i]) > "/dev/stderr"
          exit 2
        }
      }
      next
    }
    {
      sample_depth = $(col["sample_depth"]) + 0
      package_name = $(col["package"])
      input_mode = $(col["input_mode"])
      extrapolation_method = $(col["extrapolation_method"])
      is_graph_gaussian_full = (package_name == "alfak2" && input_mode == "full" && extrapolation_method == "graph_gaussian_baseline")
      keep = 0
      if (sample_depth == target_depth && class == "alfakR" && package_name == "alfakR") {
        keep = 1
      } else if (sample_depth == target_depth && class == "alfak2_ggfull" && is_graph_gaussian_full) {
        keep = 1
      } else if (sample_depth == target_depth && class == "alfak2_regular" && package_name == "alfak2" && !is_graph_gaussian_full) {
        keep = 1
      }
      if (keep) {
        print $(col["task_id"])
      }
    }
  ' "${run_index}" > "${task_map}"
}

submit_group() {
  local group="$1"
  local sample_depth="$2"
  local class="$3"
  local mem="$4"
  local task_map task_count array_spec group_job_name job_id

  task_map="${OUTPUT_DIR}/slurm/${group}.task_ids.tsv"
  write_task_map "${RUN_INDEX}" "${sample_depth}" "${class}" "${task_map}"
  task_count="$(awk 'END { print NR }' "${task_map}")"
  if [[ -z "${task_count}" || "${task_count}" == "0" ]]; then
    echo "Skipping ${group}: no matching tasks"
    return 0
  fi

  array_spec="1-${task_count}"
  if [[ -n "${ARRAY_LIMIT}" ]]; then
    array_spec="${array_spec}%${ARRAY_LIMIT}"
  fi

  group_job_name="${JOB_NAME}_${group}"
  echo "[$(date)] submitting group ${group}"
  echo "  sample_depth: ${sample_depth}"
  echo "  class:        ${class}"
  echo "  tasks:        ${task_count}"
  echo "  array:        ${array_spec}"
  echo "  cpu/task:     ${CPUS_PER_TASK}"
  echo "  mem/task:     ${mem}"
  echo "  time/task:    ${TIME_LIMIT}"

  SBATCH_ARGS=(
    "--parsable"
    "--array=${array_spec}"
    "--cpus-per-task=${CPUS_PER_TASK}"
    "--mem=${mem}"
    "--time=${TIME_LIMIT}"
    "--job-name=${group_job_name}"
    "--output=${OUTPUT_DIR}/slurm_logs/${group_job_name}_%A_%a.out"
    "--error=${OUTPUT_DIR}/slurm_logs/${group_job_name}_%A_%a.err"
    "--chdir=${ALFAK2_REPO}"
    "${SCRIPT_SELF}"
    "__run_task"
    "${ENV_FILE}"
    "${task_map}"
  )

  if [[ "${DRY_RUN}" == "true" ]]; then
    printf 'DRY-RUN: sbatch'
    printf ' %q' "${SBATCH_ARGS[@]}"
    printf '\n'
    job_id="DRY_RUN"
  else
    job_id=$(sbatch "${SBATCH_ARGS[@]}" | awk -F';' '{print $1}')
  fi

  printf '%s\n' "${job_id}" > "${OUTPUT_DIR}/slurm/${group}.jobid"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${group}" "${sample_depth}" "${class}" "${mem}" "${task_count}" "${array_spec}" "${task_map}" "${job_id}" \
    >> "${SUBMITTED_JOBS_TSV}"
  echo "Submitted ${group}: ${job_id}"
}

SUBMITTED_JOBS_TSV="${OUTPUT_DIR}/slurm/submitted_job_arrays.tsv"
printf 'group\tsample_depth\tclass\tmem\ttask_count\tarray_spec\ttask_map\tjob_id\n' > "${SUBMITTED_JOBS_TSV}"

echo "[$(date)] submitting alfa2_benchmark_ground_true in resource groups"
echo "  total tasks: ${N_TASKS}"
echo "  cpu/task:    ${CPUS_PER_TASK}"
echo "  time/task:   ${TIME_LIMIT}"
echo "  output dir:  ${OUTPUT_DIR}"
echo "  env file:    ${ENV_FILE}"
if [[ -n "${ARRAY_LIMIT}" ]]; then
  echo "  array limit: ${ARRAY_LIMIT}"
fi

submit_group "d200_k2_regular" 200 "alfak2_regular" "${MEM_DEPTH200_ALFAK2_REGULAR}"
submit_group "d200_kR" 200 "alfakR" "${MEM_DEPTH200_ALFAKR}"
submit_group "d200_k2_ggfull" 200 "alfak2_ggfull" "${MEM_DEPTH200_ALFAK2_GGFULL}"
submit_group "d1000_k2_regular" 1000 "alfak2_regular" "${MEM_DEPTH1000_ALFAK2_REGULAR}"
submit_group "d1000_kR" 1000 "alfakR" "${MEM_DEPTH1000_ALFAKR}"
submit_group "d1000_k2_ggfull" 1000 "alfak2_ggfull" "${MEM_DEPTH1000_ALFAK2_GGFULL}"

echo "Submission summary: ${SUBMITTED_JOBS_TSV}"
