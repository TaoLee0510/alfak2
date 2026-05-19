#!/usr/bin/env bash
#
# Submit the GRF-reference downsampled-input benchmark to SLURM.
#
# This is a compact benchmark pipeline for the case where alfak2 does not need a
# fitted-beta grid. If SHARED_INPUT_DIR already contains tables/input_table.tsv,
# the script reuses those GRF/ABM-derived downsampled inputs directly. Otherwise,
# it prepares those inputs once, then builds the alfak2-vs-alfakR benchmark task
# table from the same input files. alfakR is run as the comparison engine, but
# ABM/input generation is not repeated for alfakR or for separate beta levels.
#
# Usage from a login node:
#
#   bash benchmark/scr/submit_grf_downsampled_input_benchmark_slurm.sh
#
# Useful overrides:
#
#   OUTPUT_DIR=/path/to/results \
#   SHARED_INPUT_DIR=/path/to/existing/shared_inputs \
#   FIT_PM=5e-05 \
#   METHODS=empirical \
#   bash benchmark/scr/submit_grf_downsampled_input_benchmark_slurm.sh

set -euo pipefail

SCRIPT_SELF="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"

write_export() {
  local file="$1"
  local name="$2"
  local value="${3-}"
  printf 'export %s=%q\n' "${name}" "${value}" >> "${file}"
}

task_count_from_table() {
  local task_table="$1"
  if [[ ! -f "${task_table}" ]]; then
    echo "Missing task table: ${task_table}" >&2
    exit 2
  fi
  local n_tasks
  n_tasks=$(($(wc -l < "${task_table}") - 1))
  if (( n_tasks < 1 )); then
    echo "Task table has no fit tasks: ${task_table}" >&2
    exit 2
  fi
  printf "%s" "${n_tasks}"
}

task_array_spec_for_engine() {
  local task_table="$1"
  local engine="$2"
  awk -v engine="${engine}" '
    BEGIN { FS = "\t" }
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "engine") engine_col = i
      }
      if (!engine_col) exit 2
      next
    }
    $engine_col == engine {
      id = NR - 1
      n++
      if (start == "") {
        start = id
        prev = id
        next
      }
      if (id == prev + 1) {
        prev = id
        next
      }
      if (spec != "") spec = spec ","
      spec = spec (start == prev ? start : start "-" prev)
      start = id
      prev = id
    }
    END {
      if (start != "") {
        if (spec != "") spec = spec ","
        spec = spec (start == prev ? start : start "-" prev)
      }
      print spec
    }
  ' "${task_table}"
}

submit_sbatch() {
  if [[ "${DRY_RUN:-false}" == "true" ]]; then
    echo "DRY-RUN: sbatch $*" >&2
    echo "999999"
  else
    sbatch --parsable "$@" | awk -F';' '{print $1}'
  fi
}

run_wrapper_stage() {
  local env_file="${1:-}"
  if [[ -z "${env_file}" || ! -f "${env_file}" ]]; then
    echo "Usage: $0 __run_wrapper ENV_FILE" >&2
    exit 2
  fi
  source "${env_file}"
  "${BENCHMARK_WRAPPER}"
}

launch_fit_array() {
  local fit_env="${1:-}"
  local summary_env="${2:-}"
  if [[ -z "${fit_env}" || ! -f "${fit_env}" || -z "${summary_env}" || ! -f "${summary_env}" ]]; then
    echo "Usage: $0 __launch_fit_array FIT_ENV SUMMARY_ENV" >&2
    exit 2
  fi
  source "${fit_env}"

  local task_table n_tasks alfak2_array_spec alfakR_array_spec alfak2_fit_id alfakR_fit_id summary_id
  local fit_dependencies=()
  task_table="${BASE_OUTPUT_DIR}/tables/task_table.tsv"
  n_tasks="$(task_count_from_table "${task_table}")"

  alfak2_array_spec="$(task_array_spec_for_engine "${task_table}" "alfak2")"
  alfakR_array_spec="$(task_array_spec_for_engine "${task_table}" "alfakR")"

  if [[ -n "${alfak2_array_spec}" ]]; then
    alfak2_fit_id=$(submit_sbatch \
      --export=ALL \
      --array="${alfak2_array_spec}%${FIT_ARRAY_LIMIT}" \
      --cpus-per-task=1 \
      --mem="${ALFAK2_FIT_MEM}" \
      --time="${FIT_TIME}" \
      --job-name="alfak2_grf_down_fit_k2" \
      --output="${LOG_DIR}/fit_alfak2_%A_%a.out" \
      --error="${LOG_DIR}/fit_alfak2_%A_%a.err" \
      "${SCRIPT_SELF}" "__run_wrapper" "${fit_env}")
    printf "%s\n" "${alfak2_fit_id}" > "${CONTROL_DIR}/fit_alfak2_array.jobid"
    fit_dependencies+=("${alfak2_fit_id}")
  else
    alfak2_fit_id="SKIPPED"
    printf "%s\n" "${alfak2_fit_id}" > "${CONTROL_DIR}/fit_alfak2_array.jobid"
  fi

  if [[ -n "${alfakR_array_spec}" ]]; then
    alfakR_fit_id=$(submit_sbatch \
      --export=ALL \
      --array="${alfakR_array_spec}%${FIT_ARRAY_LIMIT}" \
      --cpus-per-task=1 \
      --mem="${ALFAKR_FIT_MEM}" \
      --time="${FIT_TIME}" \
      --job-name="alfak2_grf_down_fit_ar" \
      --output="${LOG_DIR}/fit_alfakR_%A_%a.out" \
      --error="${LOG_DIR}/fit_alfakR_%A_%a.err" \
      "${SCRIPT_SELF}" "__run_wrapper" "${fit_env}")
    printf "%s\n" "${alfakR_fit_id}" > "${CONTROL_DIR}/fit_alfakR_array.jobid"
    fit_dependencies+=("${alfakR_fit_id}")
  else
    alfakR_fit_id="SKIPPED"
    printf "%s\n" "${alfakR_fit_id}" > "${CONTROL_DIR}/fit_alfakR_array.jobid"
  fi

  if (( ${#fit_dependencies[@]} < 1 )); then
    echo "No alfak2 or alfakR fit tasks found in task table: ${task_table}" >&2
    exit 2
  fi
  printf "%s\n" "${fit_dependencies[@]}" > "${CONTROL_DIR}/fit_array.jobid"

  local dependency_arg
  dependency_arg="$(IFS=:; printf "%s" "${fit_dependencies[*]}")"

  summary_id=$(submit_sbatch \
    --export=ALL \
    --dependency="afterok:${dependency_arg}" \
    --cpus-per-task=1 \
    --mem="${SUMMARY_MEM}" \
    --time="${SUMMARY_TIME}" \
    --job-name="alfak2_grf_downsample_sum" \
    --output="${LOG_DIR}/summary_%j.out" \
    --error="${LOG_DIR}/summary_%j.err" \
    "${SCRIPT_SELF}" "__run_wrapper" "${summary_env}")
  printf "%s\n" "${summary_id}" > "${CONTROL_DIR}/summary.jobid"

  printf "\nSubmitted GRF downsampled benchmark fits.\n\n"
  printf "Task table:       %s\n" "${task_table}"
  printf "Fit tasks:        %s\n" "${n_tasks}"
  printf "alfak2 task spec: %s\n" "${alfak2_array_spec:-NONE}"
  printf "alfak2 fit mem:   %s\n" "${ALFAK2_FIT_MEM}"
  printf "alfak2 fit job:   %s\n" "${alfak2_fit_id}"
  printf "alfakR task spec: %s\n" "${alfakR_array_spec:-NONE}"
  printf "alfakR fit mem:   %s\n" "${ALFAKR_FIT_MEM}"
  printf "alfakR fit job:   %s\n" "${alfakR_fit_id}"
  printf "Summary job:      %s\n" "${summary_id}"
}

case "${1:-}" in
  __run_wrapper)
    run_wrapper_stage "${2:-}"
    exit
    ;;
  __launch_fit_array)
    launch_fit_array "${2:-}" "${3:-}"
    exit
    ;;
esac

ALFAK2_REPO="${ALFAK2_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfak2}"
ALFAKR_REPO="${ALFAKR_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfakR}"
OUTPUT_DIR="${OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/grf_downsampled_input_benchmark}"
SHARED_INPUT_DIR="${SHARED_INPUT_DIR:-${OUTPUT_DIR}/shared_inputs}"
BENCHMARK_OUTPUT_DIR="${BENCHMARK_OUTPUT_DIR:-${OUTPUT_DIR}/benchmark}"
BENCHMARK_WRAPPER="${BENCHMARK_WRAPPER:-${ALFAK2_REPO}/benchmark/scr/submit_grf_alfak2_vs_alfakR_benchmark_slurm.sh}"

MODULES="${MODULES:-R/4.4}"
TRUE_PM="${TRUE_PM:-${PM:-5e-05}}"
FIT_PM="${FIT_PM:-${TRUE_PM}}"

METHODS="${METHODS:-none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step}"
MINOBS="${MINOBS:-5,10,20}"
N_SIM="${N_SIM:-12}"
SHARED_N_SIM="${SHARED_N_SIM:-${N_SIM}}"
BENCH_N_SIM="${BENCH_N_SIM:-${N_SIM}}"
LAMBDAS="${LAMBDAS:-0.2,0.4,0.8,1.6}"
TIME_GAPS="${TIME_GAPS:-2,4,8}"
TIME_STARTS="${TIME_STARTS:-0}"
NBOOT="${NBOOT:-45}"
GRID_N="${GRID_N:-81}"
SAMPLE_DEPTH="${SAMPLE_DEPTH:-2000}"
TIME_MAX="${TIME_MAX:-360}"
PASSAGE_INTERVAL="${PASSAGE_INTERVAL:-45}"

ALFAK2_INPUT_POLICIES="${ALFAK2_INPUT_POLICIES:-minobs_matched}"
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
REUSE_SHARED_INPUTS="${REUSE_SHARED_INPUTS:-auto}"
SHARED_RECOMPILE_DLL="${SHARED_RECOMPILE_DLL:-true}"
RECOMPILE_DLL="${RECOMPILE_DLL:-false}"
EXTRA_ARGS="${EXTRA_ARGS:-}"

PREPARE_TIME="${PREPARE_TIME:-24:00:00}"
PREPARE_MEM="${PREPARE_MEM:-8G}"
LAUNCH_TIME="${LAUNCH_TIME:-02:00:00}"
LAUNCH_MEM="${LAUNCH_MEM:-1G}"
FIT_TIME="${FIT_TIME:-48:00:00}"
ALFAK2_FIT_MEM="${ALFAK2_FIT_MEM:-16G}"
ALFAKR_FIT_MEM="${ALFAKR_FIT_MEM:-32G}"
SUMMARY_TIME="${SUMMARY_TIME:-06:00:00}"
SUMMARY_MEM="${SUMMARY_MEM:-8G}"
FIT_ARRAY_LIMIT="${FIT_ARRAY_LIMIT:-64}"
DRY_RUN="${DRY_RUN:-false}"

CONTROL_DIR="${OUTPUT_DIR}/control"
LOG_DIR="${OUTPUT_DIR}/slurm_logs"
SHARED_INPUT_TABLE="${SHARED_INPUT_DIR}/tables/input_table.tsv"

mkdir -p "${CONTROL_DIR}" "${LOG_DIR}" "${BENCHMARK_OUTPUT_DIR}"

case "${REUSE_SHARED_INPUTS}" in
  auto|true|false) ;;
  *)
    echo "REUSE_SHARED_INPUTS must be auto, true, or false; got ${REUSE_SHARED_INPUTS}" >&2
    exit 2
    ;;
esac

skip_shared_prepare=false
if [[ "${REUSE_SHARED_INPUTS}" == "true" ]]; then
  if [[ ! -f "${SHARED_INPUT_TABLE}" ]]; then
    echo "REUSE_SHARED_INPUTS=true but missing shared input table: ${SHARED_INPUT_TABLE}" >&2
    exit 2
  fi
  skip_shared_prepare=true
elif [[ "${REUSE_SHARED_INPUTS}" == "auto" && "${FORCE_SIM}" != "true" && -f "${SHARED_INPUT_TABLE}" ]]; then
  skip_shared_prepare=true
fi

if [[ "${skip_shared_prepare}" != "true" ]]; then
  mkdir -p "${SHARED_INPUT_DIR}"
fi

COMMON_ENV="${CONTROL_DIR}/common_env.sh"
SHARED_PREPARE_ENV="${CONTROL_DIR}/shared_prepare.env"
BENCHMARK_PREPARE_ENV="${CONTROL_DIR}/benchmark_prepare.env"
BENCHMARK_FIT_ENV="${CONTROL_DIR}/benchmark_fit.env"
BENCHMARK_SUMMARY_ENV="${CONTROL_DIR}/benchmark_summary.env"

: > "${COMMON_ENV}"
for name in \
  SCRIPT_SELF ALFAK2_REPO ALFAKR_REPO OUTPUT_DIR SHARED_INPUT_DIR BENCHMARK_OUTPUT_DIR \
  BENCHMARK_WRAPPER MODULES TRUE_PM FIT_PM METHODS MINOBS LAMBDAS TIME_GAPS TIME_STARTS \
  NBOOT GRID_N SAMPLE_DEPTH TIME_MAX PASSAGE_INTERVAL ALFAK2_INPUT_POLICIES \
  ALFAK2_INPUT_DEPTH ALFAK2_OBSERVATION_MODEL ALFAK2_DM_CONCENTRATION \
  ALFAK2_EFFECTIVE_DEPTH_MODE ALFAK2_EFFECTIVE_DEPTH ALFAK2_LEGACY_WEIGHT CORRECT_EFFLUX \
  ALFAK2_LAMBDA_L_GRID ALFAK2_LAMBDA_E_GRID ALFAK2_SIGMA_OBS_GRID ALFAK2_GRAPH_EDGE_WEIGHT \
  ALFAK2_ANCHOR_COUNT_REFERENCE ALFAK2_ANCHOR_COUNT_POWER ALFAK2_LOCAL_SHELL_DEPTH \
  ALFAK2_GLOBAL_EXTRA_SHELL ALFAK2_MAX_NODES FORCE_REFIT REUSE_DIRTY_CACHE \
  REUSE_SHARED_INPUTS SHARED_INPUT_TABLE \
  EXTRA_ARGS PREPARE_TIME PREPARE_MEM LAUNCH_TIME LAUNCH_MEM FIT_TIME ALFAK2_FIT_MEM ALFAKR_FIT_MEM \
  SUMMARY_TIME SUMMARY_MEM FIT_ARRAY_LIMIT DRY_RUN CONTROL_DIR LOG_DIR
do
  write_export "${COMMON_ENV}" "${name}" "${!name}"
done

write_stage_env() {
  local file="$1"
  local mode="$2"
  local base_output_dir="$3"
  local source_input_dir="$4"
  local n_sim="$5"
  local force_sim="$6"
  local recompile_dll="$7"

  : > "${file}"
  printf 'source %q\n' "${COMMON_ENV}" >> "${file}"
  write_export "${file}" "MODE" "${mode}"
  write_export "${file}" "BASE_OUTPUT_DIR" "${base_output_dir}"
  write_export "${file}" "SOURCE_INPUT_DIR" "${source_input_dir}"
  write_export "${file}" "N_SIM" "${n_sim}"
  write_export "${file}" "PM" "${TRUE_PM}"
  write_export "${file}" "BETA_LEVELS" "${FIT_PM}"
  write_export "${file}" "FORCE_SIM" "${force_sim}"
  write_export "${file}" "RECOMPILE_DLL" "${recompile_dll}"
  write_export "${file}" "ALFAK2_ERROR_LOG_DIR" "${base_output_dir}/error_logs"
  write_export "${file}" "ALFAK2_WARNING_LOG_DIR" "${base_output_dir}/error_logs"
}

write_stage_env "${SHARED_PREPARE_ENV}" "prepare" "${SHARED_INPUT_DIR}" "" "${SHARED_N_SIM}" "${FORCE_SIM}" "${SHARED_RECOMPILE_DLL}"
write_stage_env "${BENCHMARK_PREPARE_ENV}" "prepare" "${BENCHMARK_OUTPUT_DIR}" "${SHARED_INPUT_DIR}" "${BENCH_N_SIM}" "false" "${RECOMPILE_DLL}"
write_stage_env "${BENCHMARK_FIT_ENV}" "fit-task" "${BENCHMARK_OUTPUT_DIR}" "${SHARED_INPUT_DIR}" "${BENCH_N_SIM}" "false" "${RECOMPILE_DLL}"
write_stage_env "${BENCHMARK_SUMMARY_ENV}" "summarize" "${BENCHMARK_OUTPUT_DIR}" "${SHARED_INPUT_DIR}" "${BENCH_N_SIM}" "false" "${RECOMPILE_DLL}"

if [[ "${skip_shared_prepare}" == "true" ]]; then
  shared_prepare_id="SKIPPED"
else
  shared_prepare_id=$(submit_sbatch \
    --export=ALL \
    --cpus-per-task=1 \
    --mem="${PREPARE_MEM}" \
    --time="${PREPARE_TIME}" \
    --job-name="alfak2_grf_downsample_input" \
    --output="${LOG_DIR}/shared_prepare_%j.out" \
    --error="${LOG_DIR}/shared_prepare_%j.err" \
    "${SCRIPT_SELF}" "__run_wrapper" "${SHARED_PREPARE_ENV}")
fi
printf "%s\n" "${shared_prepare_id}" > "${CONTROL_DIR}/shared_prepare.jobid"

benchmark_prepare_dependency=()
if [[ "${skip_shared_prepare}" != "true" ]]; then
  benchmark_prepare_dependency=("--dependency=afterok:${shared_prepare_id}")
fi
benchmark_prepare_id=$(submit_sbatch \
  --export=ALL \
  "${benchmark_prepare_dependency[@]}" \
  --cpus-per-task=1 \
  --mem="${PREPARE_MEM}" \
  --time="${PREPARE_TIME}" \
  --job-name="alfak2_grf_downsample_prepare" \
  --output="${LOG_DIR}/benchmark_prepare_%j.out" \
  --error="${LOG_DIR}/benchmark_prepare_%j.err" \
  "${SCRIPT_SELF}" "__run_wrapper" "${BENCHMARK_PREPARE_ENV}")
printf "%s\n" "${benchmark_prepare_id}" > "${CONTROL_DIR}/benchmark_prepare.jobid"

launch_id=$(submit_sbatch \
  --export=ALL \
  --dependency="afterok:${benchmark_prepare_id}" \
  --cpus-per-task=1 \
  --mem="${LAUNCH_MEM}" \
  --time="${LAUNCH_TIME}" \
  --job-name="alfak2_grf_downsample_launch" \
  --output="${LOG_DIR}/launch_fit_%j.out" \
  --error="${LOG_DIR}/launch_fit_%j.err" \
  "${SCRIPT_SELF}" "__launch_fit_array" "${BENCHMARK_FIT_ENV}" "${BENCHMARK_SUMMARY_ENV}")
printf "%s\n" "${launch_id}" > "${CONTROL_DIR}/launch_fit.jobid"

printf "\nSubmitted GRF-reference downsampled-input benchmark.\n\n"
printf "Output root:           %s\n" "${OUTPUT_DIR}"
printf "Shared input dir:      %s\n" "${SHARED_INPUT_DIR}"
printf "Benchmark output dir:  %s\n" "${BENCHMARK_OUTPUT_DIR}"
printf "True PM:               %s\n" "${TRUE_PM}"
printf "Fit PM:                %s\n" "${FIT_PM}"
printf "Methods:               %s\n" "${METHODS}"
printf "alfak2 fit memory:     %s\n" "${ALFAK2_FIT_MEM}"
printf "alfakR fit memory:     %s\n" "${ALFAKR_FIT_MEM}"
printf "alfak2 input policies: %s\n\n" "${ALFAK2_INPUT_POLICIES}"
if [[ "${skip_shared_prepare}" == "true" ]]; then
  printf "Shared prepare job:    SKIPPED, reusing %s\n" "${SHARED_INPUT_TABLE}"
else
  printf "Shared prepare job:    %s\n" "${shared_prepare_id}"
fi
printf "Benchmark prepare job: %s\n" "${benchmark_prepare_id}"
printf "Fit launcher job:      %s\n" "${launch_id}"
printf "Logs:                  %s\n\n" "${LOG_DIR}"
printf "The fit launcher will read:\n"
printf "  %s/tables/task_table.tsv\n" "${BENCHMARK_OUTPUT_DIR}"
printf "and submit the final fit array plus summarize job after benchmark prepare finishes.\n"
