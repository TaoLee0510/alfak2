#!/usr/bin/env bash

set -euo pipefail

STAGE="${1:-}"
ENV_FILE="${2:-}"

if [[ -z "${STAGE}" || -z "${ENV_FILE}" ]]; then
  echo "Usage: $0 STAGE ENV_FILE" >&2
  exit 2
fi
if [[ ! -f "${ENV_FILE}" ]]; then
  echo "Missing env file: ${ENV_FILE}" >&2
  exit 2
fi

source "${ENV_FILE}"

load_modules() {
  if [[ -n "${MODULES:-}" ]] && type module >/dev/null 2>&1; then
    module purge
    for module_name in ${MODULES}; do
      module load "${module_name}"
    done
  fi
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

extract_best_params_env() {
  local best_tsv="$1"
  if [[ ! -f "${best_tsv}" ]]; then
    echo "Missing best calibration params: ${best_tsv}" >&2
    exit 2
  fi
  Rscript -e '
args <- commandArgs(TRUE)
x <- read.table(args[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                check.names = FALSE, quote = "", comment.char = "", fill = TRUE)
if (!nrow(x)) stop("best_params.tsv is empty")
b <- x[1, , drop = FALSE]
field <- function(nm, default = "") {
  if (!nm %in% names(b)) return(default)
  value <- as.character(b[[nm]][[1]])
  if (!length(value) || is.na(value) || !nzchar(value)) default else value
}
emit <- function(name, value) {
  cat("export ", name, "=", shQuote(as.character(value), type = "sh"), "\n", sep = "")
}
emit("ALFAK2_INPUT_POLICIES", field("input_policy", "minobs_matched"))
emit("ALFAK2_INPUT_DEPTH", field("alfak2_input_depth", "effective"))
emit("ALFAK2_OBSERVATION_MODEL", field("alfak2_observation_model", "dirichlet_multinomial"))
emit("ALFAK2_LEGACY_WEIGHT", field("legacy_weight", "directly_informed"))
emit("CORRECT_EFFLUX", tolower(field("correct_efflux", "true")))
emit("ALFAK2_LAMBDA_L_GRID", field("lambda_l", "0.2"))
emit("ALFAK2_LAMBDA_E_GRID", field("lambda_e", "1"))
emit("ALFAK2_SIGMA_OBS_GRID", field("sigma_obs", "0.02"))
emit("ALFAK2_DM_CONCENTRATION", field("dm_concentration", "50"))
emit("ALFAK2_EFFECTIVE_DEPTH_MODE", field("effective_depth_mode", "min"))
emit("ALFAK2_LOCAL_SHELL_DEPTH", field("local_shell_depth", "0"))
emit("ALFAK2_GLOBAL_EXTRA_SHELL", field("global_extra_shell", "1"))
eff <- field("alfak2_effective_depth", "")
if (nzchar(eff) && is.finite(suppressWarnings(as.numeric(eff)))) {
  emit("ALFAK2_EFFECTIVE_DEPTH", eff)
}
' "${best_tsv}"
}

stage_shared_prepare() {
  export MODE=prepare
  export BASE_OUTPUT_DIR="${SHARED_INPUT_DIR}"
  export SOURCE_INPUT_DIR=""
  export PM="${TRUE_PM}"
  export BETA_LEVELS="${FIT_PM_LEVELS}"
  export N_SIM="${SHARED_N_SIM}"
  export FORCE_REFIT="${FORCE_REFIT}"
  export FORCE_SIM="${FORCE_SIM}"
  export REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE}"
  export RECOMPILE_DLL="${SHARED_RECOMPILE_DLL}"
  "${BENCHMARK_WRAPPER}"
}

stage_calibration_prepare() {
  export MODE=prepare
  export BASE_OUTPUT_DIR="${CAL_OUTPUT_DIR}"
  export SOURCE_INPUT_DIR="${SHARED_INPUT_DIR}"
  export PM="${FIT_PM}"
  export N_SIM="${CAL_N_SIM}"
  export FORCE_REFIT="${FORCE_REFIT}"
  export FORCE_SIM="${FORCE_SIM}"
  export REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE}"
  export RECOMPILE_DLL="${RECOMPILE_DLL}"
  "${CALIBRATION_WRAPPER}"
}

stage_launch_calibration_fit() {
  local n_tasks
  n_tasks="$(task_count_from_table "${CAL_OUTPUT_DIR}/tables/task_table.tsv")"

  export MODE=fit-task
  export BASE_OUTPUT_DIR="${CAL_OUTPUT_DIR}"
  export SOURCE_INPUT_DIR="${SHARED_INPUT_DIR}"
  export PM="${FIT_PM}"
  export N_SIM="${CAL_N_SIM}"
  export FORCE_REFIT="${FORCE_REFIT}"
  export FORCE_SIM="${FORCE_SIM}"
  export REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE}"
  export RECOMPILE_DLL="${RECOMPILE_DLL}"
  local fit_id
  fit_id=$(sbatch --parsable --export=ALL \
    --array="1-${n_tasks}%${CAL_ARRAY_LIMIT}" \
    --job-name="alfak2_cal_${FIT_PM_LABEL}" \
    --output="${LOG_DIR}/cal_${FIT_PM_LABEL}_%A_%a.out" \
    --error="${LOG_DIR}/cal_${FIT_PM_LABEL}_%A_%a.err" \
    "${CALIBRATION_WRAPPER}" | awk -F';' '{print $1}')
  printf "%s\n" "${fit_id}" > "${CONTROL_DIR}/cal_fit_${FIT_PM_LABEL}.jobid"

  export MODE=summarize
  local summary_id
  summary_id=$(sbatch --parsable --export=ALL \
    --dependency="afterok:${fit_id}" \
    --job-name="alfak2_cal_sum_${FIT_PM_LABEL}" \
    --output="${LOG_DIR}/cal_sum_${FIT_PM_LABEL}_%j.out" \
    --error="${LOG_DIR}/cal_sum_${FIT_PM_LABEL}_%j.err" \
    "${CALIBRATION_WRAPPER}" | awk -F';' '{print $1}')
  printf "%s\n" "${summary_id}" > "${CONTROL_DIR}/cal_summary_${FIT_PM_LABEL}.jobid"
  echo "Submitted calibration beta=${FIT_PM}: fit=${fit_id}, summarize=${summary_id}"
}

stage_benchmark_beta_driver() {
  eval "$(extract_best_params_env "${CAL_OUTPUT_DIR}/tables/best_params.tsv")"
  if [[ -f "${CAL_OUTPUT_DIR}/tables/best_params_cli_args.txt" ]]; then
    cp "${CAL_OUTPUT_DIR}/tables/best_params_cli_args.txt" \
      "${PIPELINE_TABLE_DIR}/best_params_cli_args_${FIT_PM_LABEL}.txt"
  fi

  export MODE=prepare
  export BASE_OUTPUT_DIR="${BENCH_OUTPUT_DIR}"
  export SOURCE_INPUT_DIR="${SHARED_INPUT_DIR}"
  export PM="${TRUE_PM}"
  export BETA_LEVELS="${FIT_PM}"
  export N_SIM="${BENCH_N_SIM}"
  export FORCE_REFIT="${FORCE_REFIT}"
  export FORCE_SIM="${FORCE_SIM}"
  export REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE}"
  export RECOMPILE_DLL="${RECOMPILE_DLL}"
  "${BENCHMARK_WRAPPER}"

  local n_tasks
  n_tasks="$(task_count_from_table "${BENCH_OUTPUT_DIR}/tables/task_table.tsv")"
  export MODE=fit-task
  local fit_id
  fit_id=$(sbatch --parsable --export=ALL \
    --array="1-${n_tasks}%${BENCH_ARRAY_LIMIT}" \
    --job-name="alfak2_bench_${FIT_PM_LABEL}" \
    --output="${LOG_DIR}/bench_${FIT_PM_LABEL}_%A_%a.out" \
    --error="${LOG_DIR}/bench_${FIT_PM_LABEL}_%A_%a.err" \
    "${BENCHMARK_WRAPPER}" | awk -F';' '{print $1}')
  printf "%s\n" "${fit_id}" > "${CONTROL_DIR}/bench_fit_${FIT_PM_LABEL}.jobid"

  export MODE=summarize
  local summary_id
  summary_id=$(sbatch --parsable --export=ALL \
    --dependency="afterok:${fit_id}" \
    --job-name="alfak2_bench_sum_${FIT_PM_LABEL}" \
    --output="${LOG_DIR}/bench_sum_${FIT_PM_LABEL}_%j.out" \
    --error="${LOG_DIR}/bench_sum_${FIT_PM_LABEL}_%j.err" \
    "${BENCHMARK_WRAPPER}" | awk -F';' '{print $1}')
  printf "%s\n" "${summary_id}" > "${CONTROL_DIR}/bench_summary_${FIT_PM_LABEL}.jobid"
  echo "Submitted benchmark beta=${FIT_PM}: fit=${fit_id}, summarize=${summary_id}"
}

stage_launch_benchmark_after_calibration() {
  local driver_ids=()
  local beta_env cal_summary_id driver_id summary_file
  while IFS= read -r beta_env; do
    [[ -z "${beta_env}" ]] && continue
    source "${beta_env}"
    summary_file="${CONTROL_DIR}/cal_summary_${FIT_PM_LABEL}.jobid"
    if [[ ! -f "${summary_file}" ]]; then
      echo "Missing calibration summary job id: ${summary_file}" >&2
      exit 2
    fi
    cal_summary_id="$(cat "${summary_file}")"
    driver_id=$(sbatch --parsable --export=ALL \
      --dependency="afterok:${cal_summary_id}" \
      --cpus-per-task=1 \
      --mem="${PREPARE_MEM}" \
      --time="${PREPARE_TIME}" \
      --job-name="alfak2_bench_launch_${FIT_PM_LABEL}" \
      --output="${LOG_DIR}/bench_launch_${FIT_PM_LABEL}_%j.out" \
      --error="${LOG_DIR}/bench_launch_${FIT_PM_LABEL}_%j.err" \
      "${STAGE_RUNNER}" "benchmark-beta-driver" "${beta_env}" | awk -F';' '{print $1}')
    printf "%s\n" "${driver_id}" > "${CONTROL_DIR}/bench_driver_${FIT_PM_LABEL}.jobid"
    driver_ids+=("${driver_id}")
  done < "${CONTROL_DIR}/beta_env_files.txt"

  local dep combine_launcher_id
  dep=$(IFS=:; echo "${driver_ids[*]}")
  combine_launcher_id=$(sbatch --parsable --export=ALL \
    --dependency="afterok:${dep}" \
    --cpus-per-task=1 \
    --mem="${LAUNCH_MEM}" \
    --time="${LAUNCH_TIME}" \
    --job-name="alfak2_beta_combine_launch" \
    --output="${LOG_DIR}/combine_launch_%j.out" \
    --error="${LOG_DIR}/combine_launch_%j.err" \
    "${STAGE_RUNNER}" "launch-combine-after-benchmark" "${CONTROL_DIR}/common_env.sh" | awk -F';' '{print $1}')
  printf "%s\n" "${combine_launcher_id}" > "${CONTROL_DIR}/combine_launcher.jobid"
  echo "Submitted benchmark drivers: ${driver_ids[*]}"
  echo "Submitted combine launcher: ${combine_launcher_id}"
}

stage_launch_combine_after_benchmark() {
  local summary_ids=()
  local beta_env summary_file
  while IFS= read -r beta_env; do
    [[ -z "${beta_env}" ]] && continue
    source "${beta_env}"
    summary_file="${CONTROL_DIR}/bench_summary_${FIT_PM_LABEL}.jobid"
    if [[ ! -f "${summary_file}" ]]; then
      echo "Missing benchmark summary job id: ${summary_file}" >&2
      exit 2
    fi
    summary_ids+=("$(cat "${summary_file}")")
  done < "${CONTROL_DIR}/beta_env_files.txt"
  local dep combine_id
  dep=$(IFS=:; echo "${summary_ids[*]}")
  combine_id=$(sbatch --parsable --export=ALL \
    --dependency="afterok:${dep}" \
    --cpus-per-task=1 \
    --mem="${COMBINE_MEM}" \
    --time="${COMBINE_TIME}" \
    --job-name="alfak2_beta_combine" \
    --output="${LOG_DIR}/combine_%j.out" \
    --error="${LOG_DIR}/combine_%j.err" \
    "${STAGE_RUNNER}" "combine" "${CONTROL_DIR}/common_env.sh" | awk -F';' '{print $1}')
  printf "%s\n" "${combine_id}" > "${CONTROL_DIR}/combine.jobid"
  echo "Submitted final combine job: ${combine_id}"
}

stage_combine() {
  load_modules
  cd "${ALFAK2_REPO}"
  Rscript "${COMBINE_SCRIPT}" \
    "--pipeline-dir=${PIPELINE_OUTPUT_DIR}" \
    "--output-dir=${COMBINED_OUTPUT_DIR}"
}

case "${STAGE}" in
  shared-prepare) stage_shared_prepare ;;
  calibration-prepare) stage_calibration_prepare ;;
  launch-calibration-fit) stage_launch_calibration_fit ;;
  benchmark-beta-driver) stage_benchmark_beta_driver ;;
  launch-benchmark-after-calibration) stage_launch_benchmark_after_calibration ;;
  launch-combine-after-benchmark) stage_launch_combine_after_benchmark ;;
  combine) stage_combine ;;
  *)
    echo "Unsupported stage: ${STAGE}" >&2
    exit 2
    ;;
esac
