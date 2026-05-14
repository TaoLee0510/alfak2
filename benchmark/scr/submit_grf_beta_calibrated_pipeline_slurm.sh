#!/usr/bin/env bash
#
# Submit the beta-calibrated GRF benchmark pipeline to SLURM.
#
# This script is meant to be run from a login node:
#
#   bash benchmark/scr/submit_grf_beta_calibrated_pipeline_slurm.sh
#
# It submits the following dependency chain:
#   shared input prepare ->
#   per-beta alfak2 calibration prepare/fit/summarize ->
#   per-beta benchmark prepare/fit/summarize with that beta's best alfak2 params ->
#   combined result table build.

set -euo pipefail

ALFAK2_REPO="${ALFAK2_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfak2}"
ALFAKR_REPO="${ALFAKR_REPO:-/share/lab_crd/lab_crd/taoli/Project/alfakR}"
PIPELINE_OUTPUT_DIR="${PIPELINE_OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/grf_alfak2_beta_calibrated_pipeline}"

TRUE_PM="${TRUE_PM:-5e-05}"
FIT_PM_LEVELS="${FIT_PM_LEVELS:-1e-05,5e-05,1e-04,1e-03,1e-02}"

SHARED_INPUT_DIR="${SHARED_INPUT_DIR:-${PIPELINE_OUTPUT_DIR}/shared_inputs}"
CALIBRATION_ROOT="${CALIBRATION_ROOT:-${PIPELINE_OUTPUT_DIR}/calibration}"
BENCHMARK_ROOT="${BENCHMARK_ROOT:-${PIPELINE_OUTPUT_DIR}/benchmark}"
COMBINED_OUTPUT_DIR="${COMBINED_OUTPUT_DIR:-${PIPELINE_OUTPUT_DIR}/combined}"

MODULES="${MODULES:-R/4.4}"
BENCHMARK_WRAPPER="${BENCHMARK_WRAPPER:-${ALFAK2_REPO}/benchmark/scr/submit_grf_alfak2_vs_alfakR_benchmark_slurm.sh}"
CALIBRATION_WRAPPER="${CALIBRATION_WRAPPER:-${ALFAK2_REPO}/benchmark/scr/submit_grf_alfak2_parameter_calibration_slurm.sh}"
COMBINE_SCRIPT="${COMBINE_SCRIPT:-${ALFAK2_REPO}/benchmark/scr/combine_grf_beta_pipeline_results.R}"
STAGE_RUNNER="${STAGE_RUNNER:-${ALFAK2_REPO}/benchmark/scr/run_grf_beta_pipeline_stage.sh}"

METHODS="${METHODS:-none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step}"
MINOBS="${MINOBS:-5,10,20}"
BENCH_N_SIM="${BENCH_N_SIM:-12}"
CAL_N_SIM="${CAL_N_SIM:-3}"
SHARED_N_SIM="${SHARED_N_SIM:-${BENCH_N_SIM}}"
LAMBDAS="${LAMBDAS:-0.2,0.4,0.8,1.6}"
TIME_GAPS="${TIME_GAPS:-2,4,8}"
TIME_STARTS="${TIME_STARTS:-0}"
NBOOT="${NBOOT:-45}"
GRID_N="${GRID_N:-81}"
SAMPLE_DEPTH="${SAMPLE_DEPTH:-2000}"
TIME_MAX="${TIME_MAX:-360}"
PASSAGE_INTERVAL="${PASSAGE_INTERVAL:-45}"

INPUT_POLICY="${INPUT_POLICY:-minobs_matched}"
ALFAK2_INPUT_DEPTH="${ALFAK2_INPUT_DEPTH:-effective}"
ALFAK2_OBSERVATION_MODEL="${ALFAK2_OBSERVATION_MODEL:-dirichlet_multinomial}"
ALFAK2_EFFECTIVE_DEPTH="${ALFAK2_EFFECTIVE_DEPTH:-}"
ALFAK2_MAX_NODES="${ALFAK2_MAX_NODES:-150000}"

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
OBJECTIVE_METRIC="${OBJECTIVE_METRIC:-mae}"
DIRECT_WEIGHT="${DIRECT_WEIGHT:-0.25}"
BIAS_WEIGHT="${BIAS_WEIGHT:-0.10}"

CAL_ARRAY_LIMIT="${CAL_ARRAY_LIMIT:-15}"
BENCH_ARRAY_LIMIT="${BENCH_ARRAY_LIMIT:-64}"
PREPARE_TIME="${PREPARE_TIME:-24:00:00}"
PREPARE_MEM="${PREPARE_MEM:-8G}"
LAUNCH_TIME="${LAUNCH_TIME:-02:00:00}"
LAUNCH_MEM="${LAUNCH_MEM:-1G}"
COMBINE_TIME="${COMBINE_TIME:-06:00:00}"
COMBINE_MEM="${COMBINE_MEM:-8G}"

FORCE_REFIT="${FORCE_REFIT:-false}"
FORCE_SIM="${FORCE_SIM:-false}"
REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE:-false}"
SHARED_RECOMPILE_DLL="${SHARED_RECOMPILE_DLL:-true}"
RECOMPILE_DLL="${RECOMPILE_DLL:-false}"
DRY_RUN="${DRY_RUN:-false}"

CONTROL_DIR="${PIPELINE_OUTPUT_DIR}/control"
PIPELINE_TABLE_DIR="${PIPELINE_OUTPUT_DIR}/tables"
LOG_DIR="${PIPELINE_OUTPUT_DIR}/slurm_logs"

mkdir -p "${CONTROL_DIR}" "${PIPELINE_TABLE_DIR}" "${LOG_DIR}" "${CALIBRATION_ROOT}" "${BENCHMARK_ROOT}" "${COMBINED_OUTPUT_DIR}"

trim_spaces() {
  local x="$1"
  x="${x#"${x%%[![:space:]]*}"}"
  x="${x%"${x##*[![:space:]]}"}"
  printf "%s" "${x}"
}

beta_label() {
  printf "%s" "$1" | sed 's/[^A-Za-z0-9]/_/g'
}

write_export() {
  local file="$1"
  local name="$2"
  local value="${3-}"
  printf 'export %s=%q\n' "${name}" "${value}" >> "${file}"
}

submit_sbatch() {
  if [[ "${DRY_RUN}" == "true" ]]; then
    echo "DRY-RUN: sbatch $*" >&2
    echo "999999"
  else
    sbatch --parsable "$@" | awk -F';' '{print $1}'
  fi
}

COMMON_ENV="${CONTROL_DIR}/common_env.sh"
: > "${COMMON_ENV}"
for name in \
  ALFAK2_REPO ALFAKR_REPO PIPELINE_OUTPUT_DIR TRUE_PM FIT_PM_LEVELS \
  SHARED_INPUT_DIR CALIBRATION_ROOT BENCHMARK_ROOT COMBINED_OUTPUT_DIR \
  MODULES BENCHMARK_WRAPPER CALIBRATION_WRAPPER COMBINE_SCRIPT STAGE_RUNNER \
  METHODS MINOBS BENCH_N_SIM CAL_N_SIM SHARED_N_SIM LAMBDAS TIME_GAPS \
  TIME_STARTS NBOOT GRID_N SAMPLE_DEPTH TIME_MAX PASSAGE_INTERVAL \
  INPUT_POLICY ALFAK2_INPUT_DEPTH ALFAK2_OBSERVATION_MODEL ALFAK2_EFFECTIVE_DEPTH \
  ALFAK2_MAX_NODES LEGACY_WEIGHTS CORRECT_EFFLUX_VALUES LAMBDA_L_VALUES \
  LAMBDA_E_VALUES SIGMA_OBS_VALUES DM_CONCENTRATIONS EFFECTIVE_DEPTH_MODES \
  LOCAL_SHELL_DEPTHS GLOBAL_EXTRA_SHELLS OBJECTIVE_SCOPE OBJECTIVE_METRIC \
  DIRECT_WEIGHT BIAS_WEIGHT CAL_ARRAY_LIMIT BENCH_ARRAY_LIMIT PREPARE_TIME \
  PREPARE_MEM LAUNCH_TIME LAUNCH_MEM COMBINE_TIME COMBINE_MEM FORCE_REFIT \
  FORCE_SIM REUSE_DIRTY_CACHE SHARED_RECOMPILE_DLL RECOMPILE_DLL CONTROL_DIR \
  PIPELINE_TABLE_DIR LOG_DIR
do
  write_export "${COMMON_ENV}" "${name}" "${!name}"
done

IFS=',' read -r -a RAW_BETAS <<< "${FIT_PM_LEVELS}"
BETA_ENV_LIST="${CONTROL_DIR}/beta_env_files.txt"
MANIFEST="${PIPELINE_TABLE_DIR}/beta_manifest.tsv"
: > "${BETA_ENV_LIST}"
printf "pm\tfit_beta_label\tcalibration_dir\tbenchmark_dir\n" > "${MANIFEST}"

for raw_beta in "${RAW_BETAS[@]}"; do
  fit_pm="$(trim_spaces "${raw_beta}")"
  if [[ -z "${fit_pm}" ]]; then
    continue
  fi
  fit_label="$(beta_label "${fit_pm}")"
  beta_env="${CONTROL_DIR}/beta_${fit_label}.env"
  cal_dir="${CALIBRATION_ROOT}/beta_${fit_label}"
  bench_dir="${BENCHMARK_ROOT}/beta_${fit_label}"
  mkdir -p "${cal_dir}" "${bench_dir}"
  : > "${beta_env}"
  printf 'source %q\n' "${COMMON_ENV}" >> "${beta_env}"
  write_export "${beta_env}" "FIT_PM" "${fit_pm}"
  write_export "${beta_env}" "FIT_PM_LABEL" "${fit_label}"
  write_export "${beta_env}" "CAL_OUTPUT_DIR" "${cal_dir}"
  write_export "${beta_env}" "BENCH_OUTPUT_DIR" "${bench_dir}"
  printf "%s\n" "${beta_env}" >> "${BETA_ENV_LIST}"
  printf "%s\t%s\t%s\t%s\n" "${fit_pm}" "${fit_label}" "${cal_dir}" "${bench_dir}" >> "${MANIFEST}"
done

if [[ ! -s "${BETA_ENV_LIST}" ]]; then
  echo "No beta levels parsed from FIT_PM_LEVELS=${FIT_PM_LEVELS}" >&2
  exit 2
fi

shared_prepare_id=$(submit_sbatch \
  --export=ALL \
  --cpus-per-task=1 \
  --mem="${PREPARE_MEM}" \
  --time="${PREPARE_TIME}" \
  --job-name="alfak2_beta_shared_prepare" \
  --output="${LOG_DIR}/shared_prepare_%j.out" \
  --error="${LOG_DIR}/shared_prepare_%j.err" \
  "${STAGE_RUNNER}" "shared-prepare" "${COMMON_ENV}")
printf "%s\n" "${shared_prepare_id}" > "${CONTROL_DIR}/shared_prepare.jobid"

cal_launch_ids=()
while IFS= read -r beta_env; do
  source "${beta_env}"
  cal_prepare_id=$(submit_sbatch \
    --export=ALL \
    --dependency="afterok:${shared_prepare_id}" \
    --cpus-per-task=1 \
    --mem="${PREPARE_MEM}" \
    --time="${PREPARE_TIME}" \
    --job-name="alfak2_cal_prepare_${FIT_PM_LABEL}" \
    --output="${LOG_DIR}/cal_prepare_${FIT_PM_LABEL}_%j.out" \
    --error="${LOG_DIR}/cal_prepare_${FIT_PM_LABEL}_%j.err" \
    "${STAGE_RUNNER}" "calibration-prepare" "${beta_env}")
  printf "%s\n" "${cal_prepare_id}" > "${CONTROL_DIR}/cal_prepare_${FIT_PM_LABEL}.jobid"

  cal_launch_id=$(submit_sbatch \
    --export=ALL \
    --dependency="afterok:${cal_prepare_id}" \
    --cpus-per-task=1 \
    --mem="${LAUNCH_MEM}" \
    --time="${LAUNCH_TIME}" \
    --job-name="alfak2_cal_launch_${FIT_PM_LABEL}" \
    --output="${LOG_DIR}/cal_launch_${FIT_PM_LABEL}_%j.out" \
    --error="${LOG_DIR}/cal_launch_${FIT_PM_LABEL}_%j.err" \
    "${STAGE_RUNNER}" "launch-calibration-fit" "${beta_env}")
  printf "%s\n" "${cal_launch_id}" > "${CONTROL_DIR}/cal_launch_${FIT_PM_LABEL}.jobid"
  cal_launch_ids+=("${cal_launch_id}")
done < "${BETA_ENV_LIST}"

cal_launch_dep=$(IFS=:; echo "${cal_launch_ids[*]}")
benchmark_launcher_id=$(submit_sbatch \
  --export=ALL \
  --dependency="afterok:${cal_launch_dep}" \
  --cpus-per-task=1 \
  --mem="${LAUNCH_MEM}" \
  --time="${LAUNCH_TIME}" \
  --job-name="alfak2_benchmark_launcher" \
  --output="${LOG_DIR}/benchmark_launcher_%j.out" \
  --error="${LOG_DIR}/benchmark_launcher_%j.err" \
  "${STAGE_RUNNER}" "launch-benchmark-after-calibration" "${COMMON_ENV}")
printf "%s\n" "${benchmark_launcher_id}" > "${CONTROL_DIR}/benchmark_launcher.jobid"

cat <<EOF
Submitted beta-calibrated GRF pipeline.

Pipeline output:       ${PIPELINE_OUTPUT_DIR}
Shared prepare job:    ${shared_prepare_id}
Calibration launchers: ${cal_launch_ids[*]}
Benchmark launcher:    ${benchmark_launcher_id}
Beta manifest:         ${MANIFEST}
Logs:                  ${LOG_DIR}

The final combined tables will be written under:
  ${COMBINED_OUTPUT_DIR}/tables
EOF
