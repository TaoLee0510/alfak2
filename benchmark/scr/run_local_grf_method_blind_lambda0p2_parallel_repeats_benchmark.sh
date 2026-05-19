#!/usr/bin/env bash
#
# Local staged/parallel benchmark for the method-blind GRF landscape generator.
# This wrapper intentionally keeps run_grf_alfak2_vs_alfakR_benchmark.R
# unchanged: fit repeats are represented as independent repeat output dirs.
#
# Run:
#   bash benchmark/scr/run_local_grf_method_blind_lambda0p2_parallel_repeats_benchmark.sh
#
# Useful overrides:
#   FIT_REPEATS=3 N_CORES=3 N_SIM=2 METHODS=empirical MINOBS=5 \
#   bash benchmark/scr/run_local_grf_method_blind_lambda0p2_parallel_repeats_benchmark.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_SELF="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
ALFAK2_REPO="${ALFAK2_REPO:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"
ALFAKR_REPO="${ALFAKR_REPO:-$(cd "${ALFAK2_REPO}/.." && pwd)/alfakR}"

if [[ ! -d "${ALFAK2_REPO}" ]]; then
  echo "Missing alfak2 repo: ${ALFAK2_REPO}" >&2
  exit 2
fi
if [[ ! -d "${ALFAKR_REPO}" ]]; then
  echo "Missing alfakR repo: ${ALFAKR_REPO}" >&2
  echo "Set ALFAKR_REPO=/path/to/alfakR if it is not a sibling of alfak2." >&2
  exit 2
fi

OUTPUT_DIR="${OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/local_grf_method_blind_lambda0p2_full_minobs_local1_global1_parallel_repeats_pm_5e_05}"
SHARED_INPUT_DIR="${SHARED_INPUT_DIR:-${SOURCE_INPUT_DIR:-${OUTPUT_DIR}/shared_inputs}}"
PREPARE_SHARED_INPUTS="${PREPARE_SHARED_INPUTS:-true}"
FIT_REPEATS="${FIT_REPEATS:-3}"
N_CORES="${N_CORES:-3}"
SEED="${SEED:-424242}"
REPEAT_SEED_STEP="${REPEAT_SEED_STEP:-1000000}"

METHODS="${METHODS:-empirical}"
MINOBS="${MINOBS:-5}"
N_SIM="${N_SIM:-2}"
LAMBDA="${LAMBDA:-0.2}"
TIME_GAPS="${TIME_GAPS:-2}"
TIME_STARTS="${TIME_STARTS:-0}"
PM="${PM:-5e-05}"
BETA_LEVELS="${BETA_LEVELS:-5e-05}"
NBOOT="${NBOOT:-10}"
GRID_N="${GRID_N:-41}"
SAMPLE_DEPTH="${SAMPLE_DEPTH:-1000}"
N_CENTROIDS="${N_CENTROIDS:-32}"
K_DIM="${K_DIM:-22}"
TIME_MAX="${TIME_MAX:-180}"
PASSAGE_INTERVAL="${PASSAGE_INTERVAL:-45}"
ABM_POP_SIZE="${ABM_POP_SIZE:-10000}"
ABM_MAX_POP="${ABM_MAX_POP:-200000}"
ABM_CULLING_SURVIVAL="${ABM_CULLING_SURVIVAL:-0.01}"

ALFAK2_INPUT_POLICIES="${ALFAK2_INPUT_POLICIES:-full,minobs_matched}"
ALFAK2_INPUT_DEPTH="${ALFAK2_INPUT_DEPTH:-effective}"
ALFAK2_OBSERVATION_MODEL="${ALFAK2_OBSERVATION_MODEL:-dirichlet_multinomial}"
ALFAK2_DM_CONCENTRATION="${ALFAK2_DM_CONCENTRATION:-50}"
ALFAK2_EFFECTIVE_DEPTH_MODE="${ALFAK2_EFFECTIVE_DEPTH_MODE:-min}"
ALFAK2_LEGACY_WEIGHT="${ALFAK2_LEGACY_WEIGHT:-directly_informed}"
ALFAK2_LAMBDA_L_GRID="${ALFAK2_LAMBDA_L_GRID:-0.2}"
ALFAK2_LAMBDA_E_GRID="${ALFAK2_LAMBDA_E_GRID:-1}"
ALFAK2_SIGMA_OBS_GRID="${ALFAK2_SIGMA_OBS_GRID:-0.02}"
ALFAK2_GRAPH_EDGE_WEIGHT="${ALFAK2_GRAPH_EDGE_WEIGHT:-mutation}"
ALFAK2_ANCHOR_COUNT_REFERENCE="${ALFAK2_ANCHOR_COUNT_REFERENCE:-minobs}"
ALFAK2_ANCHOR_COUNT_POWER="${ALFAK2_ANCHOR_COUNT_POWER:-1}"
ALFAK2_LOCAL_SHELL_DEPTH="${ALFAK2_LOCAL_SHELL_DEPTH:-1}"
ALFAK2_GLOBAL_EXTRA_SHELL="${ALFAK2_GLOBAL_EXTRA_SHELL:-1}"
ALFAK2_MAX_NODES="${ALFAK2_MAX_NODES:-150000}"

GRF_CENTROID_MODE="${GRF_CENTROID_MODE:-method_blind}"
GRF_CENTROID_MIN_CN="${GRF_CENTROID_MIN_CN:-0}"
GRF_CENTROID_MAX_CN="${GRF_CENTROID_MAX_CN:-4}"

FORCE_REFIT="${FORCE_REFIT:-false}"
FORCE_SIM="${FORCE_SIM:-false}"
REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE:-false}"
RECOMPILE_DLL="${RECOMPILE_DLL:-false}"
WRITE_NODE_ACCURACY="${WRITE_NODE_ACCURACY:-true}"
DRY_RUN="${DRY_RUN:-false}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

repeat_output_dir() {
  local repeat_index="$1"
  printf "%s/repeat_%02d" "${OUTPUT_DIR}" "${repeat_index}"
}

repeat_seed() {
  local repeat_index="$1"
  printf "%s" "$((SEED + (repeat_index - 1) * REPEAT_SEED_STEP))"
}

build_benchmark_args() {
  local output_dir="$1"
  local seed_value="$2"
  local source_input_dir="${3:-}"
  BENCHMARK_ARGS=(
    "--alfak2-repo=${ALFAK2_REPO}"
    "--alfakR-repo=${ALFAKR_REPO}"
    "--output-dir=${output_dir}"
    "--n-cores=1"
    "--seed=${seed_value}"
    "--methods=${METHODS}"
    "--minobs=${MINOBS}"
    "--n-sim=${N_SIM}"
    "--lambdas=${LAMBDA}"
    "--time-gaps=${TIME_GAPS}"
    "--time-starts=${TIME_STARTS}"
    "--pm=${PM}"
    "--beta-levels=${BETA_LEVELS}"
    "--nboot=${NBOOT}"
    "--grid-n=${GRID_N}"
    "--sample-depth=${SAMPLE_DEPTH}"
    "--n-centroids=${N_CENTROIDS}"
    "--k-dim=${K_DIM}"
    "--time-max=${TIME_MAX}"
    "--passage-interval=${PASSAGE_INTERVAL}"
    "--abm-pop-size=${ABM_POP_SIZE}"
    "--abm-max-pop=${ABM_MAX_POP}"
    "--abm-culling-survival=${ABM_CULLING_SURVIVAL}"
    "--grf-centroid-mode=${GRF_CENTROID_MODE}"
    "--grf-centroid-min-cn=${GRF_CENTROID_MIN_CN}"
    "--grf-centroid-max-cn=${GRF_CENTROID_MAX_CN}"
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
    "--recompile-dll=${RECOMPILE_DLL}"
    "--write-node-accuracy=${WRITE_NODE_ACCURACY}"
  )
  if [[ -n "${source_input_dir}" ]]; then
    BENCHMARK_ARGS+=("--source-input-dir=${source_input_dir}")
  fi
  if [[ -n "${EXTRA_ARGS:-}" ]]; then
    # shellcheck disable=SC2206
    EXTRA_ARGS_ARRAY=(${EXTRA_ARGS})
    BENCHMARK_ARGS+=("${EXTRA_ARGS_ARRAY[@]}")
  fi
}

run_benchmark_custom() {
  local mode="$1"
  local output_dir="$2"
  local seed_value="$3"
  local source_input_dir="${4:-}"
  shift 4 || true
  build_benchmark_args "${output_dir}" "${seed_value}" "${source_input_dir}"
  if [[ "${DRY_RUN}" == "true" ]]; then
    printf "Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R --mode=%q" "${mode}"
    printf " %q" "${BENCHMARK_ARGS[@]}"
    if (($#)); then
      printf " %q" "$@"
    fi
    printf "\n"
  else
    Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R "--mode=${mode}" "${BENCHMARK_ARGS[@]}" "$@"
  fi
}

run_benchmark_mode() {
  local mode="$1"
  local repeat_index="$2"
  shift 2
  run_benchmark_custom "${mode}" "$(repeat_output_dir "${repeat_index}")" "$(repeat_seed "${repeat_index}")" "${SHARED_INPUT_DIR}" "$@"
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

run_common_node_summary() {
  local repeat_dir="$1"
  Rscript benchmark/scr/summarize_grf_common_nodes_by_alfakR_scope.R "${repeat_dir}"
}

combine_repeat_tables() {
  Rscript -e '
    args <- commandArgs(TRUE)
    root <- args[[1]]
    fit_repeats <- as.integer(args[[2]])
    if (!requireNamespace("data.table", quietly = TRUE)) quit(status = 0)
    tables_dir <- file.path(root, "tables")
    dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
    read_repeat <- function(rep_idx, file_name) {
      path <- file.path(root, sprintf("repeat_%02d", rep_idx), "tables", file_name)
      if (!file.exists(path)) return(NULL)
      x <- data.table::fread(path, sep = "\t", header = TRUE, showProgress = FALSE)
      x[, fit_repeat := rep_idx]
      x[]
    }
    write_combined <- function(file_name, out_name = paste0("combined_", file_name)) {
      xs <- lapply(seq_len(fit_repeats), read_repeat, file_name = file_name)
      xs <- Filter(Negate(is.null), xs)
      if (!length(xs)) return(invisible(NULL))
      data.table::fwrite(data.table::rbindlist(xs, fill = TRUE), file.path(tables_dir, out_name), sep = "\t")
    }
    write_combined("fit_results.tsv")
    write_combined("summary_by_lambda_time_minobs_method.tsv")
    write_combined("common_node_condition_metrics_by_alfakR_scope.tsv")
    write_combined("common_node_performance_summary_by_alfakR_scope.tsv")
    condition_path <- file.path(tables_dir, "combined_common_node_condition_metrics_by_alfakR_scope.tsv")
    if (file.exists(condition_path)) {
      x <- data.table::fread(condition_path, sep = "\t", header = TRUE, showProgress = FALSE)
      x[, alfakR_group := ifelse(alfakR_method == "alfakR_none", "alfakR none", "alfakR NN-prior")]
      y <- x[support_scope %in% c("fq", "NN", "whole"),
             .(conditions = .N,
               median_delta_centered_rmse = median(delta_centered_rmse, na.rm = TRUE),
               alfak2_win_rate = mean(delta_centered_rmse < 0, na.rm = TRUE)),
             by = .(support_scope, minobs, alfakR_group)]
      data.table::setorder(y, support_scope, minobs, alfakR_group)
      data.table::fwrite(y, file.path(tables_dir, "combined_common_node_quick_summary.tsv"), sep = "\t")
      print(y)
    }
  ' "${OUTPUT_DIR}" "${FIT_REPEATS}"
}

run_fit_task_from_queue() {
  local repeat_index="$1"
  local task_index="$2"
  cd "${ALFAK2_REPO}"
  run_benchmark_mode fit-task "${repeat_index}" "--task-index=${task_index}"
}

case "${1:-}" in
  __fit_task)
    run_fit_task_from_queue "${2:-}" "${3:-}"
    exit
    ;;
esac

cd "${ALFAK2_REPO}"
mkdir -p "${OUTPUT_DIR}"

echo "Local method-blind GRF staged parallel repeat benchmark"
echo "  alfak2:       ${ALFAK2_REPO}"
echo "  alfakR:       ${ALFAKR_REPO}"
echo "  output root:  ${OUTPUT_DIR}"
echo "  shared input: ${SHARED_INPUT_DIR}"
echo "  fit_repeats:  ${FIT_REPEATS}"
echo "  n_cores:      ${N_CORES}"
echo "  lambda:       ${LAMBDA}"
echo "  pm/beta:      ${PM} / ${BETA_LEVELS}"
echo "  minobs:       ${MINOBS}"
echo "  n_sim:        ${N_SIM}"
echo "  methods:      ${METHODS}"
echo "  GRF mode:     ${GRF_CENTROID_MODE}, CN ${GRF_CENTROID_MIN_CN}-${GRF_CENTROID_MAX_CN}"

if [[ "${DRY_RUN}" == "true" ]]; then
  echo "DRY-RUN shared input prepare command:"
  run_benchmark_custom prepare "${SHARED_INPUT_DIR}" "${SEED}" ""
  echo "DRY-RUN prepare commands:"
  for repeat_index in $(seq 1 "${FIT_REPEATS}"); do
    run_benchmark_mode prepare "${repeat_index}"
  done
  echo "DRY-RUN fit queue command:"
  printf "xargs -P %q -n 2 bash %q __fit_task < %q\n" "${N_CORES}" "${SCRIPT_SELF}" "${OUTPUT_DIR}/control/fit_task_queue.tsv"
  echo "DRY-RUN summarize commands:"
  for repeat_index in $(seq 1 "${FIT_REPEATS}"); do
    run_benchmark_mode summarize "${repeat_index}"
  done
  exit 0
fi

echo "Stage 1/4: preparing shared inputs and repeat task tables."
if [[ "${PREPARE_SHARED_INPUTS}" == "true" || ! -f "${SHARED_INPUT_DIR}/tables/input_table.tsv" ]]; then
  mkdir -p "${SHARED_INPUT_DIR}"
  run_benchmark_custom prepare "${SHARED_INPUT_DIR}" "${SEED}" ""
else
  echo "Using existing shared input table: ${SHARED_INPUT_DIR}/tables/input_table.tsv"
fi

for repeat_index in $(seq 1 "${FIT_REPEATS}"); do
  mkdir -p "$(repeat_output_dir "${repeat_index}")"
  run_benchmark_mode prepare "${repeat_index}"
done

CONTROL_DIR="${OUTPUT_DIR}/control"
QUEUE_FILE="${CONTROL_DIR}/fit_task_queue.tsv"
mkdir -p "${CONTROL_DIR}"
: > "${QUEUE_FILE}"

total_tasks=0
for repeat_index in $(seq 1 "${FIT_REPEATS}"); do
  task_table="$(repeat_output_dir "${repeat_index}")/tables/task_table.tsv"
  n_tasks="$(task_count_from_table "${task_table}")"
  total_tasks=$((total_tasks + n_tasks))
  for task_index in $(seq 1 "${n_tasks}"); do
    printf "%s\t%s\n" "${repeat_index}" "${task_index}" >> "${QUEUE_FILE}"
  done
done

echo "Stage 2/4: running ${total_tasks} fit tasks with xargs -P ${N_CORES}."
xargs -P "${N_CORES}" -n 2 bash "${SCRIPT_SELF}" __fit_task < "${QUEUE_FILE}"

echo "Stage 3/4: summarizing each repeat."
for repeat_index in $(seq 1 "${FIT_REPEATS}"); do
  run_benchmark_mode summarize "${repeat_index}"
  run_common_node_summary "$(repeat_output_dir "${repeat_index}")"
done

echo "Stage 4/4: combining repeat-level tables."
combine_repeat_tables

echo "Wrote repeat benchmark outputs under: ${OUTPUT_DIR}"
