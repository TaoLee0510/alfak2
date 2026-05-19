#!/usr/bin/env bash
#
# Local mini benchmark for the method-blind GRF landscape generator.
# Assumes alfak2 and alfakR are sibling directories by default:
#   .../GitHub/alfak2
#   .../GitHub/alfakR
#
# Run:
#   bash benchmark/scr/run_local_grf_method_blind_lambda0p2_benchmark.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

OUTPUT_DIR="${OUTPUT_DIR:-${ALFAK2_REPO}/benchmark/results/local_grf_method_blind_lambda0p2_pm_5e_05}"

METHODS="${METHODS:-none,empirical,empirical_censored,empirical_censored_weighted,empirical_two_step}"
MINOBS="${MINOBS:-5,10,20}"
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

ALFAK2_INPUT_POLICIES="${ALFAK2_INPUT_POLICIES:-minobs_matched}"
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
ALFAK2_LOCAL_SHELL_DEPTH="${ALFAK2_LOCAL_SHELL_DEPTH:-0}"
ALFAK2_GLOBAL_EXTRA_SHELL="${ALFAK2_GLOBAL_EXTRA_SHELL:-1}"
ALFAK2_MAX_NODES="${ALFAK2_MAX_NODES:-150000}"

GRF_CENTROID_MODE="${GRF_CENTROID_MODE:-method_blind}"
GRF_CENTROID_MIN_CN="${GRF_CENTROID_MIN_CN:-0}"
GRF_CENTROID_MAX_CN="${GRF_CENTROID_MAX_CN:-4}"

FORCE_REFIT="${FORCE_REFIT:-false}"
FORCE_SIM="${FORCE_SIM:-false}"
REUSE_DIRTY_CACHE="${REUSE_DIRTY_CACHE:-false}"
RECOMPILE_DLL="${RECOMPILE_DLL:-false}"
DRY_RUN="${DRY_RUN:-false}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "${ALFAK2_REPO}"
mkdir -p "${OUTPUT_DIR}"

COMMON_ARGS=(
  "--mode=all"
  "--alfak2-repo=${ALFAK2_REPO}"
  "--alfakR-repo=${ALFAKR_REPO}"
  "--output-dir=${OUTPUT_DIR}"
  "--n-cores=1"
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
)

if [[ -n "${EXTRA_ARGS:-}" ]]; then
  # shellcheck disable=SC2206
  EXTRA_ARGS_ARRAY=(${EXTRA_ARGS})
  COMMON_ARGS+=("${EXTRA_ARGS_ARRAY[@]}")
fi

echo "Local method-blind GRF mini benchmark"
echo "  alfak2:       ${ALFAK2_REPO}"
echo "  alfakR:       ${ALFAKR_REPO}"
echo "  output:       ${OUTPUT_DIR}"
echo "  lambda:       ${LAMBDA}"
echo "  pm/beta:      ${PM} / ${BETA_LEVELS}"
echo "  minobs:       ${MINOBS}"
echo "  n_sim:        ${N_SIM}"
echo "  methods:      ${METHODS}"
echo "  GRF mode:     ${GRF_CENTROID_MODE}, CN ${GRF_CENTROID_MIN_CN}-${GRF_CENTROID_MAX_CN}"

if [[ "${DRY_RUN}" == "true" ]]; then
  printf "Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R"
  printf " %q" "${COMMON_ARGS[@]}"
  printf "\n"
  exit 0
fi

Rscript benchmark/scr/run_grf_alfak2_vs_alfakR_benchmark.R "${COMMON_ARGS[@]}"

Rscript benchmark/scr/summarize_grf_common_nodes_by_alfakR_scope.R "${OUTPUT_DIR}"

Rscript -e '
  tables_dir <- file.path(commandArgs(TRUE)[[1]], "tables")
  path <- file.path(tables_dir, "common_node_condition_metrics_by_alfakR_scope.tsv")
  if (!file.exists(path)) quit(status = 0)
  x <- data.table::fread(path)
  x[, alfakR_group := ifelse(alfakR_method == "alfakR_none", "alfakR none", "alfakR NN-prior")]
  y <- x[support_scope %in% c("fq", "NN", "whole"),
         .(conditions = .N,
           median_delta_centered_rmse = median(delta_centered_rmse, na.rm = TRUE),
           alfak2_win_rate = mean(delta_centered_rmse < 0, na.rm = TRUE)),
         by = .(support_scope, minobs, alfakR_group)]
  data.table::setorder(y, support_scope, minobs, alfakR_group)
  print(y)
' "${OUTPUT_DIR}"

echo "Wrote benchmark tables under: ${OUTPUT_DIR}/tables"
