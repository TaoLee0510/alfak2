#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R"))
input <- round2_arg("input", file.path(script_dir, "results", "development81_comparison"))
report <- round2_arg("report", file.path(script_dir, "development_stage_review.md"))
overall <- utils::read.csv(file.path(input, "overall_metrics.csv"), stringsAsFactors = FALSE)
family <- utils::read.csv(file.path(input, "family_metrics.csv"), stringsAsFactors = FALSE)
land <- utils::read.csv(file.path(input, "landscape_metrics.csv"), stringsAsFactors = FALSE)
availability <- utils::read.csv(file.path(input, "availability_gate.csv"), stringsAsFactors = FALSE)
distances <- 0:3
metrics <- c("rmse", "spatially_centered_rmse", "edge_gradient_rmse", "interval_score95")
point_metrics <- metrics[1:3]
candidates <- c("A", "B")
if (!all(c("v1_joint", candidates) %in% overall$model)) stop("Complete fixed V1/A/B comparison required")
required <- expand.grid(model = c("v1_joint", candidates), distance = distances)
if (!all(paste(required$model, required$distance) %in% paste(overall$model, overall$distance)))
  stop("Missing predeclared native-distance strata")

gains <- safety <- pairs <- offset <- decisions <- list()
for (candidate in candidates) {
  for (d in distances) for (metric in metrics) {
    reference <- overall[[metric]][overall$model == "v1_joint" & overall$distance == d]
    value <- overall[[metric]][overall$model == candidate & overall$distance == d]
    if (length(reference) != 1L || length(value) != 1L || any(!is.finite(c(reference, value))))
      stop("Incomplete development risk grid")
    gains[[length(gains) + 1L]] <- data.frame(candidate_id = candidate, distance = d,
      metric = metric, reference_risk = reference, candidate_risk = value,
      ratio = value / reference, percent_change = 100 * (value / reference - 1),
      observed_improvement = value < reference,
      meets_5percent_planning_target = value <= .95 * reference,
      role = if (metric == "interval_score95") "raw_uncalibrated_diagnostic" else "development_screen")
  }
  for (f in unique(family$family)) for (d in distances) for (metric in point_metrics) {
    reference <- family[[metric]][family$model == "v1_joint" & family$family == f & family$distance == d]
    value <- family[[metric]][family$model == candidate & family$family == f & family$distance == d]
    if (length(reference) != 1L || length(value) != 1L || any(!is.finite(c(reference, value))))
      stop("Incomplete family-by-distance shape safety grid")
    safety[[length(safety) + 1L]] <- data.frame(candidate_id = candidate, family = f,
      distance = d, metric = metric, reference_risk = reference, candidate_risk = value,
      ratio = value / reference, observed_within_5percent_margin = value <= 1.05 * reference,
      inference = "descriptive_screen_not_formal_noninferiority")
  }
  for (metric in point_metrics) {
    reference <- land[land$model == "v1_joint" & land$distance %in% distances,
      c("family", "landscape_id", "distance", metric)]
    value <- land[land$model == candidate & land$distance %in% distances,
      c("family", "landscape_id", "distance", metric)]
    names(reference)[4L] <- "reference_loss"; names(value)[4L] <- "candidate_loss"
    paired <- merge(reference, value, by = c("family", "landscape_id", "distance"), all = TRUE)
    if (nrow(paired) != nrow(reference) || any(!is.finite(paired$reference_loss)) || any(!is.finite(paired$candidate_loss)))
      stop("Unpaired landscape losses")
    paired$candidate_id <- candidate; paired$metric <- metric
    paired$role <- "development_only_point_shape_screen"
    paired$oof_verified <- FALSE
    paired$why_not_joint_oof <- "No frozen single winner or nested interval-calibration OOF loss vector; do not pass to joint-power certification"
    pairs[[length(pairs) + 1L]] <- paired
  }
  reference <- overall[overall$model == "v1_joint" & overall$distance == -1L, ]
  value <- overall[overall$model == candidate & overall$distance == -1L, ]
  for (x in list(reference, value))
    if (abs(x$mse - x$spatially_centered_mse - x$global_offset_squared) > 1e-10)
      stop("Whole-domain MSE shape/offset decomposition failed")
  mse_gain <- reference$mse - value$mse
  offset[[length(offset) + 1L]] <- data.frame(candidate_id = candidate,
    reference_mse = reference$mse, candidate_mse = value$mse,
    reference_shape_mse = reference$spatially_centered_mse,
    candidate_shape_mse = value$spatially_centered_mse,
    reference_offset_squared = reference$global_offset_squared,
    candidate_offset_squared = value$global_offset_squared,
    mse_reduction = mse_gain,
    shape_mse_reduction = reference$spatially_centered_mse - value$spatially_centered_mse,
    offset_mse_reduction = reference$global_offset_squared - value$global_offset_squared,
    offset_fraction_of_mse_reduction = if (mse_gain > 0) (reference$global_offset_squared - value$global_offset_squared) / mse_gain else NA_real_)
}
gains <- do.call(rbind, gains); safety <- do.call(rbind, safety)
for (candidate in candidates) {
  x <- gains[gains$candidate_id == candidate, ]; s <- safety[safety$candidate_id == candidate, ]
  a <- availability[availability$model == candidate, ]
  decisions[[length(decisions) + 1L]] <- data.frame(candidate_id = candidate,
    all_point_shape_predictions_available = all(a$all_fits_returned & a$all_predictions_available),
    all_distance_point_5percent_target = all(x$meets_5percent_planning_target[x$metric == "rmse"]),
    all_distance_shape_improved = all(x$observed_improvement[x$metric == "spatially_centered_rmse"]),
    all_distance_gradient_improved = all(x$observed_improvement[x$metric == "edge_gradient_rmse"]),
    all_descriptive_family_shape_safety = all(s$observed_within_5percent_margin),
    nested_calibration_performed = FALSE, joint_power_validated = FALSE,
    formal_promotion_test_performed = FALSE, new_calibration_generated = FALSE,
    new_test_generated = FALSE, freeze_single_winner = FALSE,
    stage_decision = "DEFER_NO_CREDIBLE_ALL_DISTANCE_POINT_SHAPE_WINNER")
}
decisions <- do.call(rbind, decisions); offset <- do.call(rbind, offset)
round2_write_csv(gains, file.path(input, "stage_overall_risk_changes.csv"))
round2_write_csv(safety, file.path(input, "stage_family_shape_safety.csv"))
round2_write_csv(do.call(rbind, pairs), file.path(input, "development_point_shape_paired_losses.csv"))
round2_write_csv(offset, file.path(input, "stage_shape_offset_decomposition.csv"))
round2_write_csv(decisions, file.path(input, "stage_decision.csv"))
sources <- file.path(input, c("overall_metrics.csv", "family_metrics.csv", "landscape_metrics.csv", "availability_gate.csv"))
round2_write_csv(data.frame(file = basename(sources), sha256 = unname(round2_sha256(sources))),
  file.path(input, "stage_inputs_sha256.csv"))

lines <- c("# Round 2 development stage decision", "",
  "The 18-landscape 81-node development screen does not identify a credible single winner for the requested all-distance point, shape and uncertainty improvement. No new calibration or test dataset is generated. This is a development stopping decision, not a formal rejection of noninferiority or an independent-test result.", "",
  "All models use the same unchanged counts, 81 target nodes and observed-reference contrast. Six independent old exact landscapes per family are used, one A1/S1 dataset per landscape. Old calibration/revealed-test data are explicitly development material. A/B each completed 18/18 fits with trusted TMB covariance. Frozen V1 summaries are reused with exact dataset/reference/truth checks.", "",
  "## Point and shape screening", "",
  "Percent changes are candidate/V1 minus one; negative is better. These are equal-family means of landscape risks. The 5% target is the prospective planning effect, not a claim that a superiority test must establish a full 5% margin.", "",
  "| Candidate | Distance | RMSE change | Shape RMSE change | Edge RMSE change | Raw interval-score change |",
  "|---|---:|---:|---:|---:|---:|")
for (candidate in candidates) for (d in distances) {
  x <- gains[gains$candidate_id == candidate & gains$distance == d, ]
  values <- x$percent_change[match(metrics, x$metric)]
  lines <- c(lines, sprintf("| %s | %d | %+.3f%% | %+.3f%% | %+.3f%% | %+.3f%% |", candidate, d, values[1], values[2], values[3], values[4]))
}
lines <- c(lines, "", "A changes point/shape risk only at numerical scale; its raw interval improvement does not establish improved landscape reconstruction. B lowers global reference-contrast error, while shape/gradient improvement is absent across several distances.", "",
  "## Descriptive family safety flags", "",
  "The following point/shape cells exceed the 5% family margin in their observed mean. They are retained as safety warnings; no p-value or formal noninferiority conclusion is inferred from this small development sample.", "",
  "| Candidate | Family | Distance | Metric | Ratio to V1 |", "|---|---|---:|---|---:|")
flag <- safety[!safety$observed_within_5percent_margin, ]
if (!nrow(flag)) lines <- c(lines, "| None | | | | |") else for (i in seq_len(nrow(flag))) {
  x <- flag[i, ]; lines <- c(lines, sprintf("| %s | %s | %d | %s | %.5f |", x$candidate_id, x$family, x$distance, x$metric, x$ratio))
}
lines <- c(lines, "", "## Reference offset versus shape", "",
  "For the complete target domain, MSE equals globally centered shape MSE plus the squared mean error. No truth-fitted slope, shift or shell centering is applied to an estimate.", "",
  "| Candidate | Total MSE change | Shape MSE change | Offset-squared change | Offset fraction of MSE reduction |", "|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(offset))) {
  x <- offset[i, ]; lines <- c(lines, sprintf("| %s | %+.9f | %+.9f | %+.9f | %.2f%% |", x$candidate_id,
    -x$mse_reduction, -x$shape_mse_reduction, -x$offset_mse_reduction, 100 * x$offset_fraction_of_mse_reduction))
}
lines <- c(lines, "", "B's common-offset reduction exceeds its total MSE gain because its centered shape MSE increases. Wider or better calibrated intervals cannot resolve that point/shape failure.", "",
  "## Calibration, power and sealed validation", "",
  "No final all-distance calibrated interval comparison or joint power estimate is claimed. Ordinary landscape-LOO scores used to choose gamma/kappa are development selection scores. An actual joint OOF power vector would require nested outer-landscape exclusion for that choice, plus complete RMSE, interval-score and shape/gradient gates. Because the prerequisite point/shape stage is not credible, this additional computation is deferred. The point/shape paired loss CSV is explicitly marked oof_verified=FALSE to prevent accidental use as the complete joint OOF vector.", "",
  "The primary eventual promotion remains one intersection-union conjunction at one-sided alpha 0.05, with all predeclared distance benefits and family-by-distance 5% noninferiority. Simultaneous/Holm individual-cell claims remain separate. There is no Bonferroni-based assertion that the task is impossible, and no enlargement of test N to chase significance. Provisional N=36 and the reserved cap N=60 remain uncertified until a suitable single winner and credible OOF loss vectors exist.", "",
  "New domain-calibration/test seeds remain unused. The 64-node generator is prepared and checked, including honest actual-distance coverage and the separate no-information d5/d6 fixture; no candidate accuracy claim in that domain is made from generator preparation alone.", "",
  "Full model/interval/family tables are in [the development comparison](development81_review.md). The [protocol](PROTOCOL.md) and machine-readable stage CSVs preserve the stop conditions. No prior weight, nullspace scale, truth generator or test seed is changed in response to this screen.", "",
  paste("Generated:", format(Sys.time(), tz = "America/New_York", usetz = TRUE)))
writeLines(lines, report)
cat(normalizePath(report), "\n")
