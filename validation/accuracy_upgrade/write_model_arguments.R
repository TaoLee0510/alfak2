#!/usr/bin/env Rscript
# Explicit development candidates; no fitting or truth access.
args <- commandArgs(TRUE)
output <- if (length(args)) args[[1L]] else tempfile("alfak2_model_arguments_")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
models <- list(
  ablation_laplace = list(pipeline = "two_layer", forward_model = "selection_transport",
                          inference_method = "laplace", extrapolation_control = list(compute_covariance = TRUE)),
  joint_endpoint = list(pipeline = "joint_count", forward_model = "selection_transport",
                        inference_method = "laplace", sdreport_mode = "all_f_current"),
  joint_division = list(pipeline = "joint_count", forward_model = "division_process",
                        inference_method = "laplace"))
for (scale in c(0.05, 0.2, 0.5, 1, 2)) {
  tag <- gsub(".", "p", as.character(scale), fixed = TRUE)
  models[[paste0("joint_endpoint_s", tag)]] <- list(
    pipeline = "joint_count", forward_model = "selection_transport", inference_method = "laplace",
    sdreport_mode = "all_f_current", joint_control = list(baseline_sd = 0.5,
      main_sd = 0.35 * scale, curvature_sd = 0.2 * scale,
      interaction_sd = 0.2 * scale, residual_sd = 0.25 * scale, edge_sd = 0.5 * scale))
}
for (basis in c(0.1, 0.2, 0.35, 0.5)) for (residual in c(0.2, 0.5, 1)) {
  tag <- function(x) gsub(".", "p", as.character(x), fixed = TRUE)
  models[[paste0("joint_endpoint_b", tag(basis), "_r", tag(residual))]] <- list(
    pipeline = "joint_count", forward_model = "selection_transport", inference_method = "laplace",
    sdreport_mode = "all_f_current", joint_control = list(baseline_sd = 0.5,
      main_sd = 0.35 * basis, curvature_sd = 0.2 * basis, interaction_sd = 0.2 * basis,
      residual_sd = 0.25 * residual, edge_sd = 0.5 * residual))
}
for (model in names(models)) saveRDS(models[[model]], file.path(output, paste0(model, ".rds")))
manifest <- data.frame(model = names(models),
                       fit_arguments_file = normalizePath(file.path(output, paste0(names(models), ".rds"))))
utils::write.csv(manifest, file.path(output, "model_arguments_manifest.csv"), row.names = FALSE)
cat(normalizePath(output), "\n")
