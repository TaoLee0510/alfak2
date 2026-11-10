#!/usr/bin/env Rscript
source("validation/accuracy_round2/pooled_calibration.R")
input <- "validation/accuracy_upgrade/results/final_exact/node_matrices.rds"
manifest_file <- "validation/accuracy_upgrade/results/final_exact/test_dataset_manifest.csv"
out <- "validation/accuracy_round2/results/development_calibration_selection"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(input)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
stopifnot(identical(manifest$dataset_id, x$dataset_ids))
for (model in c("baseline", "division_v1")) {
  z <- x$models[[model]]
  data <- do.call(rbind, lapply(seq_along(x$dataset_ids), function(j)
    data.frame(landscape_id = manifest$landscape_id[j], dataset_id = x$dataset_ids[j],
      family = manifest$family[j], split = "development", original_split = "round1_revealed_test",
      distance = x$reference_support_distance[, j], pred_relative = z$pred_relative[, j],
      sd_relative = z$sd_relative[, j], truth_relative = x$truth_relative[, j])))
  selection <- select_round2_calibration(data, list(model = model), x$scenario)
  saveRDS(selection, file.path(out, paste0(model, "_selection.rds")))
  for (name in c("selected", "grid_summary", "distance_risk", "heldout_scores"))
    write.csv(selection[[name]], file.path(out, paste0(model, "_", name, ".csv")), row.names = FALSE)
  cat(model, "development-only selection:",
      "gamma_far", selection$selected$gamma_far, "kappa", selection$selected$kappa,
      "interval score", selection$selected$macro_interval_score, "\n")
}
