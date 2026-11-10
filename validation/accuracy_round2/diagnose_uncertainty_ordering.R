#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
source(file.path(script_dir, "common.R"))
input <- round2_arg("matrices", file.path(script_dir, "..", "accuracy_upgrade", "results", "final_exact", "node_matrices.rds"))
manifest_file <- round2_arg("manifest", file.path(script_dir, "..", "accuracy_upgrade", "results", "final_exact", "test_dataset_manifest.csv"))
output <- round2_arg("output", file.path(script_dir, "results", "old_data_diagnostics"))
x <- readRDS(input); manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
mapping <- c(original_baseline = "baseline", v1_joint = "division_v1", endpoint_ablation = "endpoint_v1")
correlation <- function(a, b, method) {
  valid <- is.finite(a) & is.finite(b)
  if (sum(valid) < 5L || stats::sd(a[valid]) <= 1e-12 || stats::sd(b[valid]) <= 1e-12) return(NA_real_)
  suppressWarnings(stats::cor(a[valid], b[valid], method = method))
}
rows <- list()
for (model in names(mapping)) {
  z <- x$models[[mapping[[model]]]]
  error <- z$pred_relative - x$truth_relative
  for (landscape in unique(manifest$landscape_id)) {
    columns <- which(manifest$landscape_id == landscape)
    repeat_mean_error <- rowMeans(error[, columns, drop = FALSE])
    for (column in columns) for (d in sort(unique(x$reference_support_distance[, column]))) {
      ii <- which(x$reference_support_distance[, column] == d)
      e <- error[ii, column]; sd <- z$sd_relative[ii, column]
      centered_spatial_sqerror <- (e - mean(e))^2
      centered_repeat_sqerror <- (e - repeat_mean_error[ii])^2
      rows[[length(rows) + 1L]] <- data.frame(model = model, family = manifest$family[column],
        landscape_id = landscape, dataset_id = manifest$dataset_id[column], distance = d, n_target = length(ii),
        sd_median = stats::median(sd), sd_iqr_over_median = stats::IQR(sd) / stats::median(sd),
        sd_p90_over_p10 = round2_quantile(sd, 0.9) / round2_quantile(sd, 0.1),
        spearman_abs_error = correlation(sd, abs(e), "spearman"),
        spearman_centered_spatial_squared_error = correlation(sd, centered_spatial_sqerror, "spearman"),
        spearman_centered_repeat_squared_error = correlation(sd, centered_repeat_sqerror, "spearman"),
        pearson_abs_error = correlation(sd, abs(e), "pearson"),
        pearson_centered_spatial_squared_error = correlation(sd, centered_spatial_sqerror, "pearson"),
        pearson_centered_repeat_squared_error = correlation(sd, centered_repeat_sqerror, "pearson"))
    }
  }
}
case <- do.call(rbind, rows)
measure <- setdiff(names(case), c("model", "family", "landscape_id", "dataset_id", "distance"))
land <- aggregate(case[, measure], case[, c("model", "family", "landscape_id", "distance")], round2_mean)
family <- aggregate(land[, measure], land[, c("model", "family", "distance")], round2_mean)
round2_write_csv(case, file.path(output, "uncertainty_ordering_case.csv"))
round2_write_csv(land, file.path(output, "uncertainty_ordering_landscape.csv"))
round2_write_csv(family, file.path(output, "uncertainty_ordering_family.csv"))
print(subset(family, model == "v1_joint", select = c("family", "distance", "sd_iqr_over_median", "spearman_abs_error", "spearman_centered_spatial_squared_error", "spearman_centered_repeat_squared_error")), row.names = FALSE)
