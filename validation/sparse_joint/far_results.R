# Figures and tables for distances 0 to 4 (Section 5.4 and Supplement S5 of
# doc/alfak2_method.tex).
#
#   Rscript validation/sparse_joint/far_results.R [OUT_DIR]
#
# Distances 0-2 come from the scoreboard (metrics_<split>.csv), distances 3-4
# from far_distance.R (sweeps/far/metrics_far_<split>.csv). Both are centred on
# the same reference karyotypes. Values are averaged over the scenarios of each
# landscape and then over the three landscapes of each wavelength; error
# ratios divide the mean error by the mean error of the null predictor.

suppressMessages(library(ggplot2))
args <- commandArgs(TRUE)
out_dir <- if (length(args)) args[[1L]] else "doc/generated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results <- "benchmark/results/sparse_two_timepoint"

labels <- c(
  "zero" = "Null predictor",
  "alfak2_two_layer_fixed" = "Graph Gaussian",
  "alfak2_two_layer_nearfield_gmrf" = "Truncated GMRF",
  "alfak2_two_layer_edge_effect" = "Edge effect",
  "alfak2_two_layer_path_ensemble" = "Path ensemble",
  "alfak2_two_layer_polynomial" = "Local polynomial",
  "alfak2_two_layer_nngp" = "Nearest-neighbour GP",
  "alfak2_two_layer_delta_tree" = "Delta tree",
  "alfak2_two_layer_trend_filtering" = "Trend filtering",
  "alfakR_minobs20" = "alfakR (minobs = 20)",
  "alfak2_sparse_joint" = "Sparse joint, independent",
  "alfak2_sparse_joint_balance" = "Sparse joint, balance (default)"
)
waves <- c("w0.4", "w0.8", "w1.2", "w1.6")
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
cols <- c("rmse", "rmse_zero", "bias", "spearman", "coverage95", "prediction_fraction")

load <- function(split) {
  near <- utils::read.csv(file.path(results, sprintf("metrics_%s.csv", split)), stringsAsFactors = FALSE)
  far <- utils::read.csv(file.path(results, "sweeps", "far", sprintf("metrics_far_%s.csv", split)), stringsAsFactors = FALSE)
  near <- near[near$stratum %in% c("d0", "d1", "d2"), ]
  far <- far[far$stratum %in% c("d3", "d4"), ]
  keep <- c("dataset_id", "method", "landscape_id", "wavelength_class", "scenario", "role", "scoreboard", "stratum", cols)
  m <- rbind(near[keep], far[keep])
  m <- m[m$role == "core" & m$scoreboard & m$method %in% names(labels), ]
  cell <- stats::aggregate(m[cols], m[c("method", "wavelength_class", "landscape_id", "scenario", "stratum")], mean_na)
  land <- stats::aggregate(cell[cols], cell[c("method", "wavelength_class", "landscape_id", "stratum")], mean_na)
  wl <- stats::aggregate(land[cols], land[c("method", "wavelength_class", "stratum")], mean_na)
  wl$ratio <- wl$rmse / wl$rmse_zero
  wl$distance <- as.integer(sub("^d", "", wl$stratum))
  wl$label <- labels[wl$method]
  wl
}
held <- load("heldout")
dev <- load("dev")

## Figures: distance 0 to 4, four main methods ---------------------------------
palette <- c("Sparse joint, balance (default)" = "#2a78d6", "Sparse joint, independent" = "#eb6834",
             "Graph Gaussian" = "#1baf7a", "alfakR (minobs = 20)" = "#eda100")
shapes <- c("Sparse joint, balance (default)" = 16, "Sparse joint, independent" = 17,
            "Graph Gaussian" = 15, "alfakR (minobs = 20)" = 18)
cap <- 8
plotted <- held[held$label %in% names(palette), ]
plotted$label <- factor(plotted$label, levels = names(palette))
plotted$lambda <- factor(sub("^w", "lambda == ", plotted$wavelength_class), levels = sub("^w", "lambda == ", waves))
plotted$ratio_shown <- pmin(plotted$ratio, cap)
plotted$above <- ifelse(plotted$ratio > cap, sprintf("%.0f", plotted$ratio), "")
theme_doc <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"), strip.text = element_text(colour = "#1f1f1e"),
        legend.position = "top", legend.title = element_blank(), panel.spacing = unit(8, "pt"))
# Shaded band over the distances added by far_distance.R. On the log scale the
# band must stay within the axis limits, since -Inf has no logarithm.
far_band <- function(ymin = -Inf, ymax = Inf) annotate("rect", xmin = 2.5, xmax = 4.4, ymin = ymin, ymax = ymax, fill = "#f2f1ec")
ratio_plot <- ggplot(plotted, aes(distance, ratio_shown, colour = label, shape = label)) +
  far_band(0.2, cap * 1.9) +
  geom_hline(yintercept = 1, colour = "#8a8984", linewidth = 0.4, linetype = "22") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.7) +
  geom_text(aes(label = above), colour = "#52514e", size = 1.9, vjust = -0.9, show.legend = FALSE) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8+"),
                limits = c(0.2, cap * 1.9)) +
  scale_colour_manual(values = palette) + scale_shape_manual(values = shapes) +
  facet_wrap(~lambda, nrow = 1, labeller = label_parsed) +
  labs(x = "Distance from the nearest observed karyotype", y = "Error relative to the null predictor (log scale)") +
  theme_doc
ggsave(file.path(out_dir, "far_distance_ratio.pdf"), ratio_plot, width = 6.5, height = 2.9, device = grDevices::cairo_pdf)

rank_plot <- ggplot(plotted[is.finite(plotted$spearman), ], aes(distance, spearman, colour = label, shape = label)) +
  far_band() +
  geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.4, linetype = "22") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.7) +
  scale_colour_manual(values = palette) + scale_shape_manual(values = shapes) +
  facet_wrap(~lambda, nrow = 1, labeller = label_parsed) +
  labs(x = "Distance from the nearest observed karyotype", y = "Spearman correlation") + theme_doc
ggsave(file.path(out_dir, "far_distance_spearman.pdf"), rank_plot, width = 6.5, height = 2.9, device = grDevices::cairo_pdf)

## Tables: every method at distances 3 and 4 ---------------------------------
fmt <- function(v, digits) ifelse(is.finite(v), formatC(v, format = "f", digits = digits), "--")
table_tex <- function(wl, column, digits, title, methods = names(labels)) {
  rows <- vapply(methods, function(m) {
    vals <- unlist(lapply(waves, function(w) vapply(c(3L, 4L), function(d) {
      v <- wl[wl$method == m & wl$wavelength_class == w & wl$distance == d, column]
      if (length(v)) fmt(v, digits) else "--"
    }, "")))
    label <- sub("alfakR", "\\\\texttt{alfakR}", labels[[m]])
    paste0(label, " & ", paste(vals, collapse = " & "), " \\\\")
  }, "")
  # The minipage keeps each title on the page of its table.
  c("\\par\\medskip\\noindent\\begin{minipage}{\\linewidth}", paste0("\\textit{", title, "}"),
    "\\begin{center}", "\\footnotesize", "\\setlength{\\tabcolsep}{3.5pt}",
    "\\begin{tabular}{l*{8}{r}}", "\\toprule",
    paste0(" & \\multicolumn{2}{c}{$\\lambda=0.4$} & \\multicolumn{2}{c}{$\\lambda=0.8$}",
           " & \\multicolumn{2}{c}{$\\lambda=1.2$} & \\multicolumn{2}{c}{$\\lambda=1.6$} \\\\"),
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}",
    paste0("Method", paste(rep(" & $d=3$ & $d=4$", 4), collapse = ""), " \\\\"), "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}", "\\end{center}", "\\end{minipage}")
}
main_methods <- c("zero", "alfak2_two_layer_fixed", "alfakR_minobs20", "alfak2_sparse_joint", "alfak2_sparse_joint_balance")
writeLines(c(table_tex(held, "ratio", 2, "Error relative to the null predictor, held-out landscapes.", main_methods[-1]), "",
             table_tex(held, "spearman", 2, "Spearman correlation, held-out landscapes.", main_methods[-1])),
           file.path(out_dir, "far_main.tex"))
supp <- unlist(lapply(c("heldout", "dev"), function(split) {
  wl <- if (split == "heldout") held else dev
  name <- if (split == "heldout") "held-out" else "development"
  c(table_tex(wl, "rmse", 3, sprintf("Root mean squared error, %s landscapes.", name)), "",
    table_tex(wl, "bias", 3, sprintf("Mean error (prediction minus reference), %s landscapes.", name)), "",
    table_tex(wl, "spearman", 2, sprintf("Spearman correlation, %s landscapes.", name), names(labels)[-1]), "",
    table_tex(wl, "coverage95", 2, sprintf("Coverage of nominal 95\\%% intervals, %s landscapes.", name), names(labels)[-1]), "",
    table_tex(wl, "prediction_fraction", 2, sprintf("Fraction of karyotypes predicted, %s landscapes.", name), names(labels)[-1]), "")
}))
writeLines(supp, file.path(out_dir, "far_supp.tex"))

## Numbers quoted in the text -------------------------------------------------
options(width = 200)
show <- function(wl, methods, column, digits = 2) {
  x <- wl[wl$method %in% methods, c("label", "wavelength_class", "distance", column)]
  w <- stats::reshape(x, idvar = c("label", "distance"), timevar = "wavelength_class", direction = "wide")
  w <- w[order(w$label, w$distance), ]
  w[-(1:2)] <- round(w[-(1:2)], digits)
  print(w, row.names = FALSE)
}
for (split in c("heldout", "dev")) {
  wl <- if (split == "heldout") held else dev
  cat("==", split, "ratio\n"); show(wl, main_methods[-1], "ratio")
  cat("==", split, "spearman\n"); show(wl, main_methods[-1], "spearman")
  cat("==", split, "bias (null bias = minus mean truth)\n"); show(wl, main_methods, "bias", 3)
  cat("==", split, "coverage\n"); show(wl, c("alfak2_sparse_joint_balance", "alfakR_minobs20", "alfak2_two_layer_fixed"), "coverage95")
  cat("==", split, "near-field prediction fraction at d3/d4\n")
  print(stats::aggregate(prediction_fraction ~ label, wl[wl$distance >= 3 & grepl("two_layer", wl$method), ], max))
}
