# Benchmark results by sampling design (scenario) and wavelength, for the
# method document (Section 5.2 and Supplement S4).
#
#   Rscript validation/sparse_joint/scenario_results.R [OUT_DIR]
#
# Reads metrics_heldout.csv and metrics_dev.csv of benchmark/sparse_two_timepoint
# and writes, to OUT_DIR (default doc/generated):
#   scenario_rmse_ratio.pdf   error relative to the null predictor (held-out)
#   scenario_spearman.pdf     rank correlation with the reference (held-out)
#   scenario_main.tex         cells per karyotype and the default model's error
#                             ratio, by design and wavelength (held-out)
#   scenario_supp_<split>.tex error, rank correlation and coverage of every
#                             method by design and wavelength
#   scenario_ranges.tex       range over landscapes of the default model's ratio
# Each cell of a table holds the three landscapes of one wavelength under one
# design; ratios are the mean error over those landscapes divided by the mean
# error of the null predictor on the same karyotypes.

suppressMessages(library(ggplot2))
args <- commandArgs(TRUE)
out_dir <- if (length(args)) args[[1L]] else "doc/generated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results <- "benchmark/results/sparse_two_timepoint"

methods <- c(
  "Null predictor" = "zero",
  "Two-layer pipeline" = "alfak2_two_layer_fixed",
  "alfakR (minobs = 20)" = "alfakR_minobs20",
  "Sparse joint, independent" = "alfak2_sparse_joint",
  "Sparse joint, balance (default)" = "alfak2_sparse_joint_balance"
)
designs <- c(
  "very_sparse_50_50_180d" = "50/50, 180 d",
  "sparse_200_200_180d" = "200/200, 180 d",
  "sparse_200_200_360d" = "200/200, 360 d",
  "typical_500_700_360d" = "500/700, 360 d",
  "long_500_700_1080d" = "500/700, 1080 d",
  "bottleneck_500_700_360d" = "500/700, bottleneck",
  "unbalanced_100_3000_360d" = "100/3000, 360 d",
  "unbalanced_3000_100_360d" = "3000/100, 360 d"
)
waves <- c("w0.4", "w0.8", "w1.2", "w1.6")

load_split <- function(split) {
  m <- utils::read.csv(file.path(results, sprintf("metrics_%s.csv", split)), stringsAsFactors = FALSE)
  m <- m[m$role == "core" & m$scoreboard & m$stratum %in% c("d0", "d1", "d2") & m$method %in% methods, ]
  m$label <- names(methods)[match(m$method, methods)]
  m$design <- factor(designs[m$scenario], levels = designs)
  m$cells_per_karyotype <- (m$depth_t0 + m$depth_t1) / m$n_karyotypes
  m
}
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

# One row per method x wavelength x design x distance: means over landscapes,
# the error ratio to the null predictor, and its range over landscapes.
summarise_cells <- function(m) {
  key <- c("label", "method", "wavelength_class", "design", "stratum")
  cells <- stats::aggregate(m[c("rmse", "rmse_zero", "spearman", "coverage95")], m[key], mean_na)
  cells$ratio <- cells$rmse / cells$rmse_zero
  m$landscape_ratio <- m$rmse / m$rmse_zero
  low <- stats::aggregate(m[c("landscape_ratio", "spearman")], m[key], function(v) suppressWarnings(min(v, na.rm = TRUE)))
  high <- stats::aggregate(m[c("landscape_ratio", "spearman")], m[key], function(v) suppressWarnings(max(v, na.rm = TRUE)))
  names(low)[6:7] <- c("ratio_low", "spearman_low")
  names(high)[6:7] <- c("ratio_high", "spearman_high")
  merge(merge(cells, low, by = key), high, by = key)
}

held <- load_split("heldout")
dev <- load_split("dev")
held_cells <- summarise_cells(held)
dev_cells <- summarise_cells(dev)

## Figures -------------------------------------------------------------------

# Categorical slots 1-4 of the reference palette, in fixed order; shapes are the
# secondary encoding (the third and fourth hues are below 3:1 on white).
palette <- c("Sparse joint, balance (default)" = "#2a78d6", "Sparse joint, independent" = "#eb6834",
             "Two-layer pipeline" = "#1baf7a", "alfakR (minobs = 20)" = "#eda100")
shapes <- c("Sparse joint, balance (default)" = 16, "Sparse joint, independent" = 17,
            "Two-layer pipeline" = 15, "alfakR (minobs = 20)" = 18)
# The error axis stops at `cap`: alfakR reaches about 300 times the null
# predictor's error, which would compress every other series. Values above it
# are drawn at the cap and labelled with their ratio.
cap <- 8
plotted <- held_cells[held_cells$label %in% names(palette), ]
plotted$label <- factor(plotted$label, levels = names(palette))
plotted$lambda <- factor(sub("^w", "lambda == ", plotted$wavelength_class),
                         levels = sub("^w", "lambda == ", waves))
plotted$distance <- factor(paste0("d == ", sub("^d", "", plotted$stratum)), levels = c("d == 0", "d == 1", "d == 2"))
theme_doc <- theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "#52514e"),
        axis.text.y = element_text(colour = "#52514e"), strip.text = element_text(colour = "#1f1f1e"),
        legend.position = "top", legend.title = element_blank(), panel.spacing = unit(6, "pt"))
separators <- geom_vline(xintercept = c(4.5, 5.5, 6.5), colour = "#d7d6d0", linewidth = 0.3)
dodge <- position_dodge(width = 0.6)

plotted$ratio_shown <- pmin(plotted$ratio, cap)
plotted$low_shown <- pmin(plotted$ratio_low, cap)
plotted$high_shown <- pmin(plotted$ratio_high, cap)
plotted$above <- ifelse(plotted$ratio > cap, sprintf("%.0f", plotted$ratio), "")
ratio_plot <- ggplot(plotted, aes(design, ratio_shown, colour = label, shape = label, group = label)) +
  separators +
  geom_hline(yintercept = 1, colour = "#8a8984", linewidth = 0.4, linetype = "22") +
  geom_hline(yintercept = cap, colour = "#d7d6d0", linewidth = 0.3, linetype = "13") +
  geom_linerange(aes(ymin = low_shown, ymax = high_shown), position = dodge, linewidth = 0.45) +
  geom_point(position = dodge, size = 1.5) +
  geom_text(aes(label = above), position = dodge, colour = "#52514e", size = 1.9, vjust = -0.9,
            show.legend = FALSE) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8+"),
                limits = c(0.25, cap * 1.9)) +
  scale_colour_manual(values = palette) + scale_shape_manual(values = shapes) +
  facet_grid(lambda ~ distance, labeller = label_parsed) +
  labs(x = NULL, y = "Error relative to the null predictor (log scale)") + theme_doc
ggsave(file.path(out_dir, "scenario_rmse_ratio.pdf"), ratio_plot, width = 6.5, height = 8.2, device = grDevices::cairo_pdf)

rank_plot <- ggplot(plotted[is.finite(plotted$spearman), ], aes(design, spearman, colour = label, shape = label)) +
  separators +
  geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.4, linetype = "22") +
  geom_linerange(aes(ymin = spearman_low, ymax = spearman_high), position = dodge, linewidth = 0.45) +
  geom_point(position = dodge, size = 1.5) +
  scale_colour_manual(values = palette) + scale_shape_manual(values = shapes) +
  facet_grid(lambda ~ distance, labeller = label_parsed) +
  labs(x = NULL, y = "Spearman correlation with the reference values") + theme_doc
ggsave(file.path(out_dir, "scenario_spearman.pdf"), rank_plot, width = 6.5, height = 8.2, device = grDevices::cairo_pdf)

## Tables --------------------------------------------------------------------

fmt <- function(v) ifelse(is.finite(v), sprintf("%.2f", v), "--")
fmt3 <- function(v) ifelse(is.finite(v), sprintf("%.3f", v), "--")
header <- c("\\begin{center}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}")
lambda_head <- function(first) c(
  paste0(first, " & \\multicolumn{3}{c}{$\\lambda=0.4$} & \\multicolumn{3}{c}{$\\lambda=0.8$}",
         " & \\multicolumn{3}{c}{$\\lambda=1.2$} & \\multicolumn{3}{c}{$\\lambda=1.6$} \\\\"),
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-13}",
  paste0(paste(rep(" & $d=0$ & $d=1$ & $d=2$", 4), collapse = ""), " \\\\"), "\\midrule")
cell_value <- function(cells, label, design, column, format) {
  vapply(waves, function(w) {
    vapply(c("d0", "d1", "d2"), function(d) {
      v <- cells[cells$label == label & cells$design == design & cells$wavelength_class == w & cells$stratum == d, column]
      if (length(v)) format(v) else "--"
    }, "")
  }, character(3))
}
design_rows <- function(cells, label, column, format) {
  vapply(designs, function(design) {
    paste0(design, " & ", paste(cell_value(cells, label, design, column, format), collapse = " & "), " \\\\")
  }, "")
}

# Main text: cells per observed karyotype, then the default model's ratio.
cpk <- stats::aggregate(cells_per_karyotype ~ design + wavelength_class,
                        held[held$method == "zero" & held$stratum == "d0", ], stats::median)
cpk_rows <- vapply(designs, function(design) {
  v <- vapply(waves, function(w) cpk$cells_per_karyotype[cpk$design == design & cpk$wavelength_class == w], 0)
  paste0(design, " & ", paste(sprintf("%.1f", v), collapse = " & "), " \\\\")
}, "")
main <- c("\\noindent\\textit{Cells per observed karyotype (median over landscapes).}", header,
          "\\begin{tabular}{lrrrr}", "\\toprule",
          "Design & $\\lambda=0.4$ & $\\lambda=0.8$ & $\\lambda=1.2$ & $\\lambda=1.6$ \\\\", "\\midrule",
          cpk_rows, "\\bottomrule", "\\end{tabular}", "\\end{center}", "",
          "\\noindent\\textit{Error of the sparse joint model (balance) relative to the null predictor.}",
          header, "\\begin{tabular}{l*{12}{r}}", "\\toprule", lambda_head("Design"),
          design_rows(held_cells, "Sparse joint, balance (default)", "ratio", fmt),
          "\\bottomrule", "\\end{tabular}", "\\end{center}")
writeLines(main, file.path(out_dir, "scenario_main.tex"))

# Range over the three landscapes of the default model's ratio, two
# wavelengths per table for width.
range_value <- function(design, w, d) {
  r <- held_cells[held_cells$label == "Sparse joint, balance (default)" & held_cells$design == design &
                    held_cells$wavelength_class == w & held_cells$stratum == d, ]
  sprintf("%.2f--%.2f", r$ratio_low, r$ratio_high)
}
range_table <- function(pair) {
  c(header, "\\begin{tabular}{l*{6}{r}}", "\\toprule",
    paste0("Design", paste(sprintf(" & \\multicolumn{3}{c}{$\\lambda=%s$}", sub("^w", "", pair)), collapse = ""), " \\\\"),
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
    paste0(paste(rep(" & $d=0$ & $d=1$ & $d=2$", 2), collapse = ""), " \\\\"), "\\midrule",
    vapply(designs, function(design) {
      paste0(design, " & ", paste(unlist(lapply(pair, function(w) vapply(c("d0", "d1", "d2"),
                                                                          function(d) range_value(design, w, d), ""))),
                                  collapse = " & "), " \\\\")
    }, ""), "\\bottomrule", "\\end{tabular}", "\\end{center}")
}
writeLines(c(range_table(c("w0.4", "w0.8")), "", range_table(c("w1.2", "w1.6"))), file.path(out_dir, "scenario_ranges.tex"))

# Supplement: every method, by design, with short labels and a smaller font
# for width.
short <- c("Null predictor" = "Null predictor", "Two-layer pipeline" = "Two-layer",
           "alfakR (minobs = 20)" = "\\texttt{alfakR}", "Sparse joint, independent" = "Sparse joint, independent",
           "Sparse joint, balance (default)" = "Sparse joint, balance")
supp_table <- function(cells, column, labels, title) {
  rows <- unlist(lapply(seq_along(designs), function(k) {
    design <- designs[[k]]
    c(if (k > 1L) "\\addlinespace[2pt]" else NULL,
      paste0("\\multicolumn{13}{l}{\\textit{", design, "}} \\\\"),
      vapply(labels, function(label) {
        paste0("\\quad ", short[[label]], " & ",
               paste(cell_value(cells, label, design, column, fmt3), collapse = " & "), " \\\\")
      }, ""))
  }))
  c(paste0("\\noindent\\textit{", title, "}"), "\\begin{center}", "\\scriptsize", "\\setlength{\\tabcolsep}{2.5pt}",
    "\\begin{tabular}{l*{12}{r}}", "\\toprule",
    lambda_head("Method"), rows, "\\bottomrule", "\\end{tabular}", "\\end{center}")
}
for (split in c("heldout", "dev")) {
  cells <- if (split == "heldout") held_cells else dev_cells
  set_name <- if (split == "heldout") "held-out" else "development"
  writeLines(c(
    supp_table(cells, "rmse", names(methods), sprintf("Root mean squared error, %s landscapes.", set_name)), "",
    supp_table(cells, "spearman", names(methods)[-1L], sprintf("Spearman correlation, %s landscapes.", set_name)), "",
    supp_table(cells, "coverage95", names(methods)[-1L], sprintf("Coverage of nominal 95\\%% intervals, %s landscapes.", set_name))
  ), file.path(out_dir, sprintf("scenario_supp_%s.tex", split)))
}

## Numbers quoted in the text --------------------------------------------------
options(width = 200)
show <- function(cells, label, column) {
  x <- cells[cells$label == label, c("wavelength_class", "design", "stratum", column)]
  w <- stats::reshape(x, idvar = c("design", "stratum"), timevar = "wavelength_class", direction = "wide")
  w <- w[order(w$stratum, w$design), ]
  w[-(1:2)] <- round(w[-(1:2)], 2)
  print(w, row.names = FALSE)
}
cat("== held-out ratio, sparse joint balance\n"); show(held_cells, "Sparse joint, balance (default)", "ratio")
cat("== held-out ratio, two-layer\n"); show(held_cells, "Two-layer pipeline", "ratio")
cat("== held-out Spearman, sparse joint balance\n"); show(held_cells, "Sparse joint, balance (default)", "spearman")
cat("== held-out Spearman, two-layer\n"); show(held_cells, "Two-layer pipeline", "spearman")
cat("== held-out coverage, sparse joint balance\n"); show(held_cells, "Sparse joint, balance (default)", "coverage95")
cat("== development ratio, sparse joint balance\n"); show(dev_cells, "Sparse joint, balance (default)", "ratio")
cat("== held-out ratio, alfakR: cells below 1\n")
a <- held_cells[held_cells$label == "alfakR (minobs = 20)", ]
cat(sum(a$ratio < 1), "of", nrow(a), "; range", round(range(a$ratio), 2), "\n")
cat("== ratio range over all plotted series\n"); print(round(range(c(plotted$ratio_low, plotted$ratio_high)), 2))
cat("== held-out ratio range over landscapes, balance: max of ratio_high\n")
print(max(held_cells$ratio_high[held_cells$label == "Sparse joint, balance (default)"]))
