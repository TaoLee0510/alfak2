# The three accuracy measures by wavelength and distance (Section 5.2 of
# doc/alfak2_method.tex).
#
#   Rscript validation/sparse_joint/accuracy_by_wavelength.R [OUT_DIR]
#
# Writes to OUT_DIR (default doc/generated):
#   accuracy_by_wavelength.pdf  numerical accuracy (error relative to the null
#                               predictor), direction accuracy and Spearman
#                               correlation of the main methods, held-out
#   wavelength_direction.tex    direction accuracy table, held-out
# Values are averaged over the scenarios of each landscape and then over the
# three landscapes of each wavelength (24 datasets per point); the error ratio
# divides the mean error by the mean error of the null predictor. Direction
# accuracy comes from direction_accuracy.R.

suppressMessages({
  library(ggplot2)
  library(patchwork)
})
args <- commandArgs(TRUE)
out_dir <- if (length(args)) args[[1L]] else "doc/generated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
results <- "benchmark/results/sparse_two_timepoint"

methods <- c("Sparse joint, balance (default)" = "alfak2_sparse_joint_balance",
             "Sparse joint, independent" = "alfak2_sparse_joint",
             "Two-layer pipeline" = "alfak2_two_layer_fixed",
             "alfakR (minobs = 20)" = "alfakR_minobs20")
waves <- c("w0.4", "w0.8", "w1.2", "w1.6")
mean_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

load_split <- function(split) {
  m <- utils::read.csv(file.path(results, sprintf("metrics_%s.csv", split)), stringsAsFactors = FALSE)
  m <- m[m$role == "core" & m$scoreboard & m$stratum %in% c("d0", "d1", "d2") & m$method %in% methods, ]
  direction <- utils::read.csv(file.path(results, "sweeps", "direction", sprintf("direction_%s.csv", split)),
                               stringsAsFactors = FALSE)
  m <- merge(m, direction[c("dataset_id", "method", "stratum", "direction")], by = c("dataset_id", "method", "stratum"),
             all.x = TRUE)
  cols <- c("rmse", "rmse_zero", "direction", "spearman")
  cell <- stats::aggregate(m[cols], m[c("method", "wavelength_class", "landscape_id", "scenario", "stratum")], mean_na)
  land <- stats::aggregate(cell[cols], cell[c("method", "wavelength_class", "landscape_id", "stratum")], mean_na)
  wl <- stats::aggregate(land[cols], land[c("method", "wavelength_class", "stratum")], mean_na)
  wl$ratio <- wl$rmse / wl$rmse_zero
  wl$distance <- as.integer(sub("^d", "", wl$stratum))
  wl$label <- factor(names(methods)[match(wl$method, methods)], levels = names(methods))
  wl$lambda <- factor(sub("^w", "lambda == ", wl$wavelength_class), levels = sub("^w", "lambda == ", waves))
  wl
}
held <- load_split("heldout")
dev <- load_split("dev")

## Figure ----------------------------------------------------------------------
palette <- c("Sparse joint, balance (default)" = "#2a78d6", "Sparse joint, independent" = "#eb6834",
             "Two-layer pipeline" = "#1baf7a", "alfakR (minobs = 20)" = "#eda100")
shapes <- c("Sparse joint, balance (default)" = 16, "Sparse joint, independent" = 17,
            "Two-layer pipeline" = 15, "alfakR (minobs = 20)" = 18)
theme_doc <- theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"), strip.text = element_text(colour = "#1f1f1e"),
        legend.position = "top", legend.title = element_blank(), panel.spacing = unit(8, "pt"))
cap <- 8
held$ratio_shown <- pmin(held$ratio, cap)
held$above <- ifelse(held$ratio > cap, sprintf("%.0f", held$ratio), "")
row_plot <- function(y, reference, ylab, scale, labels = NULL, strips = FALSE, xlab = NULL) {
  p <- ggplot(held[is.finite(held[[y]]), ], aes(distance, .data[[y]], colour = label, shape = label)) +
    geom_hline(yintercept = reference, colour = "#8a8984", linewidth = 0.4, linetype = "22") +
    geom_line(linewidth = 0.5) + geom_point(size = 1.7) + scale +
    scale_x_continuous(breaks = 0:2, limits = c(-0.2, 2.2)) +
    scale_colour_manual(values = palette) + scale_shape_manual(values = shapes) +
    facet_wrap(~lambda, nrow = 1, labeller = label_parsed) +
    labs(x = xlab, y = ylab) + theme_doc
  if (!is.null(labels)) p <- p + geom_text(aes(label = .data[[labels]]), colour = "#52514e", size = 1.9, vjust = -0.9,
                                           show.legend = FALSE)
  if (!strips) p <- p + theme(strip.text = element_blank())
  p
}
p1 <- row_plot("ratio_shown", 1, "Error relative to\nthe null (log scale)",
               scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8), labels = c("0.25", "0.5", "1", "2", "4", "8+"),
                             limits = c(0.25, cap * 1.9)), labels = "above", strips = TRUE)
p2 <- row_plot("direction", 0.5, "Direction accuracy", scale_y_continuous(limits = c(0.3, 1), breaks = seq(0.3, 1, 0.1)))
p3 <- row_plot("spearman", 0, "Spearman correlation", scale_y_continuous(limits = c(-0.1, 0.6), breaks = seq(-0.1, 0.6, 0.1)),
               xlab = "Distance from the nearest observed karyotype")
figure <- (p1 / p2 / p3) + plot_layout(guides = "collect") & theme(legend.position = "top")
ggsave(file.path(out_dir, "accuracy_by_wavelength.pdf"), figure, width = 6.5, height = 6.6, device = grDevices::cairo_pdf)

## Table -----------------------------------------------------------------------
display <- c("Two-layer pipeline" = "Two-layer pipeline",
             "alfakR (minobs = 20)" = "\\texttt{alfakR} (\\texttt{minobs = 20})",
             "Sparse joint, independent" = "Sparse joint, independent",
             "Sparse joint, balance (default)" = "Sparse joint, balance")
rows <- vapply(names(display), function(label) {
  vals <- unlist(lapply(waves, function(w) vapply(c("d0", "d1", "d2"), function(d) {
    v <- held$direction[held$label == label & held$wavelength_class == w & held$stratum == d]
    if (length(v) && is.finite(v)) sprintf("%.3f", v) else "--"
  }, "")))
  paste0(display[[label]], " & ", paste(vals, collapse = " & "), " \\\\")
}, "")
writeLines(c("\\begin{center}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}", "\\begin{tabular}{l*{12}{r}}", "\\toprule",
             paste0(" & \\multicolumn{3}{c}{$\\lambda=0.4$} & \\multicolumn{3}{c}{$\\lambda=0.8$}",
                    " & \\multicolumn{3}{c}{$\\lambda=1.2$} & \\multicolumn{3}{c}{$\\lambda=1.6$} \\\\"),
             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-13}",
             paste0(paste(rep(" & $d=0$ & $d=1$ & $d=2$", 4), collapse = ""), " \\\\"), "\\midrule",
             rows, "\\bottomrule", "\\end{tabular}", "\\end{center}"),
           file.path(out_dir, "wavelength_direction.tex"))

## Numbers quoted in the text ---------------------------------------------------
for (split in c("heldout", "dev")) {
  wl <- if (split == "heldout") held else dev
  for (label in names(methods)) {
    x <- wl[wl$label == label, ]
    cat(sprintf("%-8s %-32s ratio %.2f-%.2f  direction %.2f-%.2f  spearman %.2f-%.2f\n", split, label,
                min(x$ratio), max(x$ratio), min(x$direction), max(x$direction),
                min(x$spearman, na.rm = TRUE), max(x$spearman, na.rm = TRUE)))
  }
}
