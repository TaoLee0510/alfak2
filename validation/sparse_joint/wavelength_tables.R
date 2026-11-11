# Benchmark results by wavelength class, as LaTeX tables (doc/alfak2_method.tex).
#
#   Rscript validation/sparse_joint/wavelength_tables.R METRICS_CSV OUT_TEX MODE
#
# MODE "heldout": the methods of the held-out confirmation (Section 5).
# MODE "extrapolation": the second-layer methods of the two-layer pipeline
# against the sparse joint model (Supplement S2.3).
# Values are averaged over the scenarios of each landscape and then over the
# three landscapes of each wavelength; core scenarios and patient 1 only.

args <- commandArgs(TRUE)
m <- read.csv(args[[1]], stringsAsFactors = FALSE)
mode <- args[[3]]
tt <- function(x) paste0("\\texttt{", gsub("_", "\\\\_", x), "}")
methods <- if (mode == "heldout") {
  c("Null predictor" = "zero",
    "Two-layer pipeline" = "alfak2_two_layer_fixed",
    "alfakR" = "alfakR_minobs20",
    "Sparse joint, independent" = "alfak2_sparse_joint",
    "Sparse joint, balance" = "alfak2_sparse_joint_balance")
} else {
  # Short labels; the text of Supplement S2.3 maps them to the method names.
  c("Null predictor" = "zero",
    "Graph Gaussian" = "alfak2_two_layer_fixed",
    "Truncated GMRF" = "alfak2_two_layer_nearfield_gmrf",
    "Edge effect" = "alfak2_two_layer_edge_effect",
    "Path ensemble" = "alfak2_two_layer_path_ensemble",
    "Local polynomial" = "alfak2_two_layer_polynomial",
    "Nearest-neighbour GP" = "alfak2_two_layer_nngp",
    "Delta tree" = "alfak2_two_layer_delta_tree",
    "Trend filtering" = "alfak2_two_layer_trend_filtering",
    "Sparse joint, balance" = "alfak2_sparse_joint_balance")
}
display <- function(label) {
  if (label == "alfakR") return("\\texttt{alfakR} (\\texttt{minobs = 20})")
  if (grepl("_", label)) return(tt(label))
  label
}
x <- m[m$method %in% methods & m$role == "core" & m$scoreboard & m$stratum %in% c("d0", "d1", "d2"), ]
cols <- c("rmse", "spearman", "coverage95")
f <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
cell <- aggregate(x[cols], x[c("method", "wavelength_class", "landscape_id", "scenario", "stratum")], f)

land <- aggregate(cell[cols], cell[c("method", "wavelength_class", "landscape_id", "stratum")], f)
wl <- aggregate(land[cols], land[c("method", "wavelength_class", "stratum")], f)
waves <- c("w0.4", "w0.8", "w1.2", "w1.6")
value <- function(method, w, d, metric) {
  v <- wl[wl$method == method & wl$wavelength_class == w & wl$stratum == d, metric]
  if (!length(v) || is.na(v)) "--" else sprintf("%.3f", v)
}
table_tex <- function(metric, keep, header) {
  out <- c("\\begin{center}", "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
           "\\begin{tabular}{l*{12}{r}}", "\\toprule",
           paste0(header, " & \\multicolumn{3}{c}{$\\lambda=0.4$} & \\multicolumn{3}{c}{$\\lambda=0.8$}",
                  " & \\multicolumn{3}{c}{$\\lambda=1.2$} & \\multicolumn{3}{c}{$\\lambda=1.6$} \\\\"),
           "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}\\cmidrule(lr){11-13}",
           paste0(paste(rep(" & $d=0$ & $d=1$ & $d=2$", 4), collapse = ""), " \\\\"), "\\midrule")
  for (label in keep) {
    vals <- unlist(lapply(waves, function(w) vapply(c("d0", "d1", "d2"), function(d) value(methods[[label]], w, d, metric), "")))
    out <- c(out, paste0(display(label), " & ", paste(vals, collapse = " & "), " \\\\"))
  }
  c(out, "\\bottomrule", "\\end{tabular}", "\\end{center}")
}
has <- function(metric, label) any(is.finite(wl[wl$method == methods[[label]], metric]))
keep_rank <- Filter(function(l) has("spearman", l), names(methods))
keep_cov <- Filter(function(l) has("coverage95", l), names(methods))
tex <- c("\\paragraph{Root mean squared error.}", table_tex("rmse", names(methods), ""), "",
         "\\paragraph{Spearman correlation with the reference values.}", table_tex("spearman", keep_rank, ""), "",
         "\\paragraph{Coverage of nominal 95\\% intervals.}", table_tex("coverage95", keep_cov, ""))
writeLines(tex, args[[2]])
failed <- aggregate(status ~ method + wavelength_class, unique(x[c("dataset_id", "method", "wavelength_class", "status")]),
                    function(s) mean(s != "ok"))
print(failed[failed$status > 0, ], row.names = FALSE)
