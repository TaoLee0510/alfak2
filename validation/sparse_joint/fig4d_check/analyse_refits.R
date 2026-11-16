#!/usr/bin/env Rscript
# Fig. 4d check, step 3 summary: delta f of the refits (refit_variants.R,
# copied to results/refits/) by variant and context, and the paper's Fig. 4d
# context model on each variant:
#   glmmTMB(|delta f| + 1e-6 ~ context + (1 | pdx / fi / k_id), Gamma(log))
# on all selected trajectories (no lineage bootstrap; the benchmark's 200
# draws of non-overlapping lineages are for the published comparison, this is
# for the direction under each variant). PDX_v > 0: |delta f| larger in PDX,
# as in the paper.
#
#   Rscript validation/sparse_joint/fig4d_check/analyse_refits.R
#
# Writes results/refits.md.

suppressPackageStartupMessages(library(glmmTMB))
base <- "validation/sparse_joint/fig4d_check/results"
files <- list.files(file.path(base, "refits"), pattern = "\\.rds$", full.names = TRUE)
x <- do.call(rbind, lapply(files, readRDS))
x <- x[is.finite(x$deltaf), ]
x$pdx <- sub("^(SA[0-9]+).*", "\\1", x$fi)
x$context <- stats::relevel(factor(x$context), ref = "invitro")

md <- function(d) {
  if (is.null(d) || !nrow(d)) return("(none)")
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

med <- aggregate(cbind(deltaf, absdf = abs(deltaf)) ~ variant + context, x, stats::median)
wide <- reshape(med, idvar = "variant", timevar = "context", direction = "wide")
wide$ratio_PDX_over_invitro <- wide$absdf.PDX / wide$absdf.invitro
by_src <- aggregate(deltaf ~ variant + context + source, x, stats::median)
by_src <- reshape(by_src, idvar = c("variant", "source"), timevar = "context", direction = "wide")
by_src <- by_src[order(by_src$source, by_src$variant), ]

models <- do.call(rbind, lapply(unique(x$variant), function(v) {
  d <- x[x$variant == v, ]
  d$absdf <- abs(d$deltaf) + 1e-6
  m <- tryCatch(glmmTMB(absdf ~ context + (1 | pdx / fi / k_id), family = Gamma("log"), data = d),
                error = function(e) NULL)
  if (is.null(m)) return(data.frame(variant = v, PDX_v = NA, se = NA, p = NA))
  cf <- summary(m)$coefficients$cond
  data.frame(variant = v, PDX_v = cf["contextPDX", "Estimate"], se = cf["contextPDX", "Std. Error"],
             p = cf["contextPDX", "Pr(>|z|)"])
}))
fits <- aggregate(fi ~ variant + context, unique(x[c("variant", "context", "fi")]), length)
errors <- unique(x$error[!is.na(x$error)])

writeLines(c("# Fig. 4d check, step 3: refits under variants (alfak2, two time points)", "",
             "## Median delta f and |delta f|", "", md(wide), "",
             "## Median delta f by the source of the neighbour's value", "", md(by_src), "",
             "## Fig. 4d context model (PDX_v > 0 is the paper's direction)", "", md(models), "",
             "## Trajectories per variant", "", md(fits), "",
             if (length(errors)) c("## Errors", "", errors) else "No fit failed."),
           file.path(base, "refits.md"))
cat(readLines(file.path(base, "refits.md")), sep = "\n")
