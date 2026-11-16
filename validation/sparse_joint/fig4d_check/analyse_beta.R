#!/usr/bin/env Rscript
# Fig. 4d check, step 3b summary: which missegregation rate the Salehi data
# prefer (refit_beta.R rows, copied to results/beta/), and Fig. 4d's context
# model at each rate and at each trajectory's preferred rate.
#
#   Rscript validation/sparse_joint/fig4d_check/analyse_beta.R
#
# Writes results/beta.md.

suppressPackageStartupMessages(library(glmmTMB))
base <- "validation/sparse_joint/fig4d_check/results"
files <- list.files(file.path(base, c("beta", "beta_ext")), pattern = "\\.rds$", full.names = TRUE)
x <- do.call(rbind, lapply(files, readRDS))
md <- function(d) {
  if (is.null(d) || !nrow(d)) return("(none)")
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

obj <- unique(x[is.finite(x$objective), c("fi", "context", "treated", "beta_multiple", "objective")])
best <- do.call(rbind, by(obj, obj$fi, function(d) {
  b <- d[which.min(d$objective), ]
  base <- d$objective[d$beta_multiple == 1]
  data.frame(fi = b$fi, context = b$context, best_multiple = b$beta_multiple,
             gain_vs_default = base - b$objective)
}))
pref <- as.data.frame.matrix(table(best$context, best$best_multiple))
pref <- cbind(context = rownames(pref), pref)
gain <- aggregate(gain_vs_default ~ context, best, function(v) c(median = stats::median(v), sum = sum(v)))
gain <- do.call(data.frame, gain)
# Summed objective by context and rate, relative to the default rate.
tot <- aggregate(objective ~ context + beta_multiple, obj, sum)
tot$relative <- tot$objective - ave(ifelse(tot$beta_multiple == 1, tot$objective, NA), tot$context,
                                    FUN = function(v) max(v, na.rm = TRUE))

d <- x[is.finite(x$deltaf), ]
d$pdx <- sub("^(SA[0-9]+).*", "\\1", d$fi)
d$context <- stats::relevel(factor(d$context), ref = "invitro")
d$absdf <- abs(d$deltaf) + 1e-6
model <- function(dd, label) {
  m <- tryCatch(glmmTMB(absdf ~ context + (1 | pdx / fi / k_id), family = Gamma("log"), data = dd), error = function(e) NULL)
  if (is.null(m)) return(data.frame(rate = label, PDX_v = NA, p = NA))
  cf <- summary(m)$coefficients$cond["contextPDX", ]
  data.frame(rate = label, invitro_absdf = stats::median(dd$absdf[dd$context == "invitro"]),
             PDX_absdf = stats::median(dd$absdf[dd$context == "PDX"]), PDX_v = cf[["Estimate"]], p = cf[["Pr(>|z|)"]])
}
models <- do.call(rbind, c(lapply(sort(unique(d$beta_multiple)), function(m) model(d[d$beta_multiple == m, ], paste0("x", m))),
  list(model(merge(d, best[c("fi", "best_multiple")], by = "fi")[merge(d, best[c("fi", "best_multiple")], by = "fi")$beta_multiple ==
                 merge(d, best[c("fi", "best_multiple")], by = "fi")$best_multiple, ], "preferred per trajectory"))))

writeLines(c("# Fig. 4d check, step 3b: missegregation rate (beta = 5e-5 x multiple)", "",
             "## Preferred rate per trajectory (lowest objective), count", "", md(pref), "",
             "## Gain in log marginal likelihood of the preferred rate over the default", "", md(gain), "",
             "## Summed objective by context, relative to the default rate (negative: better)", "", md(tot), "",
             "## Fig. 4d context model at each rate (PDX_v > 0 is the paper's direction)", "", md(models)),
           file.path(base, "beta.md"))
cat(readLines(file.path(base, "beta.md")), sep = "\n")
