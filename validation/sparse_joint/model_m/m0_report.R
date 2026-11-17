# Stage M, step 0: tables of what m0.R found. Sourced by m0.R --stage=report;
# writes results/m0/report.md and the tables behind it.

mm_report <- function() {
  lines <- character()
  add <- function(...) lines <<- c(lines, ...)
  fmt <- function(x, digits = 3L) {
    x[] <- lapply(x, function(v) {
      if (is.numeric(v) && any(v != round(v), na.rm = TRUE)) formatC(v, digits = digits, format = "g") else v
    })
    x
  }
  table_md <- function(x, digits = 3L) {
    x <- fmt(as.data.frame(x), digits)
    c(paste0("| ", paste(names(x), collapse = " | "), " |"),
      paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|"),
      apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  read_dir <- function(dir) {
    files <- list.files(dir, pattern = "[.]rds$", full.names = TRUE, recursive = TRUE)
    lapply(files, readRDS)
  }
  add("# Stage M, step 0", "", sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M")), "")

  # Decomposition.
  dec <- read_dir(file.path(mm_out, "decompose"))
  if (length(dec)) {
    rows <- do.call(rbind, lapply(dec, function(r) {
      top <- r$states[1L, ]
      data.frame(set = r$set, id = r$id, history = r$history, n_states = r$n_states, n_random = r$n_random,
                 refit_vs_c1 = r$objective - r$c1_objective, joint_r_vs_tmb = r$joint_r - r$joint_tmb,
                 eig1 = r$eigen[[1L]], eig2 = r$eigen[[2L]],
                 zeta_minor_share = if ("zeta_minor" %in% names(r$loading)) r$loading[["zeta_minor"]] else 0,
                 zeta_major_share = if ("zeta_major" %in% names(r$loading)) r$loading[["zeta_major"]] else 0,
                 g_share = sum(r$loading[grepl("^g_", names(r$loading))]),
                 along_prior = r$along[["prior_gn"]], along_balance2 = r$along[["balance_second"]],
                 along_expected = r$along[["likelihood_expected"]], along_residual = r$along[["likelihood_residual"]],
                 top_y0 = top$y0, top_y1 = top$y1, top_share_own = top$share_own, top_share_major = top$share_major,
                 top_share_minor = top$share_minor, top_diag_residual = top$diag_residual,
                 top_residual_formula = top$residual_formula,
                 stringsAsFactors = FALSE)
    }))
    utils::write.csv(rows, file.path(mm_out, "decompose.csv"), row.names = FALSE)
    add("## Checks", "",
        sprintf("- R translation of the joint density against TMB at the mode: largest absolute difference %.2g (%d fits).",
                max(abs(rows$joint_r_vs_tmb)), nrow(rows)),
        sprintf("- Refit objective against C1's record: largest absolute difference %.2g.",
                max(abs(rows$refit_vs_c1), na.rm = TRUE)), "")
    by_set <- do.call(rbind, lapply(split(rows, rows$set), function(x) {
      data.frame(set = x$set[[1L]], fits = nrow(x), median_eig1 = stats::median(x$eig1),
                 eig1_below_2e3 = sum(x$eig1 < 2e-3), zeta_minor_share_median = stats::median(x$zeta_minor_share),
                 zeta_minor_share_ge_0.9 = sum(x$zeta_minor_share >= 0.9),
                 residual_most_negative = sum(x$along_residual < 0 & x$along_residual <= x$along_balance2),
                 balance2_most_negative = sum(x$along_balance2 < 0 & x$along_balance2 < x$along_residual),
                 median_prior = stats::median(x$along_prior), median_expected = stats::median(x$along_expected),
                 median_balance2 = stats::median(x$along_balance2), median_residual = stats::median(x$along_residual))
    }))
    add("## The flattest direction of the inner Hessian", "",
        "Curvature along the eigenvector of the smallest eigenvalue, split by term (`prior`: priors with the balance",
        "prior at zero residual; `balance2`: the balance prior's second-order part; `expected`: likelihood with counts",
        "at their expectations; `residual`: the rest of the likelihood, sum (y - lambda) d2 log lambda).", "",
        table_md(by_set), "")
    deg <- rows[rows$set == "degenerate", ]
    if (nrow(deg)) {
      deg$y0_band <- cut(deg$top_y0, c(-1, 0, 2, Inf), c("0", "1-2", "3+"))
      deg$y1_band <- cut(deg$top_y1, c(-1, 0, 2, 5, Inf), c("0", "1-2", "3-5", "6+"))
      add("Leading non-major state of each degenerate fit (counts at the first and second time point):", "",
          table_md(as.data.frame.matrix(table(y0 = deg$y0_band, y1 = deg$y1_band))), "",
          sprintf("Its expectation at the second time point: own growth median %.2f, inflow from majors %.2f, from non-majors %.3f.",
                  stats::median(deg$top_share_own), stats::median(deg$top_share_major), stats::median(deg$top_share_minor)),
          sprintf("Residual curvature on its own diagonal against -(y1 - lambda) pi (1 - pi): largest difference %.2g.",
                  max(abs(deg$top_diag_residual - deg$top_residual_formula))), "")
    }
  }

  # Cliffs.
  cl <- read_dir(file.path(mm_out, "cliffs"))
  if (length(cl)) {
    steps <- do.call(rbind, lapply(cl, function(r) r$steps))
    utils::write.csv(steps, file.path(mm_out, "cliffs.csv"), row.names = FALSE)
    per <- do.call(rbind, lapply(split(steps, paste(steps$id, steps$history)), function(x) {
      jump <- which.max(abs(x$f_change))
      data.frame(id = x$id[[1L]], history = x$history[[1L]], max_f_change = max(abs(x$f_change), na.rm = TRUE),
                 not_finite = sum(!is.finite(x$f_change)),
                 joint_change = abs(x$joint_change[jump]), marginal_change = abs(x$marginal_change[jump]),
                 zeta_minor = x$max_zeta_minor[jump], zeta_major = x$max_zeta_major[jump], g = x$max_g[jump],
                 top_major = x$top_major[jump], top_y0 = x$top_y0[jump], top_y1 = x$top_y1[jump],
                 stringsAsFactors = FALSE)
    }))
    per$jump <- per$max_f_change > 0.05 | per$not_finite > 0
    add("## Rejected candidates of cv1: what moves across a step of 1e-3 in a scale", "",
        sprintf("%d candidates; %d jump (objective changes by more than 0.05, or is not finite next to the point).",
                nrow(per), sum(per$jump)), "",
        "At the largest change of each jumping candidate:", "",
        table_md(data.frame(
          statistic = c("objective change", "joint density change", "stage-M density change (grid rule)",
                        "largest move, non-major zeta", "largest move, major zeta", "largest move, g"),
          median = c(stats::median(per$max_f_change[per$jump]), stats::median(per$joint_change[per$jump], na.rm = TRUE),
                     stats::median(per$marginal_change[per$jump], na.rm = TRUE),
                     stats::median(per$zeta_minor[per$jump]), stats::median(per$zeta_major[per$jump]),
                     stats::median(per$g[per$jump])),
          max = c(max(per$max_f_change[per$jump]), max(per$joint_change[per$jump], na.rm = TRUE),
                  max(per$marginal_change[per$jump], na.rm = TRUE), max(per$zeta_minor[per$jump]),
                  max(per$zeta_major[per$jump]), max(per$g[per$jump])))), "",
        "The state whose initial abundance moves most:", "",
        table_md(as.data.frame.matrix(table(major = per$top_major[per$jump],
                                            counts = paste0(pmin(per$top_y0[per$jump], 3), "->",
                                                            pmin(per$top_y1[per$jump], 3))))), "",
        "(counts capped at 3.)", "")
  }

  # Quadrature.
  qd <- lapply(list.files(file.path(mm_out, "quadrature"), pattern = "[.]rds$", full.names = TRUE,
                          recursive = TRUE), readRDS)
  qd <- Filter(function(x) is.data.frame(x) && "label" %in% names(x), qd)
  if (length(qd)) {
    q <- do.call(rbind, qd)
    q <- q[q$weight != 0, ]
    rules <- intersect(names(mm_rules), names(q))
    q$y0_band <- cut(q$y0, c(-1, 0, 4, 19, 99, Inf), c("0", "1-4", "5-19", "20-99", "100+"))
    per_rule <- do.call(rbind, lapply(rules, function(r) {
      e <- abs(q[[r]])
      data.frame(rule = r, nodes = mm_rule_nodes(mm_rules[[r]]), states = length(e), max = max(e),
                 p999 = stats::quantile(e, 0.999), above_1e6 = sum(e > 1e-6),
                 worst_y0 = q$y0[which.max(e)], worst_y1 = q$y1[which.max(e)],
                 worst_u = (q$data_centre[which.max(e)] - q$m[which.max(e)]) / q$s[which.max(e)])
    }))
    add("## Quadrature on the states of the fits", "",
        sprintf("%d non-major states (%d fits, at the mode and at two moved points); error of the log integral",
                nrow(q), length(qd)),
        "against the exact rule (trapezoid, step 0.002 on +-40 prior sd). `worst_u`: where the counts put z, in prior",
        "sd from the prior mean, for the worst state.", "", table_md(per_rule), "",
        "Largest error by count at the first time point:", "",
        table_md(cbind(y0 = levels(q$y0_band), as.data.frame(lapply(stats::setNames(rules, rules), function(r) {
          as.numeric(tapply(abs(q[[r]]), q$y0_band, max))
        })))), "",
        sprintf("Data centre beyond 6 prior sd from the prior mean: %d states; beyond 8: %d.",
                sum(abs(q$data_centre - q$m) / q$s > 6), sum(abs(q$data_centre - q$m) / q$s > 8)), "")
  }
  stress_path <- file.path(mm_out, "quadrature", "stress.rds")
  if (file.exists(stress_path)) {
    st <- readRDS(stress_path)
    rules <- intersect(names(mm_rules), names(st))
    st$regime <- ifelse(abs(st$m - log(st$y0 + 0.5)) / st$s > 6, "counts far from prior", "counts near prior")
    add("## Quadrature under stress (synthetic states)", "",
        sprintf("%d states: counts 0-149 at the first time point, 0-300 at the second, prior mean -10 to 3, prior sd 1 and 3.", nrow(st)),
        sprintf("Exact rule against stats::integrate on %d of them: largest absolute difference %.2g.",
                sum(!is.na(st$integrate_check)), max(abs(st$integrate_check), na.rm = TRUE)), "",
        table_md(do.call(rbind, lapply(split(st, st$regime), function(x) {
          cbind(regime = x$regime[[1L]], states = nrow(x),
                as.data.frame(lapply(stats::setNames(rules, rules), function(r) max(abs(x[[r]])))))
        }))), "")
  }
  writeLines(lines, file.path(mm_out, "report.md"))
  cat(lines, sep = "\n")
}

mm_rule_nodes <- function(rule) {
  switch(rule$type, gh = rule$nodes, mixture = 2L * rule$nodes,
         grid = length(seq(-rule$width, rule$width, by = rule$step)), NA_integer_)
}
