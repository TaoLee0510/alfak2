#!/usr/bin/env Rscript
# Stage C, step 1: tables of the convergence diagnosis written by diagnose.R.
#
#   Rscript validation/sparse_joint/convergence/summarise.R [--build=ff13cf4] [--input-dir=DIR]
#
# One row per candidate history (results/<build>/diagnose/candidates.csv) and
# a report (report.md). Every candidate, converged or not by nlminb, is put in
# one class, by the first of these that applies:
#
#   degenerate      the inner Hessian where the fit ended has an eigenvalue
#                   below `eig_tol` (the Laplace objective is not usable there)
#   path_dependent  the objective from a fixed inner start differs from the
#                   fit's by more than `f_tol`
#   rough           the objective from a fixed inner start is not smooth:
#                   a finite-difference step of 1e-3 predicts a decrease that
#                   steps of 0.05 do not find
#   not_stationary  a coordinate step of 0.05 lowers the objective by more
#                   than `probe_tol`, or the projected finite-difference
#                   gradient exceeds `grad_tol`
#   optimum         none of these
#
# The tolerances are provisional; the report shows the distributions they are
# chosen from (smallest inner eigenvalue of converged and non-converged fits).

cv_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
args <- commandArgs(TRUE)
build <- if (any(grepl("^--build=", args))) sub("^--build=", "", grep("^--build=", args, value = TRUE)[[1L]]) else "ff13cf4"
input <- if (any(grepl("^--input-dir=", args))) {
  sub("^--input-dir=", "", grep("^--input-dir=", args, value = TRUE)[[1L]])
} else file.path(cv_script_dir, "results", build, "diagnose")

eig_tol <- 2e-3
f_tol <- 1e-3
grad_tol <- 0.05
probe_tol <- 1e-3

files <- list.files(input, pattern = "^dev_.*[.]rds$", full.names = TRUE)
records <- lapply(files, readRDS)
failed <- Filter(function(r) !is.null(r$error), records)
records <- Filter(function(r) is.null(r$error), records)

num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x[[1L]])
maxabs <- function(x) if (is.null(x) || all(is.na(x))) NA_real_ else max(abs(x), na.rm = TRUE)

rows <- do.call(rbind, lapply(records, function(r) {
  do.call(rbind, lapply(r$candidates, function(c) {
    dir <- c$warm$direction
    top <- if (is.data.frame(dir) && nrow(dir)) dir[1L, ] else NULL
    data.frame(
      id = r$id,
      scenario = sub("^dev_L[0-9]+_P[0-9]+_", "", r$id),
      landscape = sub("^dev_(L[0-9]+)_.*$", "\\1", r$id),
      history = c$history, selected = isTRUE(all.equal(c$history, r$selected$history)),
      n_states = c$n_states, n_random = c$n_random, n_free = length(c$par),
      code = num(c$nlminb$convergence), message = if (is.null(c$nlminb$message)) NA_character_ else c$nlminb$message,
      iterations = num(c$nlminb$iterations), evaluations = num(c$nlminb$evaluations),
      at_bound = any(c$at_lower) || any(c$at_upper),
      f_nlminb = num(c$nlminb$objective), f_warm = num(c$warm$f), f_fresh = num(c$fresh$f),
      eig_warm = num(c$warm$min_eig), eig2_warm = num(c$warm$second_eig),
      eig_fresh = num(c$fresh$min_eig), eig_start = num(c$start$min_eig),
      joint_warm = num(c$warm$joint), joint_fresh = num(c$fresh$joint),
      logdet_warm = num(c$warm$half_logdet), logdet_fresh = num(c$fresh$half_logdet),
      inner_grad_warm = num(c$warm$inner_grad), inner_grad_fresh = num(c$fresh$inner_grad),
      tmb_grad_warm = maxabs(c$warm$tmb_grad), tmb_grad_fresh = maxabs(c$fresh$tmb_grad),
      fd_grad = maxabs(c$fresh$fd_grad), projected_grad = maxabs(c$fresh$projected_grad),
      probe_decrease = num(c$fresh$probe_decrease),
      flat_variable = if (is.null(top)) NA_character_ else top$variable,
      flat_major = if (is.null(top)) NA else top$major,
      flat_loading = if (is.null(top)) NA_real_ else abs(top$loading),
      flat_y0 = if (is.null(top)) NA_real_ else top$y0, flat_y1 = if (is.null(top)) NA_real_ else top$y1,
      probe_seconds = num(c$probe_seconds),
      stringsAsFactors = FALSE
    )
  }))
}))
rows$converged <- rows$code == 0
rows$f_gap <- rows$f_warm - rows$f_fresh
rows$eig_min <- pmin(rows$eig_warm, rows$eig_fresh)
# Decrease that the step-1e-3 gradient predicts for a step of 0.05 against it,
# against the decrease the steps of 0.05 found.
rows$predicted_decrease <- 0.05 * rows$projected_grad
rows$class <- with(rows, ifelse(
  !is.finite(eig_warm) | eig_warm < eig_tol, "degenerate",
  ifelse(!is.finite(f_gap) | abs(f_gap) > f_tol, "path_dependent",
         ifelse(projected_grad > grad_tol & probe_decrease < 0.5 * predicted_decrease & probe_decrease <= probe_tol,
                "rough",
                ifelse(probe_decrease > probe_tol | projected_grad > grad_tol, "not_stationary", "optimum")))))
utils::write.csv(rows, file.path(input, "candidates.csv"), row.names = FALSE)

lines <- character()
add <- function(...) lines <<- c(lines, ...)
table_md <- function(x, digits = 3L) {
  x[] <- lapply(x, function(v) if (is.numeric(v) && any(v != round(v), na.rm = TRUE)) formatC(v, digits = digits, format = "g") else v)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|"),
    apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
count_table <- function(frame, by) {
  t <- as.data.frame.matrix(table(frame[[by]], factor(frame$class, c("optimum", "not_stationary", "rough",
                                                                          "path_dependent", "degenerate"))))
  cbind(stats::setNames(data.frame(rownames(t), stringsAsFactors = FALSE), by), t, total = rowSums(t))
}

add(sprintf("# Convergence diagnosis, build %s", build), "",
    sprintf("Generated %s from %d items (%d failed), %d candidate histories.", format(Sys.time(), "%Y-%m-%d %H:%M"),
            length(files), length(failed), nrow(rows)), "",
    sprintf("Provisional tolerances: eig_tol %g, f_tol %g, grad_tol %g, probe_tol %g.", eig_tol, f_tol, grad_tol,
            probe_tol), "")
if (length(failed)) add("Failed items:", "", paste0("- ", vapply(failed, function(r) paste(r$id, r$error), "")), "")

add("## Class by nlminb's verdict (all candidates)", "",
    table_md(count_table(transform(rows, nlminb = ifelse(converged, "converged", "not converged")), "nlminb")), "")
add("## Class by nlminb's verdict (selected histories)", "",
    table_md(count_table(transform(rows[rows$selected, ], nlminb = ifelse(converged, "converged", "not converged")),
                         "nlminb")), "")
add("## Class by scenario (all candidates)", "", table_md(count_table(rows, "scenario")), "")
add("## Class by history (all candidates)", "", table_md(count_table(rows, "history")), "")

bins <- c(-Inf, -10, -8, -6, -5, -4, -3, -2, -1.5, -1, Inf)
eig_table <- function(v) table(cut(log10(pmax(v, 1e-300)), bins))
eig <- data.frame(log10_eig = names(eig_table(rows$eig_warm)),
                  converged = as.integer(eig_table(rows$eig_warm[rows$converged])),
                  not_converged = as.integer(eig_table(rows$eig_warm[!rows$converged])),
                  fresh_converged = as.integer(eig_table(rows$eig_fresh[rows$converged])),
                  fresh_not_converged = as.integer(eig_table(rows$eig_fresh[!rows$converged])),
                  start = as.integer(eig_table(rows$eig_start)))
add("## Smallest inner eigenvalue where the fit ended", "",
    "`converged`/`not_converged`: warm (as the fit saw it); `fresh_*`: inner start fixed at the initial values;",
    "`start`: at the prior-centred scales.", "", table_md(eig, 0L), "")

ratio <- data.frame(statistic = c("warm eig / start eig", "fresh eig / start eig"),
                    converged = c(stats::median(rows$eig_warm[rows$converged] / rows$eig_start[rows$converged]),
                                  stats::median(rows$eig_fresh[rows$converged] / rows$eig_start[rows$converged])),
                    not_converged = c(stats::median(rows$eig_warm[!rows$converged] / rows$eig_start[!rows$converged]),
                                      stats::median(rows$eig_fresh[!rows$converged] / rows$eig_start[!rows$converged])))
add("Median ratios:", "", table_md(ratio), "")

gap <- stats::aggregate(cbind(f_gap, joint_warm - joint_fresh, logdet_warm - logdet_fresh) ~ class, rows,
                        stats::median)
names(gap) <- c("class", "f warm - fresh", "joint warm - fresh", "half logdet warm - fresh")
add("## Where the objective gap comes from (medians)", "", table_md(gap), "")

flat <- rows[rows$class == "degenerate", ]
if (nrow(flat)) {
  add("## Flattest direction of degenerate candidates", "",
      table_md(stats::aggregate(id ~ flat_variable + flat_major + flat_y0 + flat_y1, flat, length)), "")
}

time <- data.frame(statistic = c("median probe seconds", "median fit seconds per item"),
                   value = c(stats::median(rows$probe_seconds),
                             stats::median(vapply(records, function(r) r$seconds, numeric(1)))))
add("## Cost", "", table_md(time), "")

check <- vapply(records, function(r) abs(r$selected$objective - r$reference_objective), numeric(1))
add(sprintf("Selected objective against the stage-0 record: largest absolute difference %.3g.", max(check)), "")

writeLines(lines, file.path(input, "report.md"))
cat("Wrote", file.path(input, "report.md"), "and candidates.csv\n")
