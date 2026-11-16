#!/usr/bin/env Rscript
# Stage C, step 3: the convergence fix against its baseline on the stage-0
# reference set (handoff HANDOFF_CLAUDE_2026-11-13_multi_timepoint.md).
#
#   Rscript validation/sparse_joint/convergence/evaluate.R --builds=ff13cf4,cv1 \
#     [--reference-dir=DIR] [--output-dir=DIR]
#
# Reads the records reference.R wrote for both builds (by default
# multi_timepoint/results/reference; records made on HPC are synced to a
# directory of their own and given with --reference-dir). Gates:
#   - fits accepted at the first try are identical to the baseline (largest
#     relative difference at most 1e-8, convergence fields aside). The fix
#     also changes which candidate a fit keeps (only accepted ones), and
#     restarted candidates report other objectives; so where the kept history
#     (or a cohort's per-dataset histories or trend) is the baseline's, the
#     record is compared without the objectives of the other candidates, and
#     where it changed, only the kept candidate's objective is compared with
#     the baseline's objective for that candidate;
#   - at least 95% of candidate histories are accepted;
#   - fits that changed: median over items of the RMSE ratio (new / baseline)
#     at distances 0, 1 and 2 at most 1.02, median 95% coverage at least 0.90;
#   - total time at most twice the baseline's.
# Writes items.csv, scores.csv and report.md.

cv_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
mt_script_dir <- normalizePath(file.path(cv_script_dir, "..", "multi_timepoint"))
source(file.path(mt_script_dir, "reference.R"))

args <- list()
for (arg in commandArgs(trailingOnly = TRUE)) {
  body <- sub("^--", "", arg)
  args[[gsub("-", "_", sub("=.*$", "", body))]] <- if (grepl("=", body, fixed = TRUE)) sub("^[^=]*=", "", body) else TRUE
}
builds <- strsplit(if (is.null(args$builds)) "ff13cf4,cv1" else args$builds, ",", fixed = TRUE)[[1L]]
stopifnot(length(builds) == 2L)
reference_dir <- if (is.null(args$reference_dir)) file.path(mt_results, "reference") else args$reference_dir
output <- if (is.null(args$output_dir)) file.path(cv_script_dir, "results", paste(builds, collapse = "_"), "evaluate") else
  args$output_dir
dir.create(output, recursive = TRUE, showWarnings = FALSE)
parts <- c("dev", "entry", "free", "rates", "cohort", "compat")
# What the fix is allowed to change in a fit it accepted at the first try.
convergence_keys <- paste0("/seconds$|/diagnostics/(convergence|message|optimisation|trend_convergence)(/|$)|",
                           "/history_profiles?/.*convergence|/history_profile/convergence")
# Objectives of the candidates a fit did not keep.
profile_keys <- "/history_profiles?/.*objective$|/history_profile/objective$|/diagnostics/trend_objective$"

record_of <- function(build, part, item) {
  path <- file.path(reference_dir, build, part, paste0(item, ".rds"))
  if (file.exists(path)) readRDS(path) else NULL
}

# Candidate histories of a record and whether each was accepted; the inner
# fits of a rates record. Exact names: `$` would match `history` to
# `history_profiles`.
candidates_of <- function(record) {
  if (!is.null(record[["fits"]])) return(unlist(lapply(record[["fits"]], candidates_of)))
  d <- record[["diagnostics"]]
  if (!is.null(d[["history_profile"]])) return(d[["history_profile"]][["convergence"]] == 0)
  if (!is.null(d[["history_profiles"]])) {
    return(unlist(lapply(d[["history_profiles"]], function(p) p[["convergence"]] == 0)))
  }
  if (!is.null(d[["convergence"]])) d[["convergence"]] == 0 else logical()
}
convergence_of <- function(record) {
  if (!is.null(record[["fits"]])) {
    return(max(vapply(record[["fits"]], convergence_of, numeric(1))))
  }
  value <- record[["diagnostics"]][["convergence"]]
  if (is.null(value)) NA_real_ else as.numeric(value)
}
history_of <- function(record) {
  value <- record[["diagnostics"]][["history"]]
  if (is.null(value)) NA_real_ else as.numeric(value[[1L]])
}
# What a record kept: its history, or a cohort's trend and per-dataset
# histories (cohort records do not store those, so they are chosen again from
# the profiles by each build's rule: the baseline's lowest objective, the
# fix's lowest accepted one).
lowest <- function(objective, accepted = rep(TRUE, length(objective))) {
  usable <- is.finite(objective) & accepted
  if (!any(usable)) usable <- is.finite(objective)
  if (!any(usable)) NA_integer_ else which(usable)[which.min(objective[usable])]
}
selection_of <- function(record, fixed) {
  d <- record[["diagnostics"]]
  if (is.null(d)) return(NA_character_)
  if (!is.null(d[["history_profiles"]])) {
    histories <- vapply(d[["history_profiles"]], function(p) {
      p[["history"]][[lowest(p[["objective"]], if (fixed) p[["convergence"]] == 0 else TRUE)]]
    }, numeric(1))
    return(paste(c(d[["trend"]], histories), collapse = ","))
  }
  paste(c(d[["trend"]], d[["history"]]), collapse = ",")
}
# The baseline's objective for the history the new record kept.
baseline_objective_of <- function(a, b) {
  profile <- a[["diagnostics"]][["history_profile"]]
  history <- history_of(b)
  if (is.null(profile) || is.na(history)) return(NA_real_)
  value <- profile[["objective"]][profile[["history"]] == history]
  if (length(value) == 1L) value else NA_real_
}

items <- list()
for (part in parts) {
  files <- intersect(list.files(file.path(reference_dir, builds[[1L]], part)),
                     list.files(file.path(reference_dir, builds[[2L]], part)))
  for (file in files) {
    item <- sub("[.]rds$", "", file)
    a <- record_of(builds[[1L]], part, item)
    b <- record_of(builds[[2L]], part, item)
    errored <- !is.null(a$error) || !is.null(b$error)
    comparison <- if (errored) NULL else mt_compare_records(a, b, ignore = convergence_keys)
    same_selection <- !errored && identical(selection_of(a, FALSE), selection_of(b, TRUE))
    kept <- if (same_selection) mt_compare_records(a, b, ignore = paste(convergence_keys, profile_keys, sep = "|"))
    kept_objective <- if (errored) NA_real_ else baseline_objective_of(a, b)
    check <- b[["diagnostics"]][["optimisation"]]
    accepted <- candidates_of(b)
    items[[length(items) + 1L]] <- data.frame(
      part = part, item = item,
      error_baseline = if (is.null(a$error)) "" else a$error, error_new = if (is.null(b$error)) "" else b$error,
      identical = !errored && comparison$worst <= 1e-8 && !length(comparison$mismatched) &&
        !length(setdiff(comparison$missing, grep(convergence_keys, comparison$missing, value = TRUE))),
      worst_relative = if (errored) NA_real_ else comparison$worst,
      same_selection = same_selection,
      kept_identical = if (errored) NA else if (same_selection) {
        kept$worst <= 1e-8 && !length(kept$mismatched) &&
          !length(setdiff(kept$missing, grep(paste(convergence_keys, profile_keys, sep = "|"), kept$missing, value = TRUE)))
      } else if (is.finite(kept_objective)) {
        abs(b[["diagnostics"]][["objective"]] - kept_objective) <= 1e-8 * max(1, abs(kept_objective))
      } else NA,
      history_baseline = history_of(a), history_new = history_of(b),
      convergence_baseline = if (errored) NA_real_ else convergence_of(a),
      convergence_new = if (errored) NA_real_ else convergence_of(b),
      stage = if (is.null(check)) NA_character_ else check[["stage"]],
      reasons = if (is.null(check)) NA_character_ else paste(check[["reasons"]], collapse = "; "),
      candidates = length(accepted), candidates_accepted = sum(accepted),
      seconds_baseline = if (is.null(a$seconds)) NA_real_ else a$seconds,
      seconds_new = if (is.null(b$seconds)) NA_real_ else b$seconds,
      stringsAsFactors = FALSE
    )
  }
}
items <- do.call(rbind, items)
utils::write.csv(items, file.path(output, "items.csv"), row.names = FALSE)

# Scores of the development items against the truth.
as_prediction <- function(p) {
  data.frame(karyotype = p$karyotype, fitness = p$relative_fitness, fitness_sd = p$relative_sd,
             stringsAsFactors = FALSE)
}
dev <- items[items$part == "dev" & !nzchar(items$error_baseline) & !nzchar(items$error_new), ]
scores <- do.call(rbind, lapply(dev$item, function(item) {
  nodes <- mt_nodes(item)
  truth <- stt_true_growth_rate(nodes$karyotype, mt_dataset(item)$truth)
  do.call(rbind, lapply(builds, function(build) {
    record <- record_of(build, "dev", item)
    m <- stt_score_predictions(nodes, truth, as_prediction(record$prediction), 5)$metrics
    cbind(data.frame(item = item, build = build, stringsAsFactors = FALSE),
          m[m$stratum %in% c("d0", "d1", "d2"), c("stratum", "rmse", "bias", "spearman", "coverage95", "skill")])
  }))
}))
utils::write.csv(scores, file.path(output, "scores.csv"), row.names = FALSE)

lines <- character()
add <- function(...) lines <<- c(lines, ...)
table_md <- function(x, digits = 3L) {
  x[] <- lapply(x, function(v) if (is.numeric(v) && any(v != round(v), na.rm = TRUE)) formatC(v, digits = digits, format = "g") else v)
  c(paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("|", paste(rep(" --- ", ncol(x)), collapse = "|"), "|"),
    apply(x, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
}
add(sprintf("# Convergence fix: %s against %s", builds[[2L]], builds[[1L]]), "",
    sprintf("Generated %s from %s.", format(Sys.time(), "%Y-%m-%d %H:%M"), reference_dir), "")

ok <- !nzchar(items$error_baseline) & !nzchar(items$error_new)
first <- ok & items$stage %in% "first"
by_part <- do.call(rbind, lapply(split(items, items$part), function(x) {
  data.frame(part = x$part[[1L]], items = nrow(x), errors_baseline = sum(nzchar(x$error_baseline)),
             errors_new = sum(nzchar(x$error_new)), first_try = sum(x$stage %in% "first"),
             first_try_identical = sum(x$stage %in% "first" & x$kept_identical %in% TRUE),
             kept_changed = sum(!x$same_selection),
             converged_baseline = sum(x$convergence_baseline == 0, na.rm = TRUE),
             converged_new = sum(x$convergence_new == 0, na.rm = TRUE),
             candidates = sum(x$candidates), accepted = sum(x$candidates_accepted),
             seconds_baseline = sum(x$seconds_baseline, na.rm = TRUE), seconds_new = sum(x$seconds_new, na.rm = TRUE))
}))
add("## Items", "", table_md(by_part, 4L), "")
errors <- items[nzchar(items$error_baseline) | nzchar(items$error_new), c("part", "item", "error_baseline", "error_new")]
if (nrow(errors)) add("Errors:", "", table_md(errors), "")

changed <- dev$item[!dev$identical]
ratio_table <- function(set) {
  if (!length(set)) return(NULL)
  s <- scores[scores$item %in% set, ]
  w <- merge(s[s$build == builds[[1L]], ], s[s$build == builds[[2L]], ], by = c("item", "stratum"),
             suffixes = c("_baseline", "_new"))
  w$ratio <- w$rmse_new / w$rmse_baseline
  do.call(rbind, lapply(split(w, w$stratum), function(x) {
    data.frame(stratum = x$stratum[[1L]], items = nrow(x), median_rmse_ratio = stats::median(x$ratio),
               better = sum(x$ratio < 1), worse = sum(x$ratio > 1), max_ratio = max(x$ratio),
               coverage_baseline = stats::median(x$coverage95_baseline), coverage_new = stats::median(x$coverage95_new))
  }))
}
changed_table <- ratio_table(changed)
if (!is.null(changed_table)) add(sprintf("## Development items that changed (%d)", length(changed)), "",
                                 table_md(changed_table), "")
all_table <- ratio_table(dev$item)
add("## All development items", "", table_md(all_table), "")
history_moves <- dev[!is.na(dev$history_baseline) & dev$history_baseline != dev$history_new,
                     c("item", "history_baseline", "history_new", "convergence_baseline", "convergence_new", "reasons")]
if (nrow(history_moves)) add("History changed:", "", table_md(history_moves), "")
not_accepted <- items[ok & items$convergence_new %in% 1, c("part", "item", "stage", "reasons")]
if (nrow(not_accepted)) add("Not accepted (new build):", "", table_md(not_accepted), "")

gate <- data.frame(
  gate = c("first-try fits identical", "candidates accepted >= 95%", "changed fits: median RMSE ratio <= 1.02",
           "changed fits: median coverage >= 0.90", "time <= 2x"),
  value = c(sprintf("%d / %d (kept candidate unchanged: %d / %d; changed: objective identical %d / %d; not checkable %d)",
                    sum(first & items$kept_identical %in% TRUE), sum(first),
                    sum(first & items$same_selection & items$kept_identical %in% TRUE), sum(first & items$same_selection),
                    sum(first & !items$same_selection & items$kept_identical %in% TRUE),
                    sum(first & !items$same_selection), sum(first & is.na(items$kept_identical))),
            sprintf("%.1f%%", 100 * sum(items$candidates_accepted[ok]) / sum(items$candidates[ok])),
            if (is.null(changed_table)) "none changed" else paste(sprintf("%.3f", changed_table$median_rmse_ratio), collapse = " / "),
            if (is.null(changed_table)) "none changed" else paste(sprintf("%.3f", changed_table$coverage_new), collapse = " / "),
            sprintf("%.2f", sum(items$seconds_new[ok]) / sum(items$seconds_baseline[ok]))),
  pass = c(!any(items$kept_identical[first] %in% FALSE),
           sum(items$candidates_accepted[ok]) / sum(items$candidates[ok]) >= 0.95,
           is.null(changed_table) || all(changed_table$median_rmse_ratio <= 1.02),
           is.null(changed_table) || all(changed_table$coverage_new >= 0.90),
           sum(items$seconds_new[ok]) / sum(items$seconds_baseline[ok]) <= 2),
  stringsAsFactors = FALSE
)
add("## Gates", "", table_md(gate), "")
writeLines(lines, file.path(output, "report.md"))
cat(lines, sep = "\n")
