#!/usr/bin/env Rscript
# Offline evaluation of an automatic choice among the three series methods
# (user request of 2026-11-15). Uses only the existing fits of a build on
# e1-e3: a rule picks, per dataset, one of the methods already fitted, and is
# scored as that method's fit. Nothing is refitted.
#
#   Rscript validation/sparse_joint/series_methods/auto_eval.R --build=cv10 --pkg=PATH [--workers=6]
#
# --pkg: an alfak2 source tree to load (devtools::load_all) for the domain
# functions; they are those of the fits.
#
# Fixed before any rule was scored:
#
# Indicators, from the counts alone with the fits' settings (major frequency
# 0.01, major count 3, max states 5000, copy numbers 1-8), plus one from the
# anchored fit:
#   ends_major   least share, over the time points in between, of the cells
#                at a time point on the majors of the first and last time
#                points (0: anchored leaves that time point out);
#   ends_domain  the same for cells in the anchored domain (those majors and
#                their one-step neighbours);
#   continuity   least share, over consecutive pairs, of the cells at the later
#                time point within one step of the earlier time point's own
#                majors;
#   weight       the anchored fit's intermediate weight, 1 / max(1,
#                dispersion of the majors' counts): low with strong drift.
# Rules:
#   fixed        ends, anchored, pairwise, state_space;
#   switch       pairwise (or state_space) when an indicator is below a
#                threshold, anchored otherwise;
#   switch+drift below the threshold, pairwise when weight >= w, otherwise
#                state_space; anchored above it.
# A rule is acceptable when
#   e1  in at least 18 of the 24 cells the median RMSE ratio to `ends` is at
#       most 1.03 at every distance, and the median over cells at most 1.02
#       at every distance (the gate of report.R);
#   e3  its median ratio at every distance is at most anchored's + 0.01: on
#       neutral data every fitness difference is an error, and the rule must
#       not add them.
# Among acceptable rules the choice is the lowest e2 ratio (median over
# distances 0-2 of the set's medians); within 0.01 of it the simpler one
# (one indicator before two; pairwise before state_space; the coarser
# threshold). Scores use datasets where every method has a fit (the
# state-space fits that hit the time limit are left out, and counted).
# The thresholds are tuned on e1-e3; the paper benchmark (E4) is the held-out
# check. As a check on stability, the chosen rule's threshold is also chosen
# on landscapes 1-6 alone and scored on 7-12, and the other way round.

sm_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
})
source(file.path(sm_script_dir, "setup.R"))

ae_arg <- function(name, default = NULL) {
  hit <- grep(sprintf("^--%s=", name), commandArgs(TRUE), value = TRUE)
  if (length(hit)) sub(sprintf("^--%s=", name), "", hit[[1L]]) else default
}
build <- ae_arg("build", "cv10")
workers <- as.integer(ae_arg("workers", "6"))
suppressMessages(devtools::load_all(ae_arg("pkg"), quiet = TRUE, export_all = TRUE))
out_dir <- file.path(sm_default_output, build, "auto")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
settings <- list(major_frequency = 0.01, major_count = 3, max_states = 5000, min_cn = 1L, max_cn = 8L)

ae_share_in <- function(counts, k, labels) sum(counts[rownames(counts) %in% labels, k]) / sum(counts[, k])

ae_indicators <- function(counts) {
  s <- settings
  K <- ncol(counts)
  ends <- sparse_joint_domain(counts, s$major_frequency, s$major_count, s$max_states, s$min_cn, s$max_cn, c(1L, K))
  major <- ends$labels[seq_len(ends$n_major)]
  middle <- seq(2L, K - 1L)
  ends_major <- min(vapply(middle, function(k) ae_share_in(counts, k, major), numeric(1)))
  ends_domain <- min(vapply(middle, function(k) ae_share_in(counts, k, ends$labels), numeric(1)))
  continuity <- min(vapply(seq_len(K - 1L), function(k) {
    y <- counts[, c(k, k + 1L), drop = FALSE]
    y <- y[rowSums(y) > 0, , drop = FALSE]
    d <- sparse_joint_domain(y, s$major_frequency, s$major_count, s$max_states, s$min_cn, s$max_cn, 1L)
    ae_share_in(y, 2L, d$labels)
  }, numeric(1)))
  c(ends_major = ends_major, ends_domain = ends_domain, continuity = continuity)
}

# ---- Indicators and per-dataset errors ---------------------------------------

methods <- c("ends", "anchored", "pairwise", "state_space")
rows <- list()
for (set in sm_sets) {
  config <- sm_config(set, build)
  sm_link(config)
  index <- utils::read.csv(si_index_path(config), stringsAsFactors = FALSE)
  ind <- si_parallel(index$dataset_id, function(id) {
    dataset <- readRDS(stt_dataset_path(config, id))
    c(list(dataset_id = id), as.list(ae_indicators(mt_counts(dataset, "all"))))
  }, workers)
  ind <- do.call(rbind, lapply(ind, as.data.frame))
  fits <- utils::read.csv(file.path(config$output_dir, "metrics_fit.csv"), stringsAsFactors = FALSE)
  weight <- fits[fits$method == "anchored", c("dataset_id", "intermediate_weight")]
  names(weight)[2] <- "weight"
  status <- stats::reshape(fits[fits$method %in% methods, c("dataset_id", "method", "status")],
                           idvar = "dataset_id", timevar = "method", direction = "wide")
  d <- utils::read.csv(file.path(config$output_dir, "metrics_distance.csv"), stringsAsFactors = FALSE)
  d <- d[d$stratum %in% c("d0", "d1", "d2") & d$method %in% methods & d$status == "ok",
         c("dataset_id", "method", "stratum", "rmse")]
  wide <- stats::reshape(d, idvar = c("dataset_id", "stratum"), timevar = "method", direction = "wide")
  cell <- c("landscape_id", "bottleneck", "spacing", "depth", "samples", if (set == "e2") "start")
  x <- merge(merge(merge(index[c("dataset_id", intersect(cell, names(index)))], ind), weight), wide)
  x <- merge(x, status, by = "dataset_id")
  x$set <- set
  rows[[set]] <- x
}
all <- do.call(rbind, lapply(rows, function(x) { x$start <- if (is.null(x$start)) NA else x$start; x }))
complete <- stats::complete.cases(all[paste0("rmse.", methods)])
missing <- stats::aggregate(cbind(n = 1, missing = !complete) ~ set, data.frame(set = all$set, complete), sum)
utils::write.csv(all, file.path(out_dir, "indicators.csv"), row.names = FALSE)
x <- all[complete, ]

# ---- Rules ------------------------------------------------------------------

thresholds <- c(1e-9, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)
rules <- list(list(name = "ends", pick = function(x) rep("ends", nrow(x)), complexity = 0),
              list(name = "anchored", pick = function(x) rep("anchored", nrow(x)), complexity = 0),
              list(name = "pairwise", pick = function(x) rep("pairwise", nrow(x)), complexity = 0),
              list(name = "state_space", pick = function(x) rep("state_space", nrow(x)), complexity = 0))
for (ind in c("ends_major", "ends_domain", "continuity")) for (c0 in thresholds) {
  for (m in c("pairwise", "state_space")) {
    rules[[length(rules) + 1L]] <- local({
      ind <- ind; c0 <- c0; m <- m
      list(name = sprintf("%s<%g->%s", ind, c0, m), indicator = ind, threshold = c0, target = m,
           pick = function(x) ifelse(x[[ind]] < c0, m, "anchored"), complexity = 1)
    })
  }
  for (w in c(0.5, 0.8, 0.95)) {
    rules[[length(rules) + 1L]] <- local({
      ind <- ind; c0 <- c0; w <- w
      list(name = sprintf("%s<%g->pairwise|weight>=%g else state_space", ind, c0, w), indicator = ind,
           threshold = c0, target = "both",
           pick = function(x) ifelse(x[[ind]] < c0, ifelse(x$weight >= w, "pairwise", "state_space"), "anchored"),
           complexity = 2)
    })
  }
}

ae_score <- function(x, pick) {
  rmse <- x[cbind(seq_len(nrow(x)), match(paste0("rmse.", pick), names(x)))]
  x$ratio <- as.numeric(rmse) / x$rmse.ends
  x
}

ae_evaluate <- function(x, rule, anchored_e3) {
  pick <- rule$pick(x)
  s <- ae_score(x, pick)
  med <- stats::aggregate(ratio ~ set + stratum, s, stats::median)
  e1 <- s[s$set == "e1", ]
  cells <- stats::aggregate(ratio ~ bottleneck + spacing + depth + samples + stratum, e1, stats::median)
  per_cell <- stats::aggregate(ratio ~ bottleneck + spacing + depth + samples, cells, max)
  over <- stats::aggregate(ratio ~ stratum, cells, stats::median)
  e3 <- med[med$set == "e3", ]
  e3_ok <- all(e3$ratio[order(e3$stratum)] <= anchored_e3[order(names(anchored_e3))] + 0.01)
  e2 <- med$ratio[med$set == "e2"]
  share <- function(set) {
    p <- pick[x$set == set & x$stratum == "d0"]
    sprintf("a %.0f%% p %.0f%% s %.0f%%", 100 * mean(p == "anchored"), 100 * mean(p == "pairwise"),
            100 * mean(p == "state_space"))
  }
  data.frame(rule = rule$name, complexity = rule$complexity,
             e1_cells = sum(per_cell$ratio <= 1.03), e1_over = max(over$ratio),
             e1_pass = sum(per_cell$ratio <= 1.03) >= 18 && all(over$ratio <= 1.02),
             e3_max_excess = max(e3$ratio - anchored_e3[e3$stratum]), e3_pass = e3_ok,
             e2_ratio = stats::median(e2),
             e1 = paste(sprintf("%.3f", med$ratio[med$set == "e1"]), collapse = "/"),
             e2 = paste(sprintf("%.3f", e2), collapse = "/"),
             e3 = paste(sprintf("%.3f", e3$ratio), collapse = "/"),
             pick_e1 = share("e1"), pick_e2 = share("e2"), pick_e3 = share("e3"),
             stringsAsFactors = FALSE)
}

ae_table <- function(x) {
  a <- ae_score(x, rep("anchored", nrow(x)))
  anchored_e3 <- with(stats::aggregate(ratio ~ stratum, a[a$set == "e3", ], stats::median),
                      stats::setNames(ratio, stratum))
  out <- do.call(rbind, lapply(rules, ae_evaluate, x = x, anchored_e3 = anchored_e3))
  out$acceptable <- out$e1_pass & out$e3_pass
  out
}

ae_choose <- function(table) {
  ok <- table[table$acceptable, ]
  # Nothing acceptable (possible on half the landscapes): the default stays.
  if (!nrow(ok)) return(table[table$rule == "anchored", ])
  best <- min(ok$e2_ratio)
  near <- ok[ok$e2_ratio <= best + 0.01, ]
  # Simpler: fewer indicators, then pairwise before state_space, then the
  # coarser (smaller) threshold. Rule names carry the target.
  target <- ifelse(grepl("->pairwise$", near$rule), 1L, ifelse(grepl("->state_space$", near$rule), 2L, 3L))
  threshold <- suppressWarnings(as.numeric(sub(".*<([0-9.e-]+)->.*", "\\1", near$rule)))
  near[order(near$complexity, target, threshold, near$e2_ratio), ][1L, ]
}

table <- ae_table(x)
utils::write.csv(table, file.path(out_dir, "rules.csv"), row.names = FALSE)
chosen <- ae_choose(table)

# Oracle: per dataset the method with the lowest mean ratio over distances.
o <- stats::aggregate(cbind(rmse.ends, rmse.anchored, rmse.pairwise, rmse.state_space) ~ set + dataset_id, x, mean)
o$best <- methods[apply(o[paste0("rmse.", methods)][, -1L], 1L, which.min) + 1L]
oracle_pick <- o$best[match(paste(x$set, x$dataset_id), paste(o$set, o$dataset_id))]
oracle <- ae_score(x, oracle_pick)
oracle_med <- stats::aggregate(ratio ~ set + stratum, oracle, stats::median)

# Stability: choose on one half of the landscapes, score on the other.
halves <- list(first = sprintf("mt_L%02d", 1:6), second = sprintf("mt_L%02d", 7:12))
stability <- do.call(rbind, lapply(names(halves), function(h) {
  train <- ae_table(x[x$landscape_id %in% halves[[h]], ])
  pick <- ae_choose(train)
  test <- ae_table(x[!x$landscape_id %in% halves[[h]], ])
  keep <- c("rule", "e1_cells", "e1_pass", "e3_max_excess", "e3_pass", "e2_ratio", "e1", "e2", "e3")
  cbind(data.frame(tuned_on = h, train_acceptable = sum(train$acceptable),
                   train_anchored_e1_cells = train$e1_cells[train$rule == "anchored"], stringsAsFactors = FALSE),
        rbind(test[test$rule == pick$rule, keep], test[test$rule == "anchored", keep],
              test[test$rule == chosen$rule, keep]))
}))

sink(file.path(out_dir, "report.md"))
cat("# Automatic choice of the series method, build", build, "\n\n")
cat("Datasets with every method fitted (state-space time-outs left out):\n\n")
print(knitr::kable(missing))
cat("\nFixed methods and the acceptable rules (ratios d0/d1/d2):\n\n")
show <- table[table$complexity == 0 | table$acceptable, ]
show <- show[order(-show$acceptable, show$e2_ratio), ]
print(knitr::kable(show[c("rule", "acceptable", "e1_cells", "e3_max_excess", "e2_ratio", "e1", "e2", "e3",
                          "pick_e1", "pick_e2", "pick_e3")], digits = 3))
cat("\nChosen rule:", chosen$rule, "\n\n")
cat("Oracle (best fitted method per dataset), median ratio d0/d1/d2:\n\n")
print(knitr::kable(stats::reshape(oracle_med, idvar = "set", timevar = "stratum", direction = "wide"), digits = 3))
cat("\nStability (threshold chosen on one half of the landscapes, scored on the other):\n\n")
print(knitr::kable(stability, digits = 3))
cat("\nIndicator by set (median, 10% quantile):\n\n")
q <- stats::aggregate(cbind(ends_major, ends_domain, continuity, weight) ~ set + samples, all,
                      function(v) sprintf("%.2f (%.2f)", stats::median(v), stats::quantile(v, 0.1)))
# The `samples` column counts time points; the header says so.
print(knitr::kable(q, col.names = sub("^samples$", "time points", names(q))))
sink()
cat(readLines(file.path(out_dir, "report.md")), sep = "\n")
