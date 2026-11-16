#!/usr/bin/env Rscript
# Fig. 4d check, step 1: where does alfak2's large in vitro |delta f| come from?
#
# Reads the Salehi fits of the paper benchmark round 05a8855 (copied to
# results/input/ from the HPC results; read only). Selects fits as the
# benchmark's Fig. 4 does (scripts/fig4_cna_effects.R): xv > 0, per trajectory
# the minobs with the highest xv, the lowest of a tie. For each frequent
# karyotype (landscape fq) and each of its 44 one-chromosome neighbours,
# delta f = f(frequent) - f(neighbour) from landscape.Rds, as in the paper.
# Each neighbour is classed by where alfak2's value comes from:
#   observed_D   in the model domain and seen at some time point
#   unseen_D     in the model domain, never seen (balance prior + deficit)
#   observed_out seen, outside the domain (level of observed non-majors)
#   trend        neither (trend plus residual)
#
#   Rscript validation/sparse_joint/fig4d_check/decompose.R
#
# Writes results/deltaf_rows.rds and results/decompose.md.

base <- "validation/sparse_joint/fig4d_check/results"
input <- file.path(base, "input")
scores <- readRDS(file.path(input, "tables", "salehi_scores.Rds"))

select_fits <- function(design, method) {
  s <- scores[scores$design == design & scores$method == method, c("fi", "min_obs", "xv", "ntrain", "train_treat")]
  s <- s[!duplicated(paste(s$fi, s$min_obs)), ]
  s <- s[!is.na(s$xv) & s$xv > 0, ]
  s <- do.call(rbind, by(seq_len(nrow(s)), s$fi, function(i) {
    sub <- s[i, , drop = FALSE]
    sub <- sub[sub$xv == max(sub$xv), ]
    sub[order(as.numeric(sub$min_obs)), ][1, ]
  }))
  rownames(s) <- NULL
  s
}

neighbours <- function(k) {
  v <- as.integer(strsplit(k, ".", fixed = TRUE)[[1]])
  out <- character()
  for (i in seq_along(v)) for (d in c(-1L, 1L)) {
    w <- v
    w[i] <- w[i] + d
    out <- c(out, paste(w, collapse = "."))
  }
  out
}

arm_rows <- function(design, method, label) {
  sel <- select_fits(design, method)
  # Two-time-point trajectories are fitted only as ends.
  rows <- lapply(seq_len(nrow(sel)), function(i) {
    fi <- sel$fi[i]
    d <- if (design == "all" && grepl("_l_2_", fi)) "ends" else design
    m <- if (d == "ends" && grepl("^alfak2", method)) "alfak2" else method
    dir <- file.path(input, "fits_salehi", fi, d, m, paste0("minobs_", sel$min_obs[i]))
    lp <- file.path(dir, "landscape.Rds")
    if (!file.exists(lp)) return(NULL)
    ls <- readRDS(lp)
    lut <- stats::setNames(ls$mean, ls$k)
    fit <- if (file.exists(file.path(dir, "alfak2_fit.Rds"))) readRDS(file.path(dir, "alfak2_fit.Rds")) else NULL
    focal <- ls$k[ls$fq]
    out <- do.call(rbind, lapply(focal, function(fk) {
      nn <- neighbours(fk)
      data.frame(k_id = fk, nn = nn, f = unname(lut[fk]), f_nn = unname(lut[nn]), stringsAsFactors = FALSE)
    }))
    out$deltaf <- out$f - out$f_nn
    if (!is.null(fit)) {
      counts <- fit$data$counts
      seen <- rownames(counts)[rowSums(counts) > 0]
      if (is.null(seen)) seen <- fit$data$labels[rowSums(counts) > 0]
      in_d <- out$nn %in% fit$summary$karyotype
      obs <- out$nn %in% seen
      out$source <- ifelse(in_d & obs, "observed_D", ifelse(in_d, "unseen_D", ifelse(obs, "observed_out", "trend")))
      out$focal_major <- out$k_id %in% fit$summary$karyotype[fit$summary$major]
      out$nn_shell <- fit$summary$shell[match(out$nn, fit$summary$karyotype)]
      out$history <- fit$settings$history
      out$dt_total <- sum(fit$data$dt)
      out$n_time <- ncol(counts)
      out$cells_mean <- mean(colSums(counts))
      out$n_major <- sum(fit$summary$major)
    } else {
      out$source <- NA_character_
    }
    out$fi <- fi
    out$min_obs <- sel$min_obs[i]
    out$treated <- grepl("y", sel$train_treat[i], fixed = TRUE)
    out
  })
  x <- do.call(rbind, rows)
  x$arm <- label
  x$context <- ifelse(grepl("^SA906|^SA039", x$fi), "invitro", "PDX")
  x$line <- sub("^(SA[0-9]+).*", "\\1", x$fi)
  x
}

arms <- list(c("ends", "alfak2", "alfak2_2"), c("all", "alfak2", "alfak2_all"),
             c("ends", "alfakR", "alfakR_2"), c("all", "alfakR", "alfakR_all"))
x <- do.call(rbind, lapply(arms, function(a) {
  r <- arm_rows(a[1], a[2], a[3])
  for (col in setdiff(c("source", "focal_major", "nn_shell", "history", "dt_total", "n_time", "cells_mean", "n_major"),
                      names(r))) r[[col]] <- NA
  r
}))
x <- x[is.finite(x$deltaf), ]
saveRDS(x, file.path(base, "deltaf_rows.rds"))

md <- function(d) {
  d[] <- lapply(d, function(v) if (is.numeric(v)) signif(v, 3) else v)
  c(paste("|", paste(names(d), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"),
    apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
ctl <- x[!x$treated, ]
by_arm <- aggregate(cbind(deltaf, absdf = abs(deltaf)) ~ arm + context, ctl, median)
by_arm$n <- aggregate(deltaf ~ arm + context, ctl, length)$deltaf
a2 <- ctl[grepl("^alfak2", ctl$arm), ]
by_source <- aggregate(cbind(deltaf, absdf = abs(deltaf)) ~ arm + context + source, a2, median)
share <- aggregate(deltaf ~ arm + context + source, a2, length)
by_source$n <- share$deltaf
by_source$share <- ave(by_source$n, by_source$arm, by_source$context, FUN = function(v) v / sum(v))
by_line <- aggregate(cbind(deltaf, absdf = abs(deltaf), dt_total, cells_mean, n_time, history) ~ arm + context + line,
                     a2, median)
traj <- unique(a2[, c("arm", "context", "line", "fi", "dt_total", "cells_mean", "n_time", "history", "n_major")])
traj_sum <- aggregate(cbind(dt_total, cells_mean, n_time, history, n_major) ~ arm + context, traj, median)
traj_sum$fits <- aggregate(fi ~ arm + context, traj, length)$fi

# Paired with alfakR on the same trajectory, focal and neighbour.
pair <- merge(ctl[ctl$arm == "alfak2_2", c("fi", "k_id", "nn", "deltaf", "source", "context")],
              ctl[ctl$arm == "alfakR_2", c("fi", "k_id", "nn", "deltaf")], by = c("fi", "k_id", "nn"),
              suffixes = c("_alfak2", "_alfakR"))
paired <- aggregate(cbind(deltaf_alfak2, deltaf_alfakR) ~ context + source, pair, median)
paired$n <- aggregate(deltaf_alfak2 ~ context + source, pair, length)$deltaf_alfak2

writeLines(c("# Fig. 4d check, step 1: decomposition (control fits)", "",
             "## delta f by arm and context", "", md(by_arm), "",
             "## alfak2: by the source of the neighbour's value", "", md(by_source), "",
             "## alfak2: by line", "", md(by_line), "",
             "## alfak2: the fits (median per trajectory)", "", md(traj_sum), "",
             "## alfak2 vs alfakR, two time points, same trajectory, focal and neighbour", "", md(paired)),
           file.path(base, "decompose.md"))
cat(readLines(file.path(base, "decompose.md")), sep = "\n")
