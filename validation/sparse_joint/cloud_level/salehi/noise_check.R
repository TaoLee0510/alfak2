#!/usr/bin/env Rscript
# Follow-up step 2, read only: how much of the observed change in the paper's
# Fig. 2d is sampling noise? For each row (trajectory, minobs, daughter
# passage) of KEYS, the parent passage's cells are split at random into two
# halves 20 times; the Wasserstein distance between the halves is pure
# sampling noise. Computed on all cells and on the karyotypes of the alfakR
# landscape of the "all" design (as the Fig. 2d axis is), with the cell
# counts of the parent and the daughter. Run from the paper benchmark's
# folder (sources R/), with BENCH_RESULTS the round's results folder.
#   Rscript noise_check.R KEYS.tsv OUT.csv
args <- commandArgs(TRUE)
for (f in list.files("R", full.names = TRUE)) source(f)
paths <- bench_paths()
keys <- utils::read.delim(args[[1]], stringsAsFactors = FALSE, colClasses = "character")
meta <- utils::read.csv(paths$salehi_meta)
inputs <- list.files(paths$salehi_inputs)
s2m <- function(k) do.call(rbind, lapply(strsplit(k, ".", fixed = TRUE), as.numeric))
wd <- function(a, b) {
  a <- a[a > 0]; b <- b[b > 0]
  if (!length(a) || !length(b)) return(NA_real_)
  transport::wasserstein(transport::wpp(s2m(names(a)), mass = a / sum(a)), transport::wpp(s2m(names(b)), mass = b / sum(b)))
}
halves <- function(p, keep = NULL, reps = 20L) {
  if (!is.null(keep)) p <- p[names(p) %in% keep]
  if (sum(p) < 4) return(NA_real_)
  cells <- rep(names(p), p)
  mean(vapply(seq_len(reps), function(r) {
    set.seed(r)
    i <- sample.int(length(cells), length(cells) %/% 2L)
    wd(table(cells[i]), table(cells[-i]))
  }, numeric(1)))
}
rows <- lapply(seq_len(nrow(keys)), function(j) {
  k <- keys[j, ]
  fj <- paste(meta[meta$uid == as.integer(k$dj), c("datasetname", "timepoint")], collapse = "_")
  tru <- readRDS(file.path(paths$salehi_inputs, inputs[grepl(fj, inputs)][1]))$x
  p <- stats::setNames(tru[, ncol(tru) - 1L], rownames(tru)); p <- p[p > 0]
  d <- stats::setNames(tru[, ncol(tru)], rownames(tru)); d <- d[d > 0]
  lscape <- readRDS(file.path(paths$results, "fits", "salehi", k$fi, Sys.getenv("NOISE_DESIGN", "all"), "alfakR", paste0("minobs_", k$min_obs), "landscape.Rds"))
  data.frame(k, cells_parent = sum(p), cells_daughter = sum(d), noise_all = halves(p), noise = halves(p, lscape$k),
             noise_daughter = halves(d, lscape$k), stringsAsFactors = FALSE)
})
out <- do.call(rbind, rows)
utils::write.csv(out, args[[2]], row.names = FALSE)
cat("rows", nrow(out), "\n")
