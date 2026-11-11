# Do the files of alfak2::alfak_compat() have the structure of alfakR::alfak()'s?
#
#   Rscript validation/sparse_joint/alfakR_compat.R DATASET_ID [OUT_DIR]
#
# Runs both on one benchmark dataset with the same arguments and compares
# every saved object: classes, names, columns, dimensions and attributes,
# then the karyotypes of the landscape and the frequent and neighbour sets.
# Set STT_R_LIB to a library with the alfak2 build under test and
# ALFAKR_LIB to one with alfakR (branch dev).

for (variable in c("STT_R_LIB", "ALFAKR_LIB")) {
  library_dir <- Sys.getenv(variable, unset = "")
  if (nzchar(library_dir)) .libPaths(c(library_dir, .libPaths()))
}
args <- commandArgs(TRUE)
if (!length(args)) stop("Give a dataset id, e.g. dev_L10_P01_unbalanced_3000_100_360d.", call. = FALSE)
id <- args[[1L]]
out_dir <- if (length(args) > 1L) args[[2L]] else file.path(tempdir(), "alfakR_compat")
dataset <- readRDS(file.path("benchmark/results/sparse_two_timepoint/datasets", paste0(id, ".rds")))
counts <- dataset$counts
colnames(counts) <- c("0", "1")
yi <- list(x = counts, dt = dataset$interval)
paths <- file.path(out_dir, id, c("alfak2", "alfakR"))
unlink(paths, recursive = TRUE)
cat(sprintf("%s: %d karyotypes, %s cells, %g days; alfakR %s\n", id, nrow(counts),
            paste(colSums(counts), collapse = "/"), dataset$interval, packageVersion("alfakR")))

set.seed(1)
started <- proc.time()[["elapsed"]]
r2r <- alfak2::alfak_compat(yi, paths[[1L]], minobs = 20, nboot = 45, n0 = 5e4, nb = 2e6, pm = dataset$beta,
                            birth_rate = 1, min_cn = 1, max_cn = 8)
cat(sprintf("alfak2::alfak_compat %.0f s, R2R %.3f\n", proc.time()[["elapsed"]] - started, r2r))
set.seed(1)
started <- proc.time()[["elapsed"]]
r2r <- suppressWarnings(alfakR::alfak(yi, paths[[2L]], minobs = 20, nboot = 45, n0 = 5e4, nb = 2e6, pm = dataset$beta))
cat(sprintf("alfakR::alfak %.0f s, R2R %.3f\n", proc.time()[["elapsed"]] - started, r2r))

# One line per object and nested element; sizes are replaced by "#" so that
# only structure is compared.
describe <- function(x, path) {
  line <- sprintf("%s class=%s", path, paste(class(x), collapse = "/"))
  extra <- setdiff(names(attributes(x)), c("dim", "dimnames", "names", "class", "row.names"))
  if (length(extra)) line <- paste0(line, " attributes=", paste(sort(extra), collapse = ","))
  if (is.matrix(x)) line <- paste0(line, " rownames=", !is.null(rownames(x)), " colnames=", !is.null(colnames(x)))
  if (is.data.frame(x)) {
    line <- paste0(line, " columns=", paste(names(x), vapply(x, function(v) class(v)[1L], ""), sep = ":", collapse = ","))
  }
  lines <- line
  if (is.list(x) && !is.data.frame(x)) for (name in names(x)) lines <- c(lines, describe(x[[name]], paste0(path, "$", name)))
  lines
}
files <- c("bootstrap_res.Rds", "landscape.Rds", "landscape_posterior_samples.Rds", "xval.Rds")
for (file in files) {
  a <- describe(readRDS(file.path(paths[[1L]], file)), file)
  b <- describe(readRDS(file.path(paths[[2L]], file)), file)
  cat(sprintf("\n%s: %s\n", file, if (identical(a, b)) "same structure" else "differs"))
  if (!identical(a, b)) {
    cat("  alfak2 only:", setdiff(a, b), sep = "\n    ")
    cat("  alfakR only:", setdiff(b, a), sep = "\n    ")
  }
}
landscape <- lapply(paths, function(p) readRDS(file.path(p, "landscape.Rds")))
boot <- lapply(paths, function(p) readRDS(file.path(p, "bootstrap_res.Rds")))
cat(sprintf("\nlandscape rows %d / %d, same karyotypes in the same order: %s\n", nrow(landscape[[1L]]),
            nrow(landscape[[2L]]), identical(landscape[[1L]]$k, landscape[[2L]]$k)))
cat(sprintf("frequent karyotypes identical: %s; neighbours identical: %s\n",
            identical(colnames(boot[[1L]]$final_fitness), colnames(boot[[2L]]$final_fitness)),
            identical(colnames(boot[[1L]]$nn_fitness), colnames(boot[[2L]]$nn_fitness))))
