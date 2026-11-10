#!/usr/bin/env Rscript
# Default invocation prepares a resource/configuration plan only. No fitting.
args <- commandArgs(TRUE)
get_arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}
output <- get_arg("output", tempfile("alfak2_scalability_", fileext = ".csv"))
plan <- expand.grid(n_chr = c(4L, 8L, 12L, 22L), max_nodes = c(128L, 512L, 1000L, 5000L))
plan$min_cn <- 1L; plan$max_cn <- 3L; plan$shell_depth <- 2L
plan$dense_fitness_covariance_mib <- 8 * as.double(plan$max_nodes)^2 / 1024^2
plan$two_parameter_block_hessian_lower_bound_mib <- 8 * (2 * as.double(plan$max_nodes))^2 / 1024^2
plan$status <- "planned_not_run"
if ("--run-build" %in% args) {
  library_dir <- get_arg("library")
  if (is.null(library_dir)) stop("--run-build requires --library")
  library(alfak2, lib.loc = library_dir)
  plan$actual_nodes <- plan$graph_seconds <- NA_real_
  plan$truncated <- NA
  for (i in seq_len(nrow(plan))) {
    p <- plan$n_chr[i]
    k <- matrix(2L, nrow = p + 1L, ncol = p)
    for (j in seq_len(p)) k[j + 1L, j] <- 3L
    counts <- cbind(t0 = c(150L, rep(2L, p)), t1 = c(140L, rep(3L, p)))
    rownames(counts) <- apply(k, 1L, paste, collapse = ".")
    started <- proc.time()[["elapsed"]]
    graph <- tryCatch(build_karyotype_graph(prepare_alfak2_data(counts, dt = 1, beta = 0.005), beta = 0.005, shell_depth = plan$shell_depth[i],
                                           min_cn = plan$min_cn[i], max_cn = plan$max_cn[i],
                                           max_nodes = plan$max_nodes[i]), error = function(e) e)
    plan$graph_seconds[i] <- proc.time()[["elapsed"]] - started
    if (inherits(graph, "error")) plan$status[i] <- conditionMessage(graph) else {
      plan$actual_nodes[i] <- length(graph$labels)
      plan$truncated[i] <- length(graph$labels) >= plan$max_nodes[i]
      plan$status[i] <- "graph_only_completed_no_fit"
    }
  }
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(plan, output, row.names = FALSE)
print(plan)
cat("Saved", normalizePath(output), "\n")
