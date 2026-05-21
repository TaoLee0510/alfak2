#!/usr/bin/env Rscript

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "benchmark/run_full_second_layer_9method_balanced_comparison.R"
script_path <- normalizePath(file.path(dirname(this_file), "run_full_second_layer_comparison.R"), mustWork = TRUE)
source(script_path, chdir = FALSE)
