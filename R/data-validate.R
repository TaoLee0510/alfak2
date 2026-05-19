validate_count_matrix <- function(counts) {
  if (is.data.frame(counts)) counts <- as.matrix(counts)
  if (!is.matrix(counts) || ncol(counts) != 2L) {
    stop("`counts` must be a matrix or data frame with exactly two columns.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
    stop("`counts` must have non-empty karyotype row names.", call. = FALSE)
  }
  storage.mode(counts) <- "double"
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("`counts` must contain finite non-negative values.", call. = FALSE)
  }
  if (any(abs(counts - round(counts)) > 1e-8)) {
    stop("`counts` must contain integer counts.", call. = FALSE)
  }
  counts <- matrix(as.integer(round(counts)), nrow = nrow(counts),
                   dimnames = dimnames(counts))
  keep <- rowSums(counts) > 0L
  if (!any(keep)) stop("At least one karyotype must have non-zero counts.", call. = FALSE)
  counts[keep, , drop = FALSE]
}

validate_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower || x > upper) {
    stop(sprintf("`%s` must be one finite number in [%s, %s].", name, lower, upper), call. = FALSE)
  }
  invisible(x)
}

validate_positive_grid <- function(x, name, lower = .Machine$double.eps) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < lower)) {
    stop(sprintf("`%s` must contain finite values >= %s.", name, lower), call. = FALSE)
  }
  as.numeric(x)
}

match_observation_model <- function(x) {
  match.arg(x, c("multinomial", "dirichlet_multinomial"))
}

match_observation_weight_mode <- function(x) {
  match.arg(x, c("likelihood", "fractional_count"))
}

alfak2_error_log_dir <- function() {
  env_dir <- Sys.getenv("ALFAK2_ERROR_LOG_DIR", unset = NA_character_)
  opt_dir <- getOption("alfak2.error_log_dir", default = NULL)
  out <- if (!is.na(env_dir) && nzchar(env_dir)) env_dir else opt_dir
  if (is.null(out)) out <- file.path(tempdir(), "alfak2_error_logs")
  out
}

alfak2_write_condition_log <- function(log_dir, prefix, message, diagnostics = list()) {
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%OS3")
  stamp <- gsub("[^0-9A-Za-z_.-]+", "_", stamp)
  path <- file.path(log_dir, paste0("alfak2_", prefix, "_", stamp, "_", Sys.getpid(), ".log"))
  diagnostic_lines <- utils::capture.output(str(diagnostics, give.attr = FALSE, vec.len = 50, max.level = 12))
  writeLines(
    c(
      paste0("timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%OS3 %Z")),
      paste0("pid: ", Sys.getpid()),
      paste0("type: ", prefix),
      paste0("message: ", message),
      "",
      "diagnostics:",
      diagnostic_lines
    ),
    con = path,
    useBytes = TRUE
  )
  path
}

alfak2_write_error_log <- function(message, diagnostics = list()) {
  alfak2_write_condition_log(alfak2_error_log_dir(), "error", message, diagnostics)
}

alfak2_warning_log_dir <- function() {
  env_dir <- Sys.getenv("ALFAK2_WARNING_LOG_DIR", unset = NA_character_)
  opt_dir <- getOption("alfak2.warning_log_dir", default = NULL)
  out <- if (!is.na(env_dir) && nzchar(env_dir)) env_dir else opt_dir
  if (is.null(out)) out <- alfak2_error_log_dir()
  out
}

alfak2_write_warning_log <- function(message, diagnostics = list()) {
  alfak2_write_condition_log(alfak2_warning_log_dir(), "warning", message, diagnostics)
}

alfak2_warn <- function(message, diagnostics = list()) {
  log_error <- NULL
  log_path <- tryCatch(
    alfak2_write_warning_log(message, diagnostics),
    error = function(e) {
      log_error <<- conditionMessage(e)
      NA_character_
    }
  )
  warning_message <- message
  if (!is.na(log_path)) {
    warning_message <- paste0(warning_message, " Diagnostics were written to: ", log_path)
  } else if (!is.null(log_error)) {
    warning_message <- paste0(warning_message, " Diagnostic log write failed: ", log_error)
  }
  warning(warning_message, call. = FALSE)
  invisible(log_path)
}

alfak2_abort <- function(message, diagnostics = list()) {
  log_error <- NULL
  log_path <- tryCatch(
    alfak2_write_error_log(message, diagnostics),
    error = function(e) {
      log_error <<- conditionMessage(e)
      NA_character_
    }
  )
  if (!is.na(log_path)) {
    message <- paste0(message, " Diagnostics were written to: ", log_path)
  } else if (!is.null(log_error)) {
    message <- paste0(message, " Diagnostic log write failed: ", log_error)
  }
  stop(structure(
    list(
      message = message,
      call = NULL,
      diagnostics = diagnostics,
      log_path = log_path,
      log_error = log_error
    ),
    class = c("alfak2_error", "error", "condition")
  ))
}
