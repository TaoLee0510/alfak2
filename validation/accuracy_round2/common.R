round2_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[[1L]])))
}
round2_arg <- function(name, default = NULL, args = commandArgs(TRUE)) {
  i <- match(paste0("--", name), args)
  if (is.na(i)) default else args[[i + 1L]]
}
round2_sha256 <- function(files) {
  files <- normalizePath(files, mustWork = TRUE)
  command <- Sys.which("shasum")
  if (nzchar(command)) lines <- system2(command, c("-a", "256", shQuote(files)), stdout = TRUE) else {
    command <- Sys.which("sha256sum")
    if (!nzchar(command)) stop("SHA256 utility unavailable")
    lines <- system2(command, shQuote(files), stdout = TRUE)
  }
  hashes <- substring(lines, 1L, 64L)
  if (length(hashes) != length(files) || any(!grepl("^[0-9a-f]{64}$", hashes))) stop("SHA256 calculation failed")
  stats::setNames(hashes, files)
}
round2_hash_object <- function(x) {
  file <- tempfile("alfak2_round2_hash_")
  on.exit(unlink(file))
  saveRDS(x, file, version = 2)
  unname(round2_sha256(file))
}
round2_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE)
}
round2_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
round2_quantile <- function(x, probability) {
  x <- x[!is.na(x)]
  if (length(x)) unname(stats::quantile(x, probability, names = FALSE)) else NA_real_
}
round2_source_fingerprint <- function(source, dll = NULL) {
  source <- normalizePath(source)
  files <- c(list.files(file.path(source, "R"), "[.]R$", full.names = TRUE),
             list.files(file.path(source, "src"), "[.](cpp|h)$|Makevars", full.names = TRUE),
             file.path(source, "DESCRIPTION"))
  files <- sort(files[file.exists(files)])
  manifest <- data.frame(file = substring(files, nchar(source) + 2L), sha256 = unname(round2_sha256(files)))
  if (!is.null(dll)) manifest <- rbind(manifest, data.frame(file = "binary/alfak2.so", sha256 = unname(round2_sha256(dll))))
  list(manifest = manifest, sha256 = round2_hash_object(manifest))
}
