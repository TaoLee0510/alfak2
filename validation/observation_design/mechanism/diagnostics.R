reference_load <- function(repo, compile = TRUE) {
  native <- Sys.getenv("ALFAK2_REFERENCE_DLL",
    "/private/tmp/alfak2_iteration2_20261104_0518/current-lib/alfak2/libs/alfak2.so")
  if (!file.exists(native)) stop("The recorded oracle-capable DLL is required")
  requireNamespace("Rcpp"); requireNamespace("TMB"); requireNamespace("posterior")
  pkg <- new.env(parent = globalenv())
  binary <- dyn.load(native)
  for (nm in names(getDLLRegisteredRoutines(binary)$`.Call`))
    assign(nm, getDLLRegisteredRoutines(binary)$`.Call`[[nm]], pkg)
  for (f in sort(list.files(file.path(repo, "R"), "[.]R$", full.names = TRUE))) sys.source(f, pkg)
  independent <- new.env(parent = globalenv())
  sys.source(file.path(repo, "validation/accuracy_round2/process_prepare_data.R"), independent)
  cpp <- new.env(parent = globalenv())
  if (compile) Rcpp::sourceCpp(file.path(repo, "validation/posterior_reference_4to8/reference.cpp"),
    env = cpp, cacheDir = "/private/tmp/alfak2_reference_4to8_cpp", showOutput = FALSE)
  list(pkg = pkg, independent = independent, cpp = cpp, native = native)
}

reference_case <- function(n, profile, beta, growth, runtime) {
  p <- as.integer(log2(n)); stopifnot(n %in% c(4L, 8L), profile %in% c("balanced", "sparse"))
  k <- as.matrix(expand.grid(rep(list(1:2), p)))
  labels <- apply(k, 1L, paste, collapse = ".")
  graph <- runtime$pkg$build_karyotype_graph(labels, beta = beta,
    shell_depth = 0, min_cn = 1, max_cn = 2, max_nodes = n)
  prior <- runtime$pkg$build_joint_landscape_prior(graph, forward_model = "division_process",
    baseline_sd = .5, main_sd = .07, curvature_sd = .04, interaction_sd = .04,
    residual_sd = .125, edge_sd = .25, interaction_support = "isotropic")
  covariance <- solve(as.matrix(prior$precision)) + tcrossprod(sweep(prior$design, 2L, prior$coef_sd, "*"))
  y <- if (n == 4L && profile == "balanced") cbind(c(38, 36, 39, 37), c(26, 33, 42, 49)) else
    if (n == 4L) cbind(c(126, 24, 0, 0), c(112, 24, 14, 0)) else
    if (profile == "balanced") cbind(c(20, 18, 19, 17, 21, 16, 20, 19), c(12, 15, 17, 19, 21, 18, 24, 24)) else
    cbind(c(110, 25, 10, 5, 0, 0, 0, 0), c(94, 24, 11, 8, 9, 4, 0, 0))
  rownames(y) <- labels; colnames(y) <- c("t0", "t1")
  y <- y[match(graph$labels, labels), , drop = FALSE]
  basis <- stats::contr.helmert(n)
  basis <- sweep(basis, 2L, sqrt(colSums(basis^2)), "/")
  probability <- runtime$independent$process_independent_operator(graph$karyotypes, rep(1, n), beta)$probability
  spec <- list(n = n, q = 2L * n - 1L, steps = 15L, dt = 1.5,
    offset = prior$offset, lower = t(chol(covariance)), eta_basis = 2 * basis,
    known_pi0 = numeric(), probability = probability, beta_zero = beta == 0,
    y0 = as.numeric(y[, 1]), y1 = as.numeric(y[, 2]),
    growth = 1.5, growth_sd = if (growth) .1 else 0)
  distance <- as.matrix(stats::dist(graph$karyotypes, method = "manhattan"))
  id <- sprintf("n%d_%s_beta%s_%s", n, profile, if (beta == 0) "0" else "0p005", if (growth) "growth" else "counts")
  list(id = id, n = n, profile = profile, beta = beta, growth = growth,
    counts = y, graph = graph, prior = prior, covariance = covariance,
    specification = spec, weight = y[, 1] / sum(y[, 1]),
    edges = which(distance == 1 & upper.tri(distance), arr.ind = TRUE))
}

reference_fit <- function(case, runtime, collapsed = FALSE) {
  prior <- case$prior
  if (collapsed) {
    prior$design <- matrix(numeric(), case$n, 0L); prior$coef_sd <- numeric()
    prior$precision <- methods::as(Matrix::Matrix(solve(case$covariance), sparse = TRUE), "generalMatrix")
    prior$logdet <- -as.numeric(determinant(case$covariance, logarithm = TRUE)$modulus)
  }
  data <- runtime$pkg$prepare_alfak2_data(case$counts, dt = 1.5, beta = case$beta)
  runtime$pkg$fit_local_posterior(data, graph = case$graph, structured_prior = prior,
    forward_model = "division_process", forward_step = .1, eta_prior_sd = 2,
    eta_borrowed_prior_mean = 0, eta_borrowed_prior_sd = 2, eta_distance_penalty = 0,
    population_sizes = if (case$growth) c(1000, 1000 * exp(1.5)) else NULL,
    population_log_growth_sd = .1, sdreport_mode = "fitness_only",
    conditional_warmstart = "none", return_tmb_objects = TRUE)
}

reference_observables <- function(z, case) {
  n <- case$n; sp <- case$specification
  g <- sweep(z[, seq_len(n), drop = FALSE] %*% t(sp$lower), 2L, sp$offset, "+")
  birth <- exp(g); relative <- birth - as.numeric(birth %*% case$weight)
  centered <- birth - rowMeans(birth)
  eta <- if (ncol(z) > n) z[, n + seq_len(ncol(z) - n), drop = FALSE] * 2 else matrix(numeric(), nrow(z), 0L)
  edge <- birth[, case$edges[, 2], drop = FALSE] - birth[, case$edges[, 1], drop = FALSE]
  out <- cbind(g, birth, relative, centered, eta, edge)
  colnames(out) <- c(paste0("log_birth_", 1:n), paste0("birth_", 1:n),
    paste0("relative_", 1:n), paste0("centered_", 1:n),
    if (ncol(eta)) paste0("eta_contrast_", seq_len(ncol(eta))) else character(), paste0("edge_", 1:ncol(edge)))
  out
}

reference_diagnostics <- function(chains, case) {
  obs <- lapply(chains, function(x) reference_observables(x$draws, case))
  vnames <- colnames(obs[[1L]])
  arr <- array(NA_real_, c(nrow(obs[[1]]), length(obs), length(vnames)),
    dimnames = list(NULL, NULL, vnames))
  for (j in seq_along(obs)) arr[, j, ] <- obs[[j]]
  summarize <- function(x) data.frame(variable = dimnames(x)[[3]],
    mean = apply(x, 3, mean), sd = apply(x, 3, stats::sd),
    rhat = apply(x, 3, posterior::rhat), ess_bulk = apply(x, 3, posterior::ess_bulk),
    ess_tail = apply(x, 3, posterior::ess_tail), mcse_mean = apply(x, 3, posterior::mcse_mean),
    mcse_sd = apply(x, 3, posterior::mcse_sd), row.names = NULL)
  summary <- summarize(arr)
  covariance <- covariance_mcse <- product_diag <- list()
  for (type in c("log_birth", "birth", "relative", "centered", "edge", if(any(startsWith(vnames, "eta_contrast_"))) "eta_contrast" else character())) {
    cols <- which(startsWith(vnames, paste0(type, "_")))
    values <- lapply(obs, function(x) x[, cols, drop = FALSE])
    joined <- do.call(rbind, values); mu <- colMeans(joined)
    covariance[[type]] <- stats::cov(joined)
    m <- length(cols); pair <- which(upper.tri(matrix(0, m, m), diag = TRUE), arr.ind = TRUE)
    products <- array(NA_real_, c(nrow(values[[1]]), length(values), nrow(pair)),
      dimnames = list(NULL, NULL, paste0(type, "_centered_product_", pair[, 1], "_", pair[, 2])))
    for (ch in seq_along(values)) {
      v <- sweep(values[[ch]], 2L, mu, "-")
      products[, ch, ] <- v[, pair[, 1], drop = FALSE] * v[, pair[, 2], drop = FALSE]
    }
    d <- summarize(products); product_diag[[type]] <- d
    mcse <- matrix(0, m, m); mcse[pair] <- d$mcse_mean
    mcse <- mcse + t(mcse) - diag(diag(mcse), m)
    covariance_mcse[[type]] <- mcse
  }
  product_summary <- do.call(rbind, product_diag)
  all_diag <- rbind(summary, product_summary)
  varying <- is.finite(all_diag$sd) & all_diag$sd > 1e-12
  precision <- summary$sd > 1e-12
  cov_relative_mcse <- vapply(names(covariance), function(type)
    norm(covariance_mcse[[type]], "F") / norm(covariance[[type]], "F"), numeric(1))
  ok <- all(is.finite(unlist(all_diag[varying, c("rhat", "ess_bulk", "ess_tail")]))) &&
    max(all_diag$rhat[varying]) <= 1.01 && min(all_diag$ess_bulk[varying]) >= 1000 &&
    min(all_diag$ess_tail[varying]) >= 1000 &&
    max(summary$mcse_mean[precision] / summary$sd[precision]) <= .025 &&
    max(summary$mcse_sd[precision] / summary$sd[precision]) <= .025 &&
    max(cov_relative_mcse) <= .05
  list(summary = summary, product_summary = product_summary, covariance = covariance,
    covariance_mcse = covariance_mcse, covariance_relative_mcse = cov_relative_mcse,
    convergence = data.frame(case_id = case$id, passed = ok, retained_per_chain = nrow(obs[[1]]),
      max_rhat = max(all_diag$rhat[varying]), min_bulk_ess = min(all_diag$ess_bulk[varying]),
      min_tail_ess = min(all_diag$ess_tail[varying]),
      max_mean_mcse_sd = max(summary$mcse_mean[precision] / summary$sd[precision]),
      max_sd_mcse_sd = max(summary$mcse_sd[precision] / summary$sd[precision]),
      max_covariance_relative_mcse = max(cov_relative_mcse)),
    quantiles = t(apply(do.call(rbind, obs), 2L, stats::quantile, c(.025, .5, .975))))
}

reference_chains <- function(case, runtime, retained, warmup = 5000L, seed = 720000L,
                              previous = NULL, cores = 2L) {
  work <- function(ch) {
    set.seed(seed + ch)
    proposal <- case$proposal
    if (is.null(proposal)) {
      q <- case$specification$q
      proposal <- list(mean=rep(.2,q),lower=diag(sqrt(2),q),precision=diag(.5,q))
    }
    initial <- if (is.null(previous)) as.numeric(proposal$mean + proposal$lower %*% stats::rnorm(case$specification$q,sd=1.5)) else previous[[ch]]$last
    set.seed(seed + 100L + ch)
    result <- runtime$cpp$ref_preconditioned_ess_chain_cpp(case$specification, initial,
      if (is.null(previous)) warmup else 0L, retained, seed + 100L + ch, proposal$mean,proposal$lower,proposal$precision)
    if (!is.null(previous)) {
      result$draws <- rbind(previous[[ch]]$draws, result$draws)
      for (nm in c("likelihood_evaluations", "rejected_angles", "nonfinite_proposals"))
        result[[nm]] <- previous[[ch]][[nm]] + result[[nm]]
    }
    result
  }
  out <- parallel::mclapply(1:4, work, mc.cores = cores, mc.set.seed = FALSE)
  if (any(vapply(out, inherits, logical(1), "try-error"))) stop("Reference chain process failed")
  out
}
