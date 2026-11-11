test_that("neighbours follow alfakR's order and rules", {
  # Chromosome by chromosome, one copy fewer then one more; no nullisomy.
  expect_equal(alfakR_neighbours("2.1"), c("1.1", "3.1", "2.2"))
  expect_equal(alfakR_neighbours(c("2.2", "2.3")), c("1.2", "3.2", "2.1", "1.3", "3.3", "2.4"))
  expect_equal(alfakR_neighbours(character()), character())
})

test_that("a rebuilt posterior reproduces the fit", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  mode <- sparse_joint_rebuild(fit, posterior = FALSE)
  expect_equal(mode$g, fit$posterior$g, tolerance = 1e-6)
  expect_equal(mode$zeta, fit$posterior$zeta, tolerance = 1e-6)
  solver <- sparse_joint_rebuild(fit)$solver
  expect_equal(sqrt(solver$variance(which(solver$names == "g"))), fit$summary$log_birth_sd, tolerance = 1e-6)
  set.seed(4)
  draws <- solver$draw(3000)
  g <- which(solver$names == "g")[1:3]
  expect_equal(apply(draws[g, ], 1L, stats::sd), fit$summary$log_birth_sd[1:3], tolerance = 0.06)
})

test_that("as_alfakR returns alfakR's objects for the model's posterior", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  set.seed(3)
  out <- as_alfakR(fit, minobs = 20, nboot = 7)
  expect_s3_class(out, "alfak2_alfakR_output")
  fq <- rownames(fit$data$counts)[rowSums(fit$data$counts) >= 20]
  nn <- alfakR_neighbours(fq)

  boot <- out$bootstrap_res
  expect_named(boot, c("initial_fitness", "final_fitness", "initial_frequencies", "final_frequencies", "nn_fitness"))
  expect_equal(dim(boot$final_fitness), c(7L, length(fq)))
  expect_equal(colnames(boot$final_fitness), fq)
  expect_null(rownames(boot$final_fitness))
  expect_equal(colnames(boot$nn_fitness), nn)
  expect_equal(boot$initial_fitness, boot$final_fitness)
  expect_equal(unname(rowSums(boot$final_frequencies)), rep(1, 7))
  expect_named(attr(boot, "bootstrap_retry_log"), c("target_replicate", "attempt", "error_message", "warning_message"))
  expect_equal(attr(boot, "bootstrap_retry_summary")$successful_replicates, 7L)

  land <- out$landscape
  expect_named(land, c("k", "mean", "median", "sd", "fq", "nn"))
  expect_equal(land$k, unique(c(fq, nn, alfakR_neighbours(nn))))
  expect_equal(land$fq, land$k %in% fq)
  expect_equal(land$nn, land$k %in% nn)
  post <- out$landscape_posterior_samples
  expect_equal(dim(post), c(nrow(land), 7L))
  expect_null(dimnames(post))
  # The landscape summarises the draws, as in alfakR, and is the model's
  # prediction on alfakR's scale.
  expect_equal(land$mean, rowMeans(post))
  expect_equal(land$median, apply(post, 1L, stats::median))
  expect_equal(land$sd, apply(post, 1L, stats::sd))
  prediction <- predict(fit, land$k)
  expect_equal(land$mean, prediction$relative_fitness + log(1e7 / 1e5) / sim$dt, tolerance = 1e-12)
  expect_equal(land$sd, prediction$relative_sd, tolerance = 1e-12)
  expect_equal(unname(t(boot$final_fitness)), post[match(fq, land$k), , drop = FALSE])
  expect_equal(unname(t(boot$nn_fitness)), post[match(nn, land$k), , drop = FALSE])

  xval <- out$xval
  expect_named(xval, c("R2R", "predictions", "observations", "xval_data", "tmp"))
  expect_length(xval$R2R, 1L)
  expect_equal(colnames(xval$tmp), c("test_f", "est_f"))
  expect_gt(nrow(xval$tmp), 1L)
  expect_named(xval$xval_data, c("k", "observation", "prediction"))
  expect_equal(xval$xval_data$k, rownames(xval$tmp))
  expect_equal(unname(xval$observations), land$mean[match(rownames(xval$tmp), land$k)])
  expect_true(all(rownames(xval$tmp) %in% fit$summary$karyotype))
  centred <- function(v) v - mean(v)
  expect_equal(xval$R2R, 1 - sum((centred(xval$predictions) - centred(xval$observations))^2) /
                 sum(centred(xval$observations)^2))
})

test_that("as_alfakR keeps the model's level or reports birth rates on request", {
  sim <- simulate_sparse_joint_counts()
  fit <- fit_sparse_joint(sim$counts, dt = sim$dt, beta = sim$beta, max_cn = 4)
  own <- as_alfakR(fit, minobs = 20, nboot = 3, n0 = NULL, nb = NULL, xval = FALSE)
  prediction <- predict(fit, own$landscape$k)
  expect_equal(own$landscape$mean, prediction$relative_fitness + fit$trend$reference$fitness, tolerance = 1e-12)
  expect_true(is.na(own$xval$R2R))
  expect_equal(nrow(own$xval$tmp), 0L)
  births <- as_alfakR(fit, minobs = 20, nboot = 3, correct_efflux = TRUE, xval = FALSE)
  self <- 2 * sparse_joint_self_probability(parse_karyotypes(births$landscape$k), sim$beta) - 1
  expect_equal(births$landscape$mean, (prediction$relative_fitness + log(100) / sim$dt) / self, tolerance = 1e-12)
  single <- as_alfakR(fit, minobs = 20, nboot = 1, xval = FALSE)
  expect_true(all(is.na(single$landscape$sd)))
  expect_error(as_alfakR(fit, minobs = 1e6), "no frequent karyotypes")
  expect_error(as_alfakR(fit, n0 = 1e5, nb = NULL), "both")
})

test_that("alfak_compat() takes alfakR's input and writes its files", {
  sim <- simulate_sparse_joint_counts()
  x <- sim$counts
  colnames(x) <- c("0", "1")
  dir <- tempfile("alfak2_alfak_")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  set.seed(5)
  expect_invisible(r2r <- alfak_compat(list(x = x, dt = sim$dt), dir, minobs = 20, nboot = 5, pm = sim$beta,
                                landscape_data_output = TRUE, max_cn = 4))
  expect_setequal(list.files(dir), c("bootstrap_res.Rds", "landscape.Rds", "landscape_posterior_samples.Rds",
                                     "xval.Rds", "landscape_data.Rds", "alfak2_fit.Rds"))
  expect_identical(r2r, readRDS(file.path(dir, "xval.Rds"))$R2R)
  landscape <- readRDS(file.path(dir, "landscape.Rds"))
  models <- readRDS(file.path(dir, "landscape_data.Rds"))
  expect_named(models, c("mean", "median"))
  expect_equal(predict(models$mean, parse_karyotypes(landscape$k[1:4]))[, 1L], landscape$mean[1:4], tolerance = 1e-12)

  # The time between the samples comes from the column names or passage times.
  again <- tempfile("alfak2_alfak_")
  on.exit(unlink(again, recursive = TRUE), add = TRUE)
  bare <- x
  colnames(bare) <- NULL
  set.seed(5)
  alfak_compat(list(x = bare), again, passage_times = c(10, 10 + sim$dt), minobs = 20, nboot = 5, pm = sim$beta,
        max_cn = 4)
  expect_equal(readRDS(file.path(again, "landscape.Rds")), landscape)

  expect_error(alfak_compat(list(x = cbind(x, x[, 1L]), dt = 1), dir), "two samples")
  expect_error(alfak_compat(list(x = x, dt = sim$dt), dir, nn_prior = "none"), "nn_prior")
  expect_error(alfak_compat(list(x = x + 0.5, dt = sim$dt), dir), "non-integer")
})
