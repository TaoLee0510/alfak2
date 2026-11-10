#!/usr/bin/env Rscript
repo <- normalizePath(".")
directory <- file.path(repo, "validation/baseline_quadrature")
output <- file.path(directory, "results")
reference <- file.path(repo, "validation/posterior_reference_4to8/results")
source(file.path(repo, "validation/accuracy_round2/common.R"))
cases <- readRDS(file.path(reference, "frozen_cases.rds"))
status <- read.csv(file.path(output, "integration_status.csv"))
if (!setequal(status$case_id, names(cases))) stop("Incomplete candidate case set")
rows <- covariance <- scalar <- list()
for (case in cases) {
  q <- readRDS(file.path(output, paste0(case$id, "_integration.rds")))
  r <- readRDS(file.path(reference, paste0(case$id, "_summary.rds")))
  if (!isTRUE(r$reference$convergence$passed)) stop("Uncertified reference")
  if (!identical(q$status, "passed_quadrature_refinement")) {
    rows[[case$id]] <- data.frame(case_id = case$id, status = q$status,
      type = "unavailable", index = NA_integer_, reference_mean = NA_real_,
      reference_sd = NA_real_, candidate_mean = NA_real_, candidate_sd = NA_real_,
      standardized_mean_difference = NA_real_, sd_ratio = NA_real_,
      material_mean = NA, material_sd = NA)
    next
  }
  n <- case$n
  E <- matrix(0, nrow(case$edges), n)
  E[cbind(1:nrow(E), case$edges[, 1])] <- -1
  E[cbind(1:nrow(E), case$edges[, 2])] <- 1
  transforms <- list(log_birth = diag(n), birth = diag(n),
    relative = diag(n) - matrix(case$weight, n, n, byrow = TRUE),
    centered = diag(n) - matrix(1/n, n, n), edge = E)
  for (type in names(transforms)) {
    A <- transforms[[type]]
    m <- if (type == "log_birth") q$latent_mean else as.numeric(A %*% q$mean)
    C <- if (type == "log_birth") q$latent_covariance else A %*% q$covariance %*% t(A)
    ref <- r$reference$summary[match(paste0(type, "_", 1:nrow(A)), r$reference$summary$variable), ]
    sd <- sqrt(pmax(0, diag(C)))
    rows[[paste(case$id, type)]] <- data.frame(case_id = case$id, status = q$status,
      type = type, index = 1:nrow(A), reference_mean = ref$mean, reference_sd = ref$sd,
      candidate_mean = as.numeric(m), candidate_sd = sd,
      standardized_mean_difference = (m - ref$mean) / ref$sd, sd_ratio = sd/ref$sd,
      material_mean = abs(m-ref$mean) > .05*ref$sd & abs(m-ref$mean) > 4*ref$mcse_mean,
      material_sd = abs(sd-ref$sd) > .1*ref$sd & abs(sd-ref$sd) > 4*ref$mcse_sd)
    covariance[[paste(case$id, type)]] <- data.frame(case_id = case$id, type = type,
      relative_frobenius_error = norm(C-r$reference$covariance[[type]], "F") /
        norm(r$reference$covariance[[type]], "F"),
      reference_covariance_relative_mcse = r$reference$covariance_relative_mcse[[type]])
  }
  sd0 <- q$proposal$prior_sd; v <- as.numeric(solve(case$covariance, rep(sd0^2,n)))
  ref_log <- r$reference$summary[match(paste0("log_birth_", 1:n), r$reference$summary$variable), ]
  alpha_mean <- sum(v * (ref_log$mean-case$prior$offset))
  alpha_var <- q$proposal$conditional_prior_variance +
    as.numeric(crossprod(v, r$reference$covariance$log_birth %*% v))
  alpha <- vapply(q$components, `[[`, numeric(1), "intercept")
  qm <- sum(q$weights*alpha); qv <- sum(q$weights*(alpha-qm)^2)
  scalar[[case$id]] <- data.frame(case_id = case$id,
    reference_intercept_mean = alpha_mean, reference_intercept_sd = sqrt(alpha_var),
    original_laplace_intercept_mean = q$proposal$mean, original_laplace_intercept_sd = q$proposal$sd,
    candidate_intercept_mean = qm, candidate_intercept_sd = sqrt(qv),
    interpretation = "intercept moments recovered by the exact Gaussian conditional law, not observed population growth")
}
moments <- do.call(rbind, rows)
write.csv(moments, file.path(output, "moment_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, covariance), file.path(output, "covariance_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, scalar), file.path(output, "intercept_moments.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(names(cases), function(id) {
  x <- moments[moments$case_id==id & moments$type=="relative", ]
  lp <- readRDS(file.path(reference,paste0(id,"_summary.rds")))
  c <- cases[[id]]; A <- diag(c$n)-matrix(c$weight,c$n,c$n,byrow=TRUE)
  s <- lp$reference$summary[match(paste0("relative_",1:c$n),lp$reference$summary$variable),]
  old_ratio <- sqrt(diag(A%*%lp$laplace$fitness_covariance%*%t(A)))/s$sd
  data.frame(case_id=id,n=c$n,profile=c$profile,beta=c$beta,growth=c$growth,
    status=status$status[match(id,status$case_id)],
    original_sd_ratio_min=min(old_ratio),original_sd_ratio_max=max(old_ratio),
    candidate_sd_ratio_min=if(nrow(x))min(x$sd_ratio)else NA_real_,
    candidate_sd_ratio_max=if(nrow(x))max(x$sd_ratio)else NA_real_,
    original_max_abs_mean_error_sd=max(abs(as.numeric(A%*%lp$laplace$summary$fitness_mean)-s$mean)/s$sd),
    candidate_max_abs_mean_error_sd=if(nrow(x))max(abs(x$standardized_mean_difference))else NA_real_,
    any_material_mean=if(nrow(x))any(x$material_mean)else NA,
    any_material_sd=if(nrow(x))any(x$material_sd)else NA)
}))
write.csv(summary,file.path(output,"case_comparison_summary.csv"),row.names=FALSE)
print(summary,row.names=FALSE)
files <- c(file.path(directory,"analyze.R"),
  file.path(repo,"validation/posterior_reference_4to8/reference.R"),
  file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
write.csv(data.frame(file=files,sha256=unname(round2_sha256(files))),
  file.path(output,"analysis_dependencies_sha256.csv"),row.names=FALSE)
