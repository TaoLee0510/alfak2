#!/usr/bin/env Rscript
repo <- normalizePath(".")
directory <- file.path(repo,"validation/baseline_quadrature")
source(file.path(repo,"validation/posterior_reference_4to8/reference.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_posterior_reference.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
source(file.path(repo,"validation/accuracy_round2/score_models.R"))
source(file.path(directory,"quadrature.R"))
runtime <- reference_load(repo,compile=FALSE)
parity <- as.integer(Sys.getenv("ALFAK2_PILOT_PARITY","0"))
stopifnot(parity %in% 0:1)
output <- file.path(directory,"results",paste0("pilot_part",parity))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
input <- "/private/tmp/alfak2_iteration2_20261104_0518/process-data"
base <- "/private/tmp/alfak2_iteration2_20261104_0518/AB-pilot"
manifest <- read.csv(file.path(input,"dataset_manifest.csv"))
manifest <- manifest[manifest$mode=="finite_population",]
manifest <- manifest[seq_len(nrow(manifest))%%2==parity,]
source_files <- c(file.path(directory,c("quadrature.R","pilot.R","PILOT_PROTOCOL.md")),
  file.path(repo,"validation/posterior_reference_4to8/reference.R"),
  file.path(repo,"validation/accuracy_round2/score_models.R"),
  list.files(file.path(repo,"R"),"[.]R$",full.names=TRUE),runtime$native)
sha <- round2_sha256(source_files)
source_key <- round2_hash_object(sha)
write.csv(data.frame(file=names(sha),sha256=unname(sha)),file.path(output,"source_sha256.csv"),row.names=FALSE)
cache <- file.path("/private/tmp/alfak2_baseline_quad_pilot",substr(source_key,1,12))
dir.create(cache,recursive=TRUE,showWarnings=FALSE)
draws_from_components <- function(components,weight,n_draws,seed) {
  set.seed(seed); u <- runif(n_draws)
  index <- findInterval(u,cumsum(weight))+1L
  index <- pmin(length(weight),index)
  n <- length(components[[1]]$latent_mean)
  z <- matrix(rnorm(n*n_draws),nrow=n)
  draws <- matrix(NA_real_,n,n_draws)
  for(k in unique(index)) {
    ii <- which(index==k); c <- components[[k]]
    e <- eigen((c$latent_covariance+t(c$latent_covariance))/2,symmetric=TRUE)
    if(min(e$values)< -1e-7) stop("Indefinite joint component")
    L <- sweep(e$vectors,2,sqrt(pmax(0,e$values)),"*")
    draws[,ii] <- exp(sweep(L%*%z[,ii,drop=FALSE],1,c$latent_mean,"+"))
  }
  draws
}
metrics <- statuses <- input_hashes <- list()
for(i in seq_len(nrow(manifest))) {
  row <- manifest[i,]; dataset <- readRDS(row$dataset_file)
  old_file <- file.path(base,paste0("v1_joint_",dataset$dataset_id,".rds"))
  frozen <- readRDS(old_file); old <- frozen$fit$local
  input_hashes[[dataset$dataset_id]] <- data.frame(dataset_id=dataset$dataset_id,
    dataset_file=row$dataset_file,dataset_sha256=unname(round2_sha256(row$dataset_file)),
    v1_file=old_file,v1_sha256=unname(round2_sha256(old_file)))
  if(!is.null(dataset$native_target_distance)) dataset$target_distance <- dataset$native_target_distance
  graph <- old$graph; n <- length(graph$labels)
  prior <- runtime$pkg$build_joint_landscape_prior(graph,baseline_sd=.5,main_sd=.07,
    curvature_sd=.04,interaction_sd=.04,residual_sd=.125,edge_sd=.25,interaction_support="isotropic")
  K <- solve(as.matrix(prior$precision))+tcrossprod(sweep(prior$design,2,prior$coef_sd,"*"))
  case <- list(id=dataset$dataset_id,n=n,beta=dataset$config$beta,growth=FALSE,
    specification=list(dt=dataset$config$dt,growth=0),counts=dataset$counts,graph=graph,
    prior=prior,covariance=K,weight=dataset$gauge_weight[match(graph$labels,dataset$target_labels)])
  original_fit <- runtime$pkg$fit_local_posterior
  # This bridge changes only fixture defaults to the explicitly fixed deployment
  # defaults; every conditional prior vector is checked against frozen V1.
  bridge <- list(prepare_alfak2_data=runtime$pkg$prepare_alfak2_data,
    fit_local_posterior=function(...) {
      a <- list(...); a$eta_prior_sd <- 5; a$eta_borrowed_prior_mean <- -6
      a$eta_borrowed_prior_sd <- 1.5; a$eta_distance_penalty <- .75
      fit <- do.call(original_fit,a)
      stopifnot(max(abs(fit$summary$eta_prior_mean-old$summary$eta_prior_mean))<1e-12,
        max(abs(fit$summary$eta_prior_sd-old$summary$eta_prior_sd))<1e-12)
      fit
    })
  deployment <- list(pkg=bridge)
  replay <- reference_fit(case,deployment)
  stopifnot(max(abs(replay$latent_fitness_mean-old$latent_fitness_mean))<1e-6,
    max(abs(replay$latent_fitness_covariance-old$latent_fitness_covariance))<1e-6)
  case_key <- list(source=source_key,case=round2_hash_object(case),input=input_hashes[[dataset$dataset_id]])
  result_file <- file.path(output,paste0(dataset$dataset_id,"_integration.rds"))
  started <- proc.time()[["elapsed"]]
  if(file.exists(result_file)) {
    q <- readRDS(result_file)
    if(!identical(q$cache_key,case_key))stop("Pilot cache fingerprint mismatch")
  } else {
    q <- tryCatch(baseline_integrate(case,old,deployment,cache_directory=cache,cache_key=case_key),
      error=function(e)list(status="integration_error",message=conditionMessage(e)))
    q$cache_key <- case_key; q$elapsed_seconds <- proc.time()[["elapsed"]]-started
    saveRDS(q,result_file)
  }
  statuses[[dataset$dataset_id]] <- data.frame(dataset_id=dataset$dataset_id,
    family=dataset$family,landscape_id=dataset$landscape_id,status=q$status,
    final_order=if(!is.null(q$history)&&nrow(q$history))tail(q$history$order,1)else NA,
    elapsed_seconds=q$elapsed_seconds,v1_replay_passed=TRUE,
    message=if(is.null(q$message))""else q$message)
  order <- match(dataset$target_labels,graph$labels)
  for(model in c("v1_joint","baseline_quadrature")) {
    valid <- model=="v1_joint" || identical(q$status,"passed_quadrature_refinement")
    if(valid) {
      if(model=="v1_joint") {
        m <- old$summary$fitness_mean; C <- old$fitness_covariance
        components <- list(list(latent_mean=old$latent_fitness_mean,latent_covariance=old$latent_fitness_covariance));weight <- 1
      } else {m<-q$mean;C<-q$covariance;components<-q$components;weight<-q$weights}
      m<-m[order];C<-C[order,order,drop=FALSE]
      w<-dataset$gauge_weight;A<-diag(n)-matrix(w,n,n,byrow=TRUE)
      draws<-draws_from_components(components,weight,40000L,dataset$sampling_seed+71000000L)[order,,drop=FALSE]
      draws<-sweep(draws,2,as.numeric(crossprod(w,draws)),"-")
      nodes<-data.frame(dataset_id=dataset$dataset_id,landscape_id=dataset$landscape_id,
        family=dataset$family,split="development",model=model,karyotype=dataset$target_labels,
        distance=dataset$target_distance,truth_raw=dataset$truth,truth_relative=dataset$truth_relative,
        pred_raw=m,pred_relative=as.numeric(A%*%m),sd_raw=sqrt(diag(C)),
        sd_relative=sqrt(pmax(0,diag(A%*%C%*%t(A)))))
      for(level in dataset$config$interval_levels) {
        suffix<-sprintf("%02d",round(100*level));qs<-t(apply(draws,1,quantile,c((1-level)/2,(1+level)/2)))
        nodes[[paste0("relative_low_",suffix)]]<-qs[,1]
        nodes[[paste0("relative_high_",suffix)]]<-qs[,2]
      }
    } else {
      nodes<-frozen$nodes;nodes$model<-model
      for(column in c("pred_raw","pred_relative","sd_raw","sd_relative",paste0("relative_low_",c("50","80","95")),paste0("relative_high_",c("50","80","95"))))nodes[[column]]<-NA_real_
    }
    saveRDS(nodes,file.path(output,paste0(model,"_",dataset$dataset_id,"_nodes.rds")))
    metrics[[paste(model,dataset$dataset_id)]]<-round2_score_nodes(nodes,dataset)
  }
  saveRDS(list(dataset=dataset,frozen_v1=frozen$fit),file.path(output,paste0(dataset$dataset_id,"_inputs.rds")))
  write.csv(do.call(rbind,metrics),file.path(output,"case_metrics.csv"),row.names=FALSE)
  write.csv(do.call(rbind,statuses),file.path(output,"status.csv"),row.names=FALSE)
  write.csv(do.call(rbind,input_hashes),file.path(output,"input_sha256.csv"),row.names=FALSE)
  cat(dataset$dataset_id,q$status,"order",statuses[[dataset$dataset_id]]$final_order,"seconds",q$elapsed_seconds,"\n")
}
