#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
src<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_SOURCE"));out<-normalizePath(Sys.getenv("ALFAK2_CONSTRAINED_OUTPUT"))
source(file.path(repo,"validation/accuracy_round2/process_prepare_data.R"))
source(file.path(repo,"validation/accuracy_round2/common.R"))
planning<-readRDS(file.path(out,"planning_context.rds"));context<-planning$context;n<-context$n
destination<-file.path(out,"noise");dir.create(destination,showWarnings=FALSE)
stopifnot(!file.exists(file.path(destination,"COMPLETE")))
initial<-round(1200*planning$nominal$A);stopifnot(sum(initial)==1200,max(abs(initial-1200*planning$nominal$A))<1e-10)
summaries<-nodes<-seeds<-list();idx<-0L
for(profile in names(planning$profiles)) {
  idx<-idx+1L;birth<-exp(planning$profiles[[profile]])
  independent<-process_independent_operator(context$graph$karyotypes,birth,.005)
  expected<-process_deterministic_mean(independent$generator,initial,1.5);qdet<-expected/sum(expected)
  cached<-planning$kernels[[paste0(profile,"_1")]]$M%*%initial
  stopifnot(max(abs(qdet-as.numeric(cached/sum(cached))))<1e-10)
  populations<-matrix(NA_real_,100,n);events<-integer(100)
  for(rep in 1:100) {
    seed<-8300000L+idx*1000L+rep
    sim<-process_exact_gillespie(context$graph$karyotypes,birth,.005,1.5,initial,seed)
    populations[rep,]<-sim$population;events[rep]<-sim$events
    seeds[[length(seeds)+1L]]<-data.frame(profile=profile,replicate=rep,seed=seed,events=sim$events)
    if(rep%%25==0) {cat("Noise",profile,rep,"/100\n");flush.console()}
  }
  Q<-populations/rowSums(populations);qm<-colMeans(Q);sampling<-matrix(0,n,n)
  for(rep in 1:100)sampling<-sampling+(diag(Q[rep,])-tcrossprod(Q[rep,]))/(300*100)
  process<-cov(Q);total<-sampling+process;deterministic<-(diag(qdet)-tcrossprod(qdet))/300
  stopifnot(max(abs(rowSums(total)))<1e-10,min(eigen(total,symmetric=TRUE,only.values=TRUE)$values)>-1e-10)
  bootstrap<-matrix(NA_real_,1000,3);set.seed(8390000L+idx)
  for(b in 1:1000) {
    q<-Q[sample.int(100,100,replace=TRUE),,drop=FALSE];meanq<-colMeans(q)
    vs<-(meanq-colMeans(q*q))/300;vp<-apply(q,2,var)
    for(d in 0:2) {
      ii<-planning$nominal$planned_distance==d
      bootstrap[b,d+1L]<-mean((vs+vp)[ii])/mean(diag(deterministic)[ii])
    }
  }
  for(d in 0:2) {
    ii<-planning$nominal$planned_distance==d;den<-mean(diag(deterministic)[ii]);ci<-quantile(bootstrap[,d+1L],c(.025,.975))
    summaries[[length(summaries)+1L]]<-data.frame(profile=profile,planned_distance=d,paths=100,
      deterministic_sampling_variance=den,stochastic_sampling_variance=mean(diag(sampling)[ii]),
      branching_process_variance=mean(diag(process)[ii]),total_variance=mean(diag(total)[ii]),
      variance_ratio=mean(diag(total)[ii])/den,bootstrap_lo=ci[1],bootstrap_hi=ci[2],
      process_fraction=mean(diag(process)[ii])/mean(diag(total)[ii]),
      scope="endpoint fraction variance at fixed initial population; not fitness posterior coverage")
  }
  popmean<-colMeans(populations);popse<-apply(populations,2,sd)/10
  nodes[[profile]]<-data.frame(profile=profile,node=context$graph$labels,planned_distance=planning$nominal$planned_distance,
    deterministic_fraction=qdet,mean_process_fraction=qm,process_fraction_mcse=apply(Q,2,sd)/10,
    mean_population=popmean,expected_population=expected,mean_population_mcse=popse,
    mean_population_z=ifelse(popse>0,(popmean-expected)/popse,NA_real_),
    sampling_variance=diag(sampling),process_variance=diag(process),total_variance=diag(total),
    deterministic_sampling_variance=diag(deterministic))
  saveRDS(list(profile=profile,initial=initial,populations=populations,Q=Q,expected=expected,
    qdet=qdet,sampling_covariance=sampling,process_covariance=process,total_covariance=total,
    deterministic_covariance=deterministic,bootstrap_variance_ratios=bootstrap),file.path(destination,paste0(profile,".rds")))
}
write.csv(do.call(rbind,summaries),file.path(destination,"variance_decomposition.csv"),row.names=FALSE)
write.csv(do.call(rbind,nodes),file.path(destination,"node_diagnostics.csv"),row.names=FALSE)
write.csv(do.call(rbind,seeds),file.path(destination,"seed_ledger.csv"),row.names=FALSE)
files<-c(file.path(src,"noise.R"),file.path(repo,"validation/accuracy_round2/process_prepare_data.R"),
  file.path(repo,"validation/accuracy_round2/common.R"),file.path(out,"planning_context.rds"))
write.csv(data.frame(file=files,sha256=unname(round2_sha256(files))),file.path(destination,"source_sha256.csv"),row.names=FALSE)
writeLines("300 fixed-initial-population trajectories complete; no posterior calibration or likelihood change.",file.path(destination,"COMPLETE"))
cat("Independent process-noise diagnostic complete\n")
