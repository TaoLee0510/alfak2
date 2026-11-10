#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."));root<-file.path(repo,"validation/observation_design")
out<-normalizePath(Sys.getenv("ALFAK2_IS_OUTPUT"));dest<-file.path(out,"mechanism")
stopifnot(!dir.exists(dest));dir.create(dest);dir.create(file.path(dest,"data"));dir.create(file.path(dest,"inference"));dir.create(file.path(dest,"fits"))
source(file.path(root,"results/computation_source_v4/data.R"))
source(file.path(repo,"validation/accuracy_round2/process_prepare_data.R"))
ctxfile<-file.path(root,"results/context.rds");ctx<-readRDS(ctxfile)
old<-read.csv(file.path(root,"results/constrained_v2/all_fit_status.csv"));old<-old[old$design=="opt_single",]
pre_manifest<-normalizePath(Sys.getenv("ALFAK2_IS_PREFLIGHT_MANIFEST",file.path(out,"preflight_manifest.csv")))
pre<-read.csv(pre_manifest);rows<-audit<-list()
for(i in seq_len(nrow(old))) {
 base<-readRDS(old$file[i]);pathsfile<-file.path(root,"results/constrained_v2/latent_paths",paste0(base$id,".rds"))
 path<-readRDS(pathsfile)[[1]];initial<-path$initial
 stopifnot(sum(initial)==1200,identical(as.numeric(path$populations[,1]),as.numeric(initial)))
 set.seed(base$sampling_seeds[1]);uniform<-matrix(runif(600L*4L),600L,4L)
 stopifnot(identical(as.numeric(design_count_uniform(initial,uniform[,1],240)),as.numeric(base$y[,1])),
  identical(as.numeric(design_count_uniform(path$populations[,2],uniform[,2],360)),as.numeric(base$y[,2])))
 independent<-process_independent_operator(ctx$graph$karyotypes,base$truth,.005)
 meanpop<-process_deterministic_mean(independent$generator,initial,1.5)
 mean_counts<-design_count_uniform(meanpop,uniform[,2],360)
 for(arm in c("A_finite_unknown","B_mean_unknown","C_finite_oracle","D_mean_oracle")) {
  c<-base;c$design<-arm;c$mechanism_arm<-arm
  if(arm!="A_finite_unknown")c$id<-paste0(base$unit_id,"_",arm)
  if(arm%in%c("B_mean_unknown","D_mean_oracle"))c$y[,2]<-mean_counts
  if(arm%in%c("C_finite_oracle","D_mean_oracle"))c$known_pi0<-as.numeric(initial/sum(initial))
  distance<-as.matrix(dist(ctx$graph$karyotypes,method="manhattan"))
  c$actual_distance<-apply(distance[,which(rowSums(c$y)>0),drop=FALSE],1,min)
  c$pairing_scope<-"same realized initial integer population, t0 observations, endpoint sampling uniforms and original-control distance/reference"
  f<-file.path(dest,"data",paste0(c$id,".rds"));saveRDS(c,f)
  inference<-c[intersect(c("id","y","times","trajectory","reference_weight","target_edges","known_pi0"),names(c))]
  inf<-file.path(dest,"inference",paste0(c$id,".rds"));saveRDS(inference,inf)
  output<-file.path(dest,"fits",c$id);seed<-10600000L+i*10000L+match(arm,c("A_finite_unknown","B_mean_unknown","C_finite_oracle","D_mean_oracle"))-1L
  reuse<-FALSE
  if(arm=="A_finite_unknown"&&c$id%in%pre$case_id) {
   p<-pre[pre$case_id==c$id,];stopifnot(identical(inference,readRDS(p$file)),seed==p$seed)
   fitted<-readRDS(file.path(p$output,"result.rds"));resource<-read.csv(file.path(p$output,"resource.csv"))
   stopifnot(all(fitted$flags),fitted$numerical_seed==seed,tolower(as.character(resource$resource_pass))=="true")
   ledger<-read.csv(file.path(p$output,"source_sha256.csv"))
   sources<-file.path(normalizePath(Sys.getenv("ALFAK2_IS_SOURCE")),c("model.R","integrate.R","importance.cpp","importance_model.cpp","positive_taylor_forward.h","run_case.R"))
   stopifnot(all(sources%in%ledger$file),all(vapply(sources,digest::digest,character(1),algo="sha256",file=TRUE)==ledger$sha256[match(sources,ledger$file)]))
   output<-p$output;reuse<-TRUE
  }
  rows[[length(rows)+1L]]<-data.frame(case_id=c$id,unit_id=c$unit_id,landscape_id=c$landscape_id,family=c$family,
   arm=arm,oracle=grepl("oracle",arm),file=inf,scoring_file=f,context=ctxfile,output=output,seed=seed,
   reuse_preflight=reuse,sha256=digest::digest(file=inf,algo="sha256"),scoring_sha256=digest::digest(file=f,algo="sha256"))
 }
 audit[[i]]<-data.frame(unit_id=base$unit_id,original_file=old$file[i],paths_file=pathsfile,
  initial_population=sum(initial),typed_t0=sum(base$y[,1]),typed_endpoint=sum(mean_counts),
  finite_uniform_replay=TRUE,known_pi0_zero_states=sum(initial==0),sampling_seed=base$sampling_seeds[1])
}
manifest<-do.call(rbind,rows);stopifnot(nrow(manifest)==48,!anyDuplicated(manifest$case_id),sum(manifest$reuse_preflight)==2)
write.csv(manifest,file.path(dest,"manifest.csv"),row.names=FALSE)
write.csv(manifest[!manifest$reuse_preflight,],file.path(dest,"run_manifest.csv"),row.names=FALSE)
write.csv(do.call(rbind,audit),file.path(dest,"pairing_audit.csv"),row.names=FALSE)
files<-unique(c(list.files(file.path(dest,"data"),full.names=TRUE),list.files(file.path(dest,"inference"),full.names=TRUE),
 ctxfile,pre_manifest,old$file,vapply(audit,`[[`,character(1),"paths_file"),normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])),
 file.path(root,"results/computation_source_v4/data.R"),file.path(repo,"validation/accuracy_round2/process_prepare_data.R")))
write.csv(data.frame(file=files,sha256=vapply(files,digest::digest,character(1),file=TRUE,algo="sha256")),file.path(dest,"generation_source_sha256.csv"),row.names=FALSE)
writeLines("48 fixed-layout development cases frozen before fitting; oracle arms are mechanism probes",file.path(dest,"DATA_FROZEN"))
