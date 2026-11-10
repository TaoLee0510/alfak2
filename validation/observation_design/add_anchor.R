#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
src<-Sys.getenv("ALFAK2_DESIGN_SOURCE",file.path(repo,"validation/observation_design"))
source(file.path(src,"data.R"));source(file.path(repo,"validation/accuracy_round2/common.R"))
out<-file.path(repo,"validation/observation_design/results")
context<-readRDS(file.path(out,"context.rds"));nominal<-design_nominal_initial(context)
path<-file.path(out,"dataset_manifest.csv");manifest<-read.csv(path)
original<-file.path(out,"dataset_manifest_initial48.csv")
if(!file.exists(original))file.copy(path,original)
manifest<-manifest[manifest$design!="composition_anchor600",]
p<-nominal$B;k<-context$graph$karyotypes
center<-which(rowSums(abs(k-2))==0);lo<-which(rowSums(abs(k-1))==0);hi<-which(rowSums(abs(k-3))==0)
p[center]<-.15;p[c(lo,hi)]<-.35
stopifnot(abs(sum(p)-1)<1e-12)
contrast_rank<-function(supports) {
  rows<-lapply(supports,function(s){M<-matrix(0,length(s)-1L,context$n)
    M[cbind(seq_len(nrow(M)),s[-1])]<-1;M[,s[1]]<- -1;M})
  qr(do.call(rbind,rows))$rank
}
rank_old<-contrast_rank(list(which(nominal$A>0),which(nominal$B>0)))
rank_new<-contrast_rank(list(which(nominal$A>0),which(p>0)))
stopifnot(rank_old==13L,rank_new==14L)
write.csv(data.frame(arm=c("composition600","composition_anchor600"),
  union_support_nodes=15L,beta_zero_contrast_rank=c(rank_old,rank_new),
  contrast_components=c(2L,1L)),file.path(out,"support_connectivity.csv"),row.names=FALSE)
rows<-list()
for(unit_id in unique(manifest$unit_id)) {
  unit<-readRDS(file.path(out,"data",paste0(unit_id,"_latent.rds")))
  old<-unit$cases$composition600
  seed<-old$process_seeds[1]+3L
  B<-design_simulate_path(k,unit$truth$birth,.005,p,600L,seed)
  set.seed(old$sampling_seeds[2]);uB<-matrix(runif(150L*3L),150L,3L)
  case<-old;case$id<-sub("composition600$","composition_anchor600",old$id)
  case$design<-"composition_anchor600"
  case$y[,3]<-design_count_uniform(B$populations[,1],uB[,1],150L)
  case$y[,4]<-design_count_uniform(B$populations[,3],uB[,3],150L)
  case$actual_distance<-apply(nominal$distance[,which(rowSums(case$y)>0),drop=FALSE],1,min)
  case$process_seeds<-c(old$process_seeds,seed)
  filename<-file.path(out,"data",paste0(case$id,".rds"))
  if(file.exists(filename)&&!identical(readRDS(filename),case))stop("Anchored design data changed")
  if(!file.exists(filename))saveRDS(case,filename)
  latent_file<-file.path(out,"data",paste0(unit_id,"_anchor_latent.rds"))
  if(!file.exists(latent_file))saveRDS(B,latent_file)
  rows[[case$id]]<-data.frame(case_id=case$id,unit_id=case$unit_id,landscape_id=case$landscape_id,
    family=case$family,biological_replicate=case$biological_replicate,design=case$design,
    typed_cells=600L,initial_cells=1200L,cultures=2L,observations=4L,
    actual_max_distance=max(case$actual_distance),baseline_max_distance=max(case$baseline_distance),
    file=normalizePath(filename),sha256=unname(round2_sha256(filename)))
}
write.csv(rbind(manifest,do.call(rbind,rows)),path,row.names=FALSE)
cat("Shared-anchor arm added; all original 48 datasets retained, 60 development cases total\n")
