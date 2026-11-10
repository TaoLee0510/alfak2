#!/usr/bin/env Rscript
repo<-normalizePath(Sys.getenv("ALFAK2_REPO","."))
src<-Sys.getenv("ALFAK2_DESIGN_SOURCE",file.path(repo,"validation/observation_design"))
source(file.path(src,"model.R"));source(file.path(src,"data.R"))
source(file.path(repo,"validation/accuracy_round2/oracle_profile_fisher.R"))
runtime<-design_load(repo)
output<-Sys.getenv("ALFAK2_DESIGN_OUTPUT",file.path(repo,"validation/observation_design/results"))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
context<-design_context(as.matrix(expand.grid(rep(list(1:3),3))),.005,runtime)
nominal<-design_nominal_initial(context);n<-context$n
P<-matrix(0,n,n);op<-context$operator;P[cbind(op$from0+1L,op$to0+1L)]<-op$weight
x<-context$graph$karyotypes-2
profiles<-list(prior=context$prior$offset,
  smooth=.08*x[,1]-.06*x[,2]+.03*x[,3]-.04*rowSums(x^2),
  interaction=.16*x[,1]*x[,2]-.12*x[,2]*x[,3]+.12*x[,1]*x[,2]*x[,3])
designs<-list(control300=list(times=c(0,1.5),trajectory=c(0,0),depth=c(150,150)),
  depth600=list(times=c(0,1.5),trajectory=c(0,0),depth=c(300,300)),
  time600=list(times=c(0,.75,1.5),trajectory=c(0,0,0),depth=rep(200,3)),
  composition600=list(times=c(0,1.5,0,1.5),trajectory=c(0,0,1,1),depth=rep(150,4)))
anchor<-nominal$B
anchor[which(rowSums(abs(context$graph$karyotypes-2))==0)]<-.15
anchor[which(rowSums(abs(context$graph$karyotypes-1))==0|rowSums(abs(context$graph$karyotypes-3))==0)]<-.35
designs$composition_anchor600<-list(times=c(0,1.5,0,1.5),trajectory=c(0,0,1,1),depth=rep(150,4))
initial<-lapply(list(nominal$A,nominal$B,anchor),function(p){p<-pmax(p,1e-4);p/sum(p)})
A<-diag(n)-matrix(nominal$A,n,n,byrow=TRUE)
summary<-risks<-list();matrices<-list()
for(profile in names(profiles))for(scale in c(.8,1,1.2)) {
  g<-profiles[[profile]]+log(scale);b<-exp(g)
  jacobians<-lapply(initial,function(p)lapply(c(0,.75,1.5),function(t)
    oracle_fisher_jacobian(p,g,P,t,epsilon=5e-6)))
  point<-paste(profile,scale,sep="_")
  for(name in names(designs)) {
    d<-designs[[name]];nr<-max(d$trajectory)+1L;m<-length(d$times)
    Eg<-matrix(0,n*m,n);Ee<-matrix(0,n*m,(n-1L)*nr)
    for(o in seq_len(m)) {
      r<-d$trajectory[o]+1L;t<-match(d$times[o],c(0,.75,1.5))
      initial_index<-if(name=="composition_anchor600"&&r==2L)3L else r
      j<-jacobians[[initial_index]][[t]];rows<-(o-1L)*n+seq_len(n)
      eta<-if(d$times[o]==0)j$j0_eta else j$j1_eta
      p<-if(d$times[o]==0)j$pi0 else j$pi1
      w<-sqrt(d$depth[o]/p)
      Eg[rows,]<-j$j1_birth*w
      Ee[rows,(r-1L)*(n-1L)+seq_len(n-1L)]<-(eta%*%context$eta_basis)*w
    }
    s<-svd(Ee,nu=min(dim(Ee)),nv=0)
    rank<-sum(s$d>max(dim(Ee))*max(1,s$d)*1e-10)
    residual<-Eg
    if(rank) {Q<-s$u[,seq_len(rank),drop=FALSE];residual<-Eg-Q%*%crossprod(Q,Eg)}
    F<-crossprod(residual)
    white<-residual%*%t(chol(context$covariance))
    values<-pmax(0,eigen(crossprod(white),symmetric=TRUE,only.values=TRUE)$values)
    Cg<-solve(solve(context$covariance)+F)
    Cb<-diag(b)%*%Cg%*%diag(b);Cr<-A%*%Cb%*%t(A)
    matrices[[paste(point,name)]]<-F
    summary[[paste(point,name)]]<-data.frame(point=point,profile=profile,rate_scale=scale,
      design=name,typed_cells=sum(d$depth),eta_nuisance_rank=rank,
      effective_directions=sum(values/(1+values)),
      numerical_fisher_rank=sum(values>max(1,max(values))*n*1e-9),
      relative_mean_variance=mean(diag(Cr)),
      shape_mean_variance=mean(diag((diag(n)-matrix(1/n,n,n))%*%Cb%*%(diag(n)-matrix(1/n,n,n)))),
      interpretation="local expected multinomial information; eta prior omitted in profiling; not accuracy evidence")
    for(distance in 0:2)risks[[paste(point,name,distance)]]<-data.frame(point=point,
      design=name,original_planned_distance=distance,n_nodes=sum(nominal$planned_distance==distance),
      local_relative_variance=mean(diag(Cr)[nominal$planned_distance==distance]))
  }
  stopifnot(norm(matrices[[paste(point,"depth600")]]-2*matrices[[paste(point,"control300")]],"F")<1e-6)
}
write.csv(do.call(rbind,summary),file.path(output,"fisher_planning.csv"),row.names=FALSE)
write.csv(do.call(rbind,risks),file.path(output,"fisher_planned_distance_variance.csv"),row.names=FALSE)
saveRDS(list(initial=initial,profiles=profiles,scales=c(.8,1,1.2),fisher=matrices,
  initial_smoothing_floor=1e-4,scope="planning_only_no_RMSE_or_coverage_claim"),file.path(output,"fisher_details.rds"))
cat("Fisher planning across nine parameter points complete; depth scaling check passed\n")
