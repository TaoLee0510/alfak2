design_nominal_initial <- function(context) {
  k<-context$graph$karyotypes;D<-as.matrix(dist(k,method="manhattan"));n<-nrow(k)
  center<-which(rowSums(abs(k-2L))==0);lo<-which(rowSums(abs(k-1L))==0)
  hi<-which(rowSums(abs(k-3L))==0)
  pA<-pB<-numeric(n);pA[center]<-.85;pA[D[center,]==1]<-.15/sum(D[center,]==1)
  pB[c(lo,hi)]<-.425
  near<-which(D[lo,]==1|D[hi,]==1);pB[near]<-.15/length(near)
  list(A=pA,B=pB,distance=D,edges=which(D==1&upper.tri(D),arr.ind=TRUE),
    planned_distance=apply(D[,which(pA>0),drop=FALSE],1,min))
}

design_simulate_path <- function(k,birth,beta,p0,N,seed,times=c(0,.75,1.5),max_events=1000000L) {
  stopifnot(all(diff(times)>0),times[1]==0,all(birth>0),sum(p0)>0)
  set.seed(seed);pop<-as.numeric(rmultinom(1,N,p0));initial<-pop
  snapshot<-matrix(0,length(pop),length(times));snapshot[,1]<-pop
  labels<-apply(k,1,paste,collapse=".")
  next_snapshot<-2L;t<-0;events<-0L
  repeat {
    rates<-pop*birth;total<-sum(rates)
    next_time<-if(total>0)t+rexp(1,total)else Inf
    while(next_snapshot<=length(times)&&next_time>times[next_snapshot]) {
      snapshot[,next_snapshot]<-pop;next_snapshot<-next_snapshot+1L
    }
    if(next_snapshot>length(times))break
    t<-next_time;parent<-sample.int(length(pop),1,prob=rates)
    x<-k[parent,];pop[parent]<-pop[parent]-1
    errors<-rbinom(1,sum(x),beta)
    if(!errors)pop[parent]<-pop[parent]+2 else {
      copies<-rep(seq_along(x),x);picked<-copies[sample.int(length(copies),errors)]
      signs<-sample(c(-1L,1L),errors,replace=TRUE)
      delta<-tabulate(picked[signs>0],length(x))-tabulate(picked[signs<0],length(x))
      for(daughter in list(x+delta,x-delta)) {
        j<-match(paste(daughter,collapse="."),labels)
        if(!is.na(j))pop[j]<-pop[j]+1
      }
    }
    events<-events+1L
    if(events>max_events)stop("Design process exceeded max_events; no seed replacement")
  }
  list(populations=snapshot,times=times,events=events,initial=initial,seed=seed)
}

design_count_uniform <- function(population,uniform,depth) {
  if(sum(population)<=0)stop("Design trajectory extinct: no count-only result available")
  cumulative<-pmin(1,cumsum(population/sum(population)));cumulative[length(cumulative)]<-1
  as.numeric(tabulate(findInterval(uniform[seq_len(depth)],cumulative)+1L,nbins=length(population)))
}

design_generate_unit <- function(context,family,landscape,replicate,truth_function) {
  fi<-match(family,c("smooth","rough","interaction"));seed<-5100000L+1000L*fi+landscape
  lattice<-list(karyotypes=context$graph$karyotypes)
  truth<-truth_function("design27_bounded",list(n_chr=3L,min_cn=1L,max_cn=3L),family,seed,lattice)
  nominal<-design_nominal_initial(context)
  ps<-5300000L+10000L*fi+100L*landscape+10L*replicate
  A1<-design_simulate_path(lattice$karyotypes,truth$birth,.005,nominal$A,600L,ps+1L)
  A2<-design_simulate_path(lattice$karyotypes,truth$birth,.005,nominal$A,600L,ps+2L)
  B<-design_simulate_path(lattice$karyotypes,truth$birth,.005,nominal$B,600L,ps+3L)
  pooled<-A1$populations+A2$populations
  ss<-5400000L+10000L*fi+100L*landscape+10L*replicate
  set.seed(ss);uA<-matrix(runif(300L*3L),300L,3L)
  set.seed(ss+1L);uB<-matrix(runif(150L*3L),150L,3L)
  make<-function(name,pop,times,traj,depths,uniforms) {
    y<-vapply(seq_along(times),function(o)design_count_uniform(pop[[o]],uniforms[[o]],depths[o]),numeric(context$n))
    actual<-apply(nominal$distance[,which(rowSums(y)>0),drop=FALSE],1,min)
    list(id=sprintf("design_%s_L%02d_R%d_%s",family,landscape,replicate,name),
      unit_id=sprintf("design_%s_L%02d_R%d",family,landscape,replicate),
      landscape_id=sprintf("design_%s_L%02d",family,landscape),family=family,
      landscape_replicate=landscape,biological_replicate=replicate,design=name,
      y=y,times=times,trajectory=as.integer(traj),reference_weight=nominal$A,
      truth=truth$birth,truth_relative=truth$birth-sum(nominal$A*truth$birth),
      actual_distance=actual,planned_distance=nominal$planned_distance,
      target_labels=context$graph$labels,target_edges=nominal$edges,
      landscape_seed=seed,process_seeds=ps+1:3,sampling_seeds=ss+0:1,
      total_typed_cells=sum(depths),total_initial_cells=1200L,
      cultures=max(traj)+1L,sampling_occasions=length(times),split="development")
  }
  counts<-list(
    control300=make("control300",list(pooled[,1],pooled[,3]),c(0,1.5),c(0,0),c(150L,150L),list(uA[,1],uA[,3])),
    depth600=make("depth600",list(pooled[,1],pooled[,3]),c(0,1.5),c(0,0),c(300L,300L),list(uA[,1],uA[,3])),
    time600=make("time600",list(pooled[,1],pooled[,2],pooled[,3]),c(0,.75,1.5),c(0,0,0),rep(200L,3),list(uA[,1],uA[,2],uA[,3])),
    composition600=make("composition600",list(A1$populations[,1],A1$populations[,3],B$populations[,1],B$populations[,3]),
      c(0,1.5,0,1.5),c(0,0,1,1),rep(150L,4),list(uA[,1],uA[,3],uB[,1],uB[,3])))
  baseline<-counts$control300$actual_distance
  for(name in names(counts))counts[[name]]$baseline_distance<-baseline
  list(cases=counts,latent_paths=list(A1=A1,A2=A2,B=B),truth=truth)
}
