constrained_kernel<-function(P,g,t=1.5,epsilon=5e-6) {
  M<-oracle_fisher_propagator(P,g,t)
  D<-lapply(seq_along(g),function(j) {
    d<-numeric(length(g));d[j]<-epsilon
    (oracle_fisher_propagator(P,g+d,t)-oracle_fisher_propagator(P,g-d,t))/(2*epsilon)
  })
  list(M=M,D=D,g=g)
}

constrained_information<-function(p,kernel,basis,depth,eta_prior=TRUE) {
  n<-length(p);m<-as.numeric(kernel$M%*%p);total<-sum(m);q<-m/total
  de<-diag(p)-tcrossprod(p);de<-de%*%basis
  de1<-kernel$M%*%de/total
  de1<-de1-tcrossprod(q,colSums(de1))
  dg<-vapply(kernel$D,function(D)as.numeric(D%*%p)/total,numeric(n))
  dg<-dg-tcrossprod(q,colSums(dg))
  E0<-de*sqrt(depth[1]/p);E1<-de1*sqrt(depth[2]/q)
  G<-dg*sqrt(depth[2]/q)
  if(eta_prior) {
    H<-diag(ncol(basis))+crossprod(E0)+crossprod(E1)
    cross<-crossprod(E1,G)
    F<-crossprod(G)-crossprod(cross,solve(H,cross))
  } else {
    E<-rbind(E0,E1);fullG<-rbind(matrix(0,n,n),G)
    s<-svd(E,nu=min(dim(E)),nv=0)
    rank<-sum(s$d>max(dim(E))*max(1,s$d)*1e-10)
    if(rank) {Q<-s$u[,seq_len(rank),drop=FALSE];fullG<-fullG-Q%*%crossprod(Q,fullG)}
    F<-crossprod(fullG)
  }
  (F+t(F))/2
}

constrained_risk<-function(initial,depths,kernel,context,A,eta_prior=TRUE) {
  F<-matrix(0,context$n,context$n)
  for(j in seq_along(initial))F<-F+constrained_information(initial[[j]],kernel,context$eta_basis,depths[[j]],eta_prior)
  Cg<-chol2inv(chol(context$landscape_precision+F))
  b<-exp(kernel$g);Cb<-Cg*tcrossprod(b)
  list(variance=diag(A%*%Cb%*%t(A)),covariance=Cb,fisher=F)
}

constrained_grid<-function(context,nominal) {
  shell<-rowSums(abs(context$graph$karyotypes-2));count<-tabulate(shell+1L,4)
  rows<-list();spec<-list();i<-0L
  for(a in 0:20)for(b in 0:(20-a))for(c in 0:(20-a-b)) {
    mass<-c(a,b,c,20-a-b-c)/20
    if(sum(mass[3:4])<.05-1e-12)next
    p<-mass[shell+1L]/count[shell+1L]
    for(frac in c(.4,.5,.6)) {
      i<-i+1L;id<-sprintf("single_%02d_%02d_%02d_%02d_t%02d",a,b,c,20-a-b-c,round(100*frac))
      spec[[id]]<-list(id=id,class="single",initial=list(p),depths=list(as.integer(600*c(frac,1-frac))),
        initial_cells=1200L,shell_mass=mass,anchor_mass=NA_real_)
      rows[[i]]<-data.frame(id=id,class="single",center=mass[1],shell1=mass[2],shell2=mass[3],shell3=mass[4],
        anchor=NA_real_,budgetA=600L,budgetB=0L,fractionA=frac,fractionB=NA_real_)
    }
  }
  for(budget in c(400L,450L))for(h in c(.15,.30,.45))for(fA in c(.4,.5,.6))for(fB in c(.4,.5,.6)) {
    p<-nominal$B;p[shell==0]<-h
    corners<-rowSums(abs(context$graph$karyotypes-1))==0|rowSums(abs(context$graph$karyotypes-3))==0
    p[corners]<-(.85-h)/2
    id<-sprintf("dual_A%d_h%02d_t%02d_%02d",budget,round(100*h),round(100*fA),round(100*fB));i<-i+1L
    spec[[id]]<-list(id=id,class="dual",initial=list(nominal$A,p),
      depths=list(as.integer(round(budget*c(fA,1-fA))),as.integer(round((600-budget)*c(fB,1-fB)))),
      initial_cells=c(600L,600L),shell_mass=NULL,anchor_mass=h)
    rows[[i]]<-data.frame(id=id,class="dual",center=NA_real_,shell1=NA_real_,shell2=NA_real_,shell3=NA_real_,
      anchor=h,budgetA=budget,budgetB=600L-budget,fractionA=fA,fractionB=fB)
  }
  list(spec=spec,table=do.call(rbind,rows))
}
