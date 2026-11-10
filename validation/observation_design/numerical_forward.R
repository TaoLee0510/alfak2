# Numerical-only additions to the frozen model helpers.
design_tmb_data_v4<-design_tmb_data
design_tmb_data<-function(case,context,shift=NULL,max_step=.1) {
  data<-design_tmb_data_v4(case,context,shift,max_step)
  data$taylor_order<-if(!is.null(shift)&&shift>5.5)64L else if(!is.null(shift)&&shift>3.5)32L else 8L
  data
}

direct_forward_r<-function(pi0,birth,operator,dt,steps,taylor_order) {
  n<-length(pi0);P<-matrix(0,n,n)
  P[cbind(operator$from0+1L,operator$to0+1L)]<-operator$weight
  A<-sweep(2*t(P)-diag(n),2,birth,"*");shift<-min(diag(A));B<-A-diag(shift,n)
  state<-as.numeric(pi0);growth<-0;h<-dt/steps
  for(s in seq_len(steps)) {
    term<-total<-state
    for(order in seq_len(taylor_order)){term<-as.numeric(B%*%term)*h/order;total<-total+term}
    state<-total;growth<-growth+h*shift
    z<-sum(state);growth<-growth+log(z);state<-state/z
  }
  list(pi1=state,log_total_growth=growth,generator=A)
}

direct_exact_forward<-function(generator,initial,dt,runtime) {
  n<-nrow(generator);shift<-max(colSums(generator))
  B<-generator-diag(shift,n)
  blocks<-max(1L,ceiling(dt*norm(B,"1")/20));h<-dt/blocks
  state<-as.numeric(initial);growth<-0
  for(block in seq_len(blocks)) {
    mass<-runtime$independent$process_deterministic_mean(B,state,h)
    total<-sum(mass)
    if(!all(is.finite(mass))||!is.finite(total)||total<=0)stop("Nonfinite independently normalized mean action")
    growth<-growth+log(total)+shift*h;state<-mass/total
  }
  list(pi1=state,log_total_growth=growth,blocks=blocks)
}
