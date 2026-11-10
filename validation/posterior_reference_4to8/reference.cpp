#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Independent diagnostic implementation. This file includes no package/TMB
// forward code and samples the Gaussian-prior likelihood target, not a fitted
// Gaussian posterior approximation.
struct RefTarget {
  int n, q, steps;
  double dt, growth, growth_sd;
  NumericVector offset, y0, y1, known;
  NumericMatrix lower, eta_basis, probability;
  bool beta_zero;
  RefTarget(List a) : n(as<int>(a["n"])), q(as<int>(a["q"])),
    steps(as<int>(a["steps"])), dt(as<double>(a["dt"])),
    growth(as<double>(a["growth"])), growth_sd(as<double>(a["growth_sd"])),
    offset(a["offset"]), y0(a["y0"]), y1(a["y1"]), known(a["known_pi0"]),
    lower(a["lower"]), eta_basis(a["eta_basis"]), probability(a["probability"]),
    beta_zero(as<bool>(a["beta_zero"])) {}

  double evaluate(const std::vector<double>& z, std::vector<double>* pout=nullptr,
                  double* gout=nullptr) const {
    std::vector<double> b(n), pi0(n), pi1(n), eta(n,0);
    for(int i=0;i<n;i++) {
      double g=offset[i];
      for(int j=0;j<n;j++) g+=lower(i,j)*z[j];
      b[i]=std::exp(g);
      if(!std::isfinite(b[i])) return R_NegInf;
    }
    if(known.size()==n) for(int i=0;i<n;i++) pi0[i]=known[i];
    else {
      for(int i=0;i<n;i++) for(int j=0;j<q-n;j++) eta[i]+=eta_basis(i,j)*z[n+j];
      double mx=*std::max_element(eta.begin(),eta.end()), total=0;
      for(int i=0;i<n;i++) total+=(pi0[i]=std::exp(eta[i]-mx));
      for(int i=0;i<n;i++) pi0[i]/=total;
    }
    double log_growth=0;
    if(beta_zero) {
      double mx=R_NegInf;
      for(int i=0;i<n;i++) {
        pi1[i]=std::log(pi0[i])+dt*b[i];
        mx=std::max(mx,pi1[i]);
      }
      if(!std::isfinite(mx)) return R_NegInf;
      double total=0;
      for(int i=0;i<n;i++) total+=(pi1[i]=std::exp(pi1[i]-mx));
      log_growth=mx+std::log(total);
      for(int i=0;i<n;i++) pi1[i]/=total;
    } else {
      std::vector<double> d(n),half(n),term(n),acc(n),next(n);
      double shift=R_NegInf, h=dt/steps;
      for(int i=0;i<n;i++) {d[i]=(2*probability(i,i)-1)*b[i];shift=std::max(shift,d[i]);}
      for(int i=0;i<n;i++) {half[i]=std::exp(.5*h*(d[i]-shift));pi1[i]=pi0[i];}
      for(int s=0;s<steps;s++) {
        for(int i=0;i<n;i++) term[i]=acc[i]=(pi1[i]*half[i]);
        for(int order=1;order<=4;order++) {
          std::fill(next.begin(),next.end(),0);
          for(int from=0;from<n;from++) for(int to=0;to<n;to++) if(from!=to)
            next[to]+=2*probability(from,to)*b[from]*term[from]*h/order;
          for(int i=0;i<n;i++) {acc[i]+=next[i];term[i]=next[i];}
        }
        for(int i=0;i<n;i++) pi1[i]=acc[i]*half[i];
        log_growth+=h*shift;
        if((s+1)%16==0 || s+1==steps) {
          double total=0;for(int i=0;i<n;i++) total+=pi1[i];
          if(!(total>0) || !std::isfinite(total)) return R_NegInf;
          log_growth+=std::log(total);
          for(int i=0;i<n;i++) pi1[i]/=total;
        }
      }
    }
    double ll=0;
    for(int i=0;i<n;i++) ll+=y0[i]*std::log(pi0[i]+1e-14)+y1[i]*std::log(pi1[i]+1e-14);
    if(growth_sd>0) ll+=R::dnorm(log_growth,growth,growth_sd,true);
    if(pout) *pout=pi1;
    if(gout) *gout=log_growth;
    return std::isfinite(ll) ? ll : R_NegInf;
  }
};

// [[Rcpp::export]]
List ref_likelihood_cpp(NumericVector z, List specification) {
  RefTarget target(specification);
  if(z.size()!=target.q) stop("Reference coordinate length mismatch");
  std::vector<double> state(z.begin(),z.end()),p;
  double growth=NA_REAL;
  double ll=target.evaluate(state,&p,&growth);
  return List::create(_["log_likelihood"]=ll,_["pi1"]=wrap(p),_["log_growth"]=growth);
}

// [[Rcpp::export]]
List ref_ess_chain_cpp(List specification, NumericVector initial, int warmup,
                       int retained, int seed) {
  RefTarget target(specification);
  if(initial.size()!=target.q || retained<1 || warmup<0) stop("Invalid sampler dimensions");
  std::vector<double> state(initial.begin(),initial.end()),proposal(state),nu(target.q);
  double current=target.evaluate(state), angle=0;
  if(!std::isfinite(current)) stop("Nonfinite initial reference likelihood");
  NumericMatrix draws(retained,target.q);
  long evaluations=1, rejected=0, nonfinite=0;
  const double two_pi=2*std::acos(-1.0);
  for(int iteration=0;iteration<warmup+retained;iteration++) {
    if(iteration%1000==0) checkUserInterrupt();
    for(int block=0;block<(target.q>target.n?2:1);block++) {
      int first=block==0?0:target.n, last=block==0?target.n:target.q;
      for(int j=first;j<last;j++) nu[j]=R::rnorm(0,1);
      double threshold=current+std::log(R::runif(0,1));
      angle=R::runif(0,two_pi);
      double lo=angle-two_pi,hi=angle;
      int trials=0;
      while(true) {
        for(int j=first;j<last;j++) proposal[j]=state[j]*std::cos(angle)+nu[j]*std::sin(angle);
        double ll=target.evaluate(proposal);evaluations++;
        if(!std::isfinite(ll)) nonfinite++;
        if(ll>threshold) {
          for(int j=first;j<last;j++) state[j]=proposal[j];
          current=ll;break;
        }
        rejected++;trials++;
        if(trials>10000) stop("Elliptical slice bracket failed to close");
        if(angle<0) lo=angle;else hi=angle;
        angle=R::runif(lo,hi);
      }
      // Keep every coordinate outside the next block at its current value.
      proposal=state;
    }
    if(iteration>=warmup) for(int j=0;j<target.q;j++) draws(iteration-warmup,j)=state[j];
  }
  return List::create(_["draws"]=draws,_["last"]=wrap(state),
    _["likelihood_evaluations"]=(double)evaluations,_["rejected_angles"]=(double)rejected,
    _["nonfinite_proposals"]=(double)nonfinite,_["seed"]=seed);
}
