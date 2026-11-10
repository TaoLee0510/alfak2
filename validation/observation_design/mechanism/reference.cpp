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
      // Independent full-generator action; no production/TMB forward header.
      std::vector<double> diagonal(n),B(n*n),term(n),acc(n),next(n);
      double shift=R_PosInf,norm=0;
      for(int i=0;i<n;i++) {diagonal[i]=(2*probability(i,i)-1)*b[i];shift=std::min(shift,diagonal[i]);}
      for(int from=0;from<n;from++) {
        double column=0;
        for(int to=0;to<n;to++) {
          double a=2*probability(from,to)*b[from];
          if(from==to) a-=b[from]+shift;
          if(a< -1e-12) stop("Independent Metzler shift failed");
          B[to*n+from]=a;column+=a;
        }
        norm=std::max(norm,column);
      }
      if(!std::isfinite(dt*norm)||dt*norm>1000000)stop("Independent mean-action numerical capacity exceeded; reference uncertified");
      int blocks=std::max(1,(int)std::ceil(dt*norm));
      double h=dt/blocks;pi1=pi0;
      for(int s=0;s<blocks;s++) {
        term=acc=pi1;bool converged=false;
        for(int order=1;order<=64;order++) {
          std::fill(next.begin(),next.end(),0);
          double maximum=0,scale=0;
          for(int to=0;to<n;to++) {
            for(int from=0;from<n;from++)next[to]+=B[to*n+from]*term[from]*h/order;
            acc[to]+=next[to];maximum=std::max(maximum,std::abs(next[to]));scale=std::max(scale,std::abs(acc[to]));
          }
          term=next;
          if(maximum<=1e-15*std::max(1.0,scale)){converged=true;break;}
        }
        if(!converged)stop("Independent mean series did not converge");
        double total=0;for(int i=0;i<n;i++)total+=acc[i];
        if(!(total>0)||!std::isfinite(total))return R_NegInf;
        log_growth+=h*shift+std::log(total);
        for(int i=0;i<n;i++)pi1[i]=acc[i]/total;
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

// [[Rcpp::export]]
List ref_preconditioned_ess_chain_cpp(List specification, NumericVector initial,
    int warmup, int retained, int seed, NumericVector center,
    NumericMatrix proposal_lower, NumericMatrix proposal_precision) {
  RefTarget target(specification);
  if(initial.size()!=target.q || center.size()!=target.q || retained<1 || warmup<0)
    stop("Invalid corrected-reference dimensions");
  std::vector<double> state(initial.begin(),initial.end()),candidate(state),nu(target.q),direction(target.q);
  auto residual=[&](const std::vector<double>& x) {
    double out=target.evaluate(x);
    if(!std::isfinite(out)) return R_NegInf;
    double prior=0,reference=0;
    for(int j=0;j<target.q;j++) {
      prior+=x[j]*x[j];
      for(int k=0;k<target.q;k++)reference+=(x[j]-center[j])*proposal_precision(j,k)*(x[k]-center[k]);
    }
    return out-.5*prior+.5*reference;
  };
  double current=residual(state);
  if(!std::isfinite(current))stop("Nonfinite corrected-reference starting density");
  NumericMatrix draws(retained,target.q);long evaluations=1,rejected=0,nonfinite=0;
  const double two_pi=2*std::acos(-1.0);
  for(int iteration=0;iteration<warmup+retained;iteration++) {
    if(iteration%1000==0)checkUserInterrupt();
    for(int j=0;j<target.q;j++)nu[j]=R::rnorm(0,1);
    for(int j=0;j<target.q;j++) {
      direction[j]=0;
      for(int k=0;k<target.q;k++)direction[j]+=proposal_lower(j,k)*nu[k];
    }
    double threshold=current+std::log(R::runif(0,1));
    double angle=R::runif(0,two_pi),lo=angle-two_pi,hi=angle;int trials=0;
    while(true) {
      for(int j=0;j<target.q;j++)candidate[j]=center[j]+(state[j]-center[j])*std::cos(angle)+direction[j]*std::sin(angle);
      double value=residual(candidate);evaluations++;
      if(!std::isfinite(value))nonfinite++;
      if(value>threshold){state=candidate;current=value;break;}
      rejected++;trials++;
      if(trials>10000)stop("Corrected elliptical slice bracket did not close");
      if(angle<0)lo=angle;else hi=angle;
      angle=R::runif(lo,hi);
    }
    if(iteration>=warmup)for(int j=0;j<target.q;j++)draws(iteration-warmup,j)=state[j];
  }
  return List::create(_["draws"]=draws,_["last"]=wrap(state),
    _["likelihood_evaluations"]=(double)evaluations,_["rejected_angles"]=(double)rejected,
    _["nonfinite_proposals"]=(double)nonfinite,_["seed"]=seed);
}
