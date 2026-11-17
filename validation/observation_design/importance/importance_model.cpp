#include <TMB.hpp>
#include "positive_taylor_forward.h"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(y0);DATA_VECTOR(y1);DATA_SCALAR(dt);DATA_INTEGER(steps);
  DATA_INTEGER(beta_zero);DATA_VECTOR(offset);DATA_MATRIX(lower);
  DATA_MATRIX(eta_basis);DATA_VECTOR(known_pi0);
  DATA_IVECTOR(trans_from);DATA_IVECTOR(trans_to);DATA_VECTOR(trans_weight);
  PARAMETER_VECTOR(z_g);PARAMETER_VECTOR(z_eta);
  int n=z_g.size();Type nll=0;
  vector<Type> g=offset+lower*z_g,b(n),p0(n),p1(n);
  for(int i=0;i<n;i++)b(i)=exp(g(i));
  if(known_pi0.size()==n) p0=known_pi0;
  else {
    vector<Type> eta=eta_basis*z_eta;Type mx=eta.maxCoeff(),total=0;
    for(int i=0;i<n;i++){p0(i)=exp(eta(i)-mx);total+=p0(i);}
    p0/=total;
  }
  Type growth=0;
  if(beta_zero) {
    int first=0;while(first<n&&asDouble(p0(first))==0)first++;
    if(first==n)error("Empty starting composition");
    Type mx=dt*b(first);
    for(int i=first+1;i<n;i++)if(asDouble(p0(i))>0)mx=CppAD::CondExpGt(dt*b(i),mx,dt*b(i),mx);
    Type total=0;
    for(int i=0;i<n;i++) {
      if(known_pi0.size()==n&&asDouble(known_pi0(i))==0)p1(i)=0;
      else p1(i)=p0(i)*exp(dt*b(i)-mx);
      total+=p1(i);
    }
    p1/=total;growth=mx+log(total);
  } else {
    DivisionForwardResult<Type> action=division_forward(p0,b,trans_from,trans_to,trans_weight,dt,steps,8);
    p1=action.frequency;growth=action.log_total;
  }
  for(int i=0;i<n;i++)nll-=y0(i)*log(p0(i)+Type(1e-14))+y1(i)*log(p1(i)+Type(1e-14));
  for(int i=0;i<z_g.size();i++)nll-=dnorm(z_g(i),Type(0),Type(1),true);
  for(int i=0;i<z_eta.size();i++)nll-=dnorm(z_eta(i),Type(0),Type(1),true);
  REPORT(g);REPORT(b);REPORT(p0);REPORT(p1);REPORT(growth);
  return nll;
}
